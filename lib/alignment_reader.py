"""
Class to read alignment files and create position dictionaries or interlap dictionaries
for read counting.
"""


from dataclasses import dataclass
from numbers import Integral
from pathlib import Path
from urllib.parse import quote

import pysam
import interlap


import lib.messaging as msg


LIBRARY_MINIMUM_KEY = "__library_minimum__"
WIG_FILENAME_SAFE = ".:^$@!+_?-|"


@dataclass(frozen=True, slots=True)
class AlignmentPolicy:
    """Filtering applied before an alignment record contributes to any count.

    Counts always use alignment records as their unit. The defaults preserve
    ORFBounder's historical behavior: accept MAPQ 0, retain duplicate-marked
    records, and exclude NH-tagged multimappers.
    """

    min_mapq: int = 0
    include_duplicates: bool = True
    multimappers: str = "exclude"

    def __post_init__(self):
        if (isinstance(self.min_mapq, bool) or not isinstance(self.min_mapq, Integral)
                or not 0 <= self.min_mapq <= 255):
            raise ValueError("min_mapq must be an integer from 0 through 255.")
        if not isinstance(self.include_duplicates, bool):
            raise ValueError("include_duplicates must be true or false.")
        if self.multimappers not in {"exclude", "include"}:
            raise ValueError("multimappers must be 'exclude' or 'include'.")
        # Normalize integer-like values (for example numpy integers) while
        # retaining a frozen, predictable public value.
        object.__setattr__(self, "min_mapq", int(self.min_mapq))


DEFAULT_ALIGNMENT_POLICY = AlignmentPolicy()
_EXCLUSION_REASONS = (
    "unmapped", "secondary", "supplementary", "qc_fail",
    "below_min_mapq", "duplicate", "multimapper",
)


def _resolve_alignment_policy(policy):
    if policy is None:
        return DEFAULT_ALIGNMENT_POLICY
    if not isinstance(policy, AlignmentPolicy):
        raise TypeError("alignment_policy must be an AlignmentPolicy instance.")
    return policy


def _initialize_alignment_diagnostics(diagnostics, policy):
    """Initialize or validate an accumulator used by the streaming iterator."""
    if diagnostics is None:
        diagnostics = {}
    policy_record = {
        "min_mapq": policy.min_mapq,
        "include_duplicates": policy.include_duplicates,
        "multimappers": policy.multimappers,
    }
    if "policy" in diagnostics and diagnostics["policy"] != policy_record:
        raise ValueError("Cannot combine alignment diagnostics from different policies.")
    diagnostics.setdefault("policy", policy_record)
    diagnostics.setdefault("count_unit", "alignment_record")
    for field in (
        "total_records", "accepted_records", "excluded_records",
        "records_without_nh", "paired_records", "duplicate_records",
        "multimapped_records",
    ):
        diagnostics.setdefault(field, 0)
    exclusions = diagnostics.setdefault("exclusions", {})
    for reason in _EXCLUSION_REASONS:
        exclusions.setdefault(reason, 0)
    return diagnostics


def _exclude_alignment(diagnostics, reason):
    diagnostics["excluded_records"] += 1
    diagnostics["exclusions"][reason] += 1


def validate_normalization_scope(normalization_scope, normalization=None):
    """Check the shared scope for positional scaling and expression denominators."""
    if normalization_scope not in {"contig", "library"}:
        raise ValueError("normalization_scope must be contig or library.")


def alignment_query_length(read):
    """Return query length, inferring it from CIGAR when SAM omits SEQ."""
    query_length = read.query_length
    if not query_length:
        query_length = read.infer_query_length()
    return query_length


def iter_primary_alignments(alignment_file_path, alignment_policy=None, diagnostics=None):
    """Stream policy-accepted alignment records without requiring an index.

    NH is optional in SAM. With the default policy, NH-tagged multimappers are
    excluded while records without NH are retained. Duplicate-marked records and
    MAPQ 0 are also retained by default. Paired mates remain separate alignment
    records; this function does not infer or count fragments.

    If supplied, ``diagnostics`` is updated in place. Every rejected record is
    assigned one exclusion reason, in flag/filter order, so accepted plus
    excluded equals the total number of records seen.
    """
    policy = _resolve_alignment_policy(alignment_policy)
    diagnostics = _initialize_alignment_diagnostics(diagnostics, policy)
    try:
        alignment_file = pysam.AlignmentFile(alignment_file_path)
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read alignment file {alignment_file_path}: {exc}") from exc
    try:
        for read in alignment_file.fetch(until_eof=True):
            diagnostics["total_records"] += 1
            has_nh = read.has_tag("NH")
            if not has_nh:
                diagnostics["records_without_nh"] += 1
            if read.is_paired:
                diagnostics["paired_records"] += 1
            if read.is_duplicate:
                diagnostics["duplicate_records"] += 1
            if has_nh:
                nh = read.get_tag("NH")
                if isinstance(nh, bool) or not isinstance(nh, Integral) or nh < 1:
                    raise ValueError(
                        f"Alignment {read.query_name!r} has an invalid NH tag {nh!r}; "
                        "expected a positive integer."
                    )
                is_multimapped = nh > 1
            else:
                is_multimapped = False
            if is_multimapped:
                diagnostics["multimapped_records"] += 1

            reason = None
            if read.is_unmapped:
                reason = "unmapped"
            elif read.is_secondary:
                reason = "secondary"
            elif read.is_supplementary:
                reason = "supplementary"
            elif read.is_qcfail:
                reason = "qc_fail"
            elif read.mapping_quality < policy.min_mapq:
                reason = "below_min_mapq"
            elif not policy.include_duplicates and read.is_duplicate:
                reason = "duplicate"
            elif policy.multimappers == "exclude" and is_multimapped:
                reason = "multimapper"

            if reason is not None:
                _exclude_alignment(diagnostics, reason)
                continue
            if read.cigartuples is None:
                raise ValueError(
                    f"Alignment {read.query_name!r} is marked as mapped but has no CIGAR; "
                    "genomic coverage cannot be determined."
                )
            diagnostics["accepted_records"] += 1
            yield read
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read alignment file {alignment_file_path}: {exc}") from exc
    finally:
        alignment_file.close()


class PositionReader:
    """
    Read sam/bam file and create a position dict for requested read lengths
    read_lengths : list (e.g. [30,32,40])
    mapping_mode : threeprime, fiveprime, centered, global
    """
    def __init__(self, alignment_file_path, read_length_dict, mapping_mode, offset_dict,
                 alignment_policy=None, circular_contigs=None):
        self.alignment_file_path = Path(alignment_file_path)
        self.read_length_dict = read_length_dict
        self.mapping_mode = mapping_mode
        self.alignment_policy = _resolve_alignment_policy(alignment_policy)
        self.circular_contigs = set(circular_contigs or ())
        self.alignment_diagnostics = {}
        if mapping_mode not in {"fiveprime", "threeprime", "centered", "global"}:
            raise ValueError(f"Unsupported mapping mode: {mapping_mode}")

        self.wildcard = self.alignment_file_path.stem.split('_')[0]

        if not self.read_length_dict:
            self.read_lengths = None
        elif self.wildcard in self.read_length_dict:
            self.read_lengths = self.read_length_dict[self.wildcard]
        elif "default" in self.read_length_dict:
            self.read_lengths = self.read_length_dict["default"]
        else:
            msg.warning("Warning no default value given for read-lengths. Using all read lengths.")
            self.read_lengths = None

        if self.wildcard in offset_dict:
            self.offset_dict = offset_dict[self.wildcard]
        elif "default" in offset_dict:
            self.offset_dict = offset_dict["default"]
        else:
            raise ValueError(f"Error: Offsets and default value missing for wildcard: {self.wildcard}")

        self.reads_position_dict = {}
        self.no_accepted_reads_dict = {}

        msg.message(f"Reading read positions from alignment file: {self.alignment_file_path.name}")
        self._read_alignment_file()

    def _read_alignment_file(self):
        """Map aligned query bases to zero-based genomic coordinates.

        Offsets preserve the historical convention: subtract on +, add on -.
        Centered mode discards 11 query bases at each end and distributes one
        read across the remaining aligned bases. CIGAR gaps receive no signal.
        """
        for read in iter_primary_alignments(
                self.alignment_file_path, self.alignment_policy, self.alignment_diagnostics):
            chrom = read.reference_name
            read_length = alignment_query_length(read)
            if self.read_lengths is not None and str(read_length) not in self.read_lengths:
                continue
            if str(read_length) in self.offset_dict:
                offset = self.offset_dict[str(read_length)]
            elif read_length in self.offset_dict:
                offset = self.offset_dict[read_length]
            elif "default" in self.offset_dict:
                offset = self.offset_dict["default"]
            else:
                raise ValueError(f"Error: Offsets and default value missing for (wildcard, readlength): {self.wildcard}, {read_length}")

            aligned_positions = read.get_reference_positions()
            if not aligned_positions:
                continue
            strand = "-" if read.is_reverse else "+"
            shift = offset if read.is_reverse else -offset
            self.no_accepted_reads_dict[chrom] = self.no_accepted_reads_dict.get(chrom, 0) + 1
            if self.mapping_mode == "fiveprime":
                positions = [aligned_positions[-1] if read.is_reverse else aligned_positions[0]]
            elif self.mapping_mode == "threeprime":
                positions = [aligned_positions[0] if read.is_reverse else aligned_positions[-1]]
            elif self.mapping_mode == "centered":
                positions = [position for position in read.get_reference_positions(full_length=True)[11:-11]
                             if position is not None]
            else:
                positions = aligned_positions
            if not positions:
                continue

            weight = 1 / len(positions) if self.mapping_mode == "centered" else 1
            reference_length = read.header.get_reference_length(chrom)
            for position in positions:
                position += shift
                if chrom in self.circular_contigs:
                    position %= reference_length
                elif not 0 <= position < reference_length:
                    continue
                counts = self.reads_position_dict.setdefault((chrom, strand), {})
                counts[position] = counts.get(position, 0) + weight


    def normalize_read_counts(self, normalization, min_read_count_dict, normalization_scope="contig"):
        """Normalize positional counts, validating all factors before mutation."""
        validate_normalization_scope(normalization_scope, normalization)
        if normalization not in {"raw", "min", "mil"}:
            raise ValueError(f"Error: Given normalization method is not supported: {normalization}. Supported methods: (raw, min, mil)")
        if normalization == "raw":
            msg.message("Skipping normalization!")
            return
        if normalization == "min":
            targets = ([LIBRARY_MINIMUM_KEY] if normalization_scope == "library" else
                       sorted({chrom for chrom, _ in self.reads_position_dict}))
            for target in targets:
                if not min_read_count_dict or target not in min_read_count_dict:
                    label = "the library" if target == LIBRARY_MINIMUM_KEY else target
                    raise ValueError(f"No minimum read count found for {label}. Supply counts for every analyzed sample.")
                if min_read_count_dict[target] <= 0:
                    raise ValueError(
                        "Error: You chose min as normalization factor, but the provided "
                        f"minimal read count is not valid: {min_read_count_dict[target]}!"
                    )
        library_total = sum(self.no_accepted_reads_dict.values())
        for (chrom, _), counts in self.reads_position_dict.items():
            numerator = (min_read_count_dict[
                LIBRARY_MINIMUM_KEY if normalization_scope == "library" else chrom
            ] if normalization == "min" else 1000000)
            denominator = library_total if normalization_scope == "library" else self.no_accepted_reads_dict[chrom]
            factor = numerator / denominator
            for position in counts:
                counts[position] *= factor

    def to_wig(self, file_path):
        """
        Create two wiggle format files (+/-) for the used positions.
        """

        for (chrom, strand) in self.reads_position_dict:
            orientation = "reverse" if strand == "-" else "forward"
            # SAM reference names may legally contain slashes and ``..``. Keep
            # the exact name in the WIG declaration, but percent-encode it in
            # the filename so it cannot create nested or escaping paths.
            filename_chrom = quote(str(chrom), safe=WIG_FILENAME_SAFE)
            cur_file = file_path / (
                f"{self.alignment_file_path.stem}_{filename_chrom}_{orientation}.wig"
            )
            cur_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cur_file, "w", encoding="utf-8") as f:
                f.write(f"track type=wiggle_0 name={cur_file.name}\nvariableStep chrom={chrom} span=1\n")

                for position in sorted(self.reads_position_dict[(chrom,strand)].keys()):
                    value = self.reads_position_dict[(chrom,strand)][position]
                    f.write(f"{int(position)+1} {float(value)}\n")


    def output(self):
        """
        return the position dictionary and number of accepted reads dictionary"""
        return self.reads_position_dict, self.no_accepted_reads_dict

class IntervalReader():
    """
    Read sam/bam file and create an interlap dictionary for all reads or a selected amount of read lengths
    read_lengths : list (e.g. [30,32,40])
    rpkm_all_reads : if True, all mapped reads are used for calculating the RPKM
                     else, only the mapped reads of the given read-lengths are used
    """
    def __init__(self, alignment_file_path, read_length_dict, rpkm_all_reads,
                 alignment_policy=None):
        self.alignment_file_path = Path(alignment_file_path)
        self.read_length_dict = read_length_dict
        self.rpkm_all_reads = rpkm_all_reads
        self.alignment_policy = _resolve_alignment_policy(alignment_policy)
        self.alignment_diagnostics = {}

        self.wildcard = self.alignment_file_path.stem.split('_')[0]
        if self.read_length_dict is None:
            self.read_lengths = None
        elif self.wildcard in self.read_length_dict:
            self.read_lengths = self.read_length_dict[self.wildcard]
        elif "default" in self.read_length_dict:
            self.read_lengths = self.read_length_dict["default"]
        else:
            msg.warning("Warning no default value given for read-lengths. Using all read lengths.")
            self.read_lengths = None

        if self.read_lengths is None:
            self.rpkm_all_reads = True

        self.no_accepted_reads_dict = {}
        self.reads_interlap_dict = {}

        self._read_alignment_file()

    def _read_alignment_file(self):
        """Store aligned CIGAR blocks, with one identifier per alignment.

        Identifiers let read counting deduplicate overlapping blocks of the same
        read without collapsing distinct alignments with identical coordinates.
        """
        tmp_dict = {}
        for read_id, read in enumerate(iter_primary_alignments(
                self.alignment_file_path, self.alignment_policy, self.alignment_diagnostics)):
            chrom = read.reference_name
            if not self.rpkm_all_reads:
                read_length = alignment_query_length(read)
                if str(read_length) not in self.read_lengths:
                    continue
            blocks = read.get_blocks()
            if not blocks:
                continue
            strand = "-" if read.is_reverse else "+"
            self.no_accepted_reads_dict[chrom] = self.no_accepted_reads_dict.get(chrom, 0) + 1
            intervals = tmp_dict.setdefault((chrom, strand), [])
            intervals.extend((start, end - 1, read_id) for start, end in blocks if start < end)
        self.reads_interlap_dict = {key: interlap.InterLap(intervals) for key, intervals in tmp_dict.items()}

    def output(self):
        """
        return the interlap dictionary and number of accepted reads dictionary
        """
        return self.reads_interlap_dict, self.no_accepted_reads_dict
