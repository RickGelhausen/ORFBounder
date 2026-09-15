#!/usr/bin/env python

from pathlib import Path
from numbers import Integral
from typing import Optional
from urllib.parse import quote, unquote
from dataclasses import dataclass


import csv
import json
import collections
import re
import pandas as pd
import pysam

from Bio import SeqIO


import lib.messaging as msg
from lib.alignment_reader import LIBRARY_MINIMUM_KEY
from lib.coordinates import parse_coordinate_id


GFF3_SEQID_SAFE = ".:^*$@!+_?-|"
_MAX_EXACT_BINARY64_INTEGER = 1 << 53
_MAX_EXPANDED_READ_LENGTHS = 100_000


def _unique_json_object(pairs):
    """Reject duplicate configuration keys instead of silently taking the last."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key {key!r}; use one value per sample or read length.")
        result[key] = value
    return result


@dataclass(frozen=True, slots=True)
class ReferenceContext:
    """Immutable parsed reference inputs safe to reuse within a batch grid.

    Genome sequences and annotation features are immutable values. Callers get
    fresh dictionaries/sets through the accessors below so run-local code cannot
    mutate cached topology state.
    """

    cache_key: tuple[tuple[str, str], tuple[str, str]]
    genome_items: tuple[tuple[str, str], ...]
    declared_reference_lengths: tuple[tuple[str, int], ...]
    circular_contigs: frozenset[str]
    annotation_features: tuple

    def genome_dict(self) -> dict[str, str]:
        return dict(self.genome_items)

    def annotation_topology(self) -> tuple[dict[str, int], set[str]]:
        return dict(self.declared_reference_lengths), set(self.circular_contigs)


def _fingerprinted_file_key(path: Path, input_records: list[dict]) -> tuple[str, str]:
    """Return the canonical path/content key already verified by provenance code."""
    canonical = str(Path(path).resolve(strict=True))
    for record in input_records:
        if record.get("path") == canonical:
            return canonical, str(record["sha256"])
    raise ValueError(f"No verified input fingerprint found for reference file: {path}")


def load_reference_context(
    genome_path: Path,
    annotation_path: Path,
    input_records: list[dict],
    cache: dict | None = None,
) -> ReferenceContext:
    """Load or reuse immutable FASTA/GFF parsing under verified content keys."""
    cache_key = (
        _fingerprinted_file_key(genome_path, input_records),
        _fingerprinted_file_key(annotation_path, input_records),
    )
    if cache is not None and cache_key in cache:
        context = cache[cache_key]
        if not isinstance(context, ReferenceContext):
            raise TypeError("Reference-context cache contains an invalid value.")
        return context

    from lib import misc

    genome = generate_genome_dict(genome_path)
    declared_lengths, circular_contigs = misc.parse_reference_topology(annotation_path)
    topology = (dict(declared_lengths), set(circular_contigs))
    annotation_features = tuple(
        misc.parse_annotation_features(annotation_path, reference_topology=topology)
    )
    context = ReferenceContext(
        cache_key=cache_key,
        genome_items=tuple(genome.items()),
        declared_reference_lengths=tuple(declared_lengths.items()),
        circular_contigs=frozenset(circular_contigs),
        annotation_features=annotation_features,
    )
    # Cache only a fully parsed, internally consistent reference pair.
    validate_reference_inputs(genome, annotation_path, [], reference_context=context)
    if cache is not None:
        cache[cache_key] = context
    return context


def _escape_gff3_seqid(seqid: str) -> str:
    """Escape a decoded sequence identifier for GFF3 column one."""
    return quote(str(seqid), safe=GFF3_SEQID_SAFE)


def reference_topology_from_dataframe(
    dataframe: pd.DataFrame,
) -> tuple[dict[str, int], set[str]]:
    """Read validated per-reference topology carried by ORF result rows."""
    reference_lengths = {
        str(chrom): int(length)
        for chrom, length in dataframe.attrs.get("reference_lengths", {}).items()
    }
    circular_contigs = set(map(str, dataframe.attrs.get("circular_contigs", ())))
    topology_by_chrom = {
        chrom: (length, chrom in circular_contigs)
        for chrom, length in reference_lengths.items()
    }
    topology_columns = {"Reference_length", "Is_circular"}.intersection(dataframe.columns)
    if not topology_columns:
        return reference_lengths, circular_contigs
    if topology_columns != {"Reference_length", "Is_circular"}:
        raise ValueError("Result tables must provide both Reference_length and Is_circular.")
    records = dataframe.to_dict(orient="records")
    if (records and all(pd.isna(row["Reference_length"]) for row in records)
            and all(isinstance(row["Is_circular"], (bool,)) and not row["Is_circular"]
                    for row in records)):
        # Merged legacy linear tables retain the new schema but cannot recover
        # lengths which were never present in their inputs.
        return reference_lengths, circular_contigs
    for row in records:
        chrom = str(row.get("Genome", unquote(str(row.get("chromosome", "")))))
        length, circular = row["Reference_length"], row["Is_circular"]
        if (isinstance(length, bool) or pd.isna(length)
                or int(length) != float(length) or int(length) < 1):
            raise ValueError(f"Invalid Reference_length for result sequence {chrom!r}.")
        if not isinstance(circular, (bool,)):
            # pandas generally exposes numpy booleans through records; avoid a
            # hard numpy dependency here while rejecting truthy strings.
            if type(circular).__name__ != "bool_":
                raise ValueError(f"Is_circular must be boolean for result sequence {chrom!r}.")
        length, circular = int(length), bool(circular)
        topology = (length, circular)
        previous_topology = topology_by_chrom.setdefault(chrom, topology)
        if previous_topology != topology:
            raise ValueError(f"Conflicting reference topology for result sequence {chrom!r}.")
        reference_lengths[chrom] = length
        if circular:
            circular_contigs.add(chrom)
        if {"Start", "Stop"}.issubset(row):
            start, stop = row["Start"], row["Stop"]
            if (pd.isna(start) or pd.isna(stop) or int(start) != float(start)
                    or int(stop) != float(stop)):
                raise ValueError(f"Invalid coordinates for result sequence {chrom!r}.")
            start, stop = int(start), int(stop)
            if not 1 <= start <= length or stop < start:
                raise ValueError(f"Invalid coordinates for result sequence {chrom!r}.")
            if (circular and stop - start + 1 > length) or (not circular and stop > length):
                raise ValueError(f"Coordinates exceed result sequence {chrom!r}.")
    return reference_lengths, circular_contigs


def gff3_topology_lines(
    reference_lengths: dict[str, int] | None,
    circular_contigs: set[str] | None,
) -> list[str]:
    """Return sequence-region and circular-region declarations for GFF3."""
    reference_lengths = dict(reference_lengths or {})
    circular_contigs = set(circular_contigs or ())
    lines = []
    for chrom in sorted(circular_contigs):
        if chrom not in reference_lengths or reference_lengths[chrom] < 1:
            raise ValueError(f"Missing reference length for circular sequence {chrom!r}.")
        length = int(reference_lengths[chrom])
        seqid = _escape_gff3_seqid(chrom)
        identifier = quote(f"{chrom}:region", safe=":.-_")
        lines.append(f"##sequence-region {seqid} 1 {length}\n")
        lines.append(
            f"{seqid}\tORFBounder\tregion\t1\t{length}\t.\t.\t.\t"
            f"ID={identifier};Is_circular=true\n"
        )
    return lines


def gff3_feature_segments(
    start: int,
    stop: int,
    strand: str,
    reference_length: int | None = None,
    circular: bool = False,
) -> tuple[tuple[int, int, str], ...]:
    """Convert a table interval into physical GFF3 rows and CDS phases.

    Result tables use virtual circular coordinates, but GFF3 feature columns
    stay on the declared sequence. Discontinuous CDS rows share their logical
    ID; phase is zero at the biological start and continues across the origin.
    """
    start, stop = int(start), int(stop)
    if strand not in {"+", "-"} or start < 1 or stop < start:
        raise ValueError(f"Invalid GFF3 feature coordinates: {start}-{stop}:{strand}")
    if reference_length is None:
        if circular:
            raise ValueError("Circular GFF3 features require a reference length.")
        return ((start, stop, "0"),)
    reference_length = int(reference_length)
    if reference_length < 1 or start > reference_length:
        raise ValueError(f"Invalid GFF3 feature coordinates: {start}-{stop}:{strand}")
    if not circular:
        if stop > reference_length:
            raise ValueError("A linear GFF3 feature cannot exceed its reference length.")
        return ((start, stop, "0"),)
    if stop - start + 1 > reference_length:
        raise ValueError("A circular GFF3 feature cannot span more than one revolution.")
    if stop <= reference_length:
        return ((start, stop, "0"),)

    high_segment = (start, reference_length)
    low_segment = (1, stop - reference_length)
    if strand == "+":
        continuation_phase = str((3 - ((reference_length - start + 1) % 3)) % 3)
        return ((*high_segment, "0"), (*low_segment, continuation_phase))
    continuation_phase = str((3 - ((stop - reference_length) % 3)) % 3)
    return ((*low_segment, "0"), (*high_segment, continuation_phase))


def validate_reference_inputs(
    genome: dict[str, str],
    annotation: Path,
    alignments: list[Path],
    *,
    reference_context: ReferenceContext | None = None,
) -> None:
    """Reject assembly/sequence-name mismatches before coverage or output work."""
    from lib.misc import parse_annotation_features, parse_reference_topology

    if not genome or any(not sequence for sequence in genome.values()):
        raise ValueError("Genome FASTA must contain nonempty sequences.")
    if reference_context is None:
        declared_lengths, circular_contigs = parse_reference_topology(annotation)
        annotation_features = parse_annotation_features(annotation)
    else:
        if not isinstance(reference_context, ReferenceContext):
            raise TypeError("reference_context must be a ReferenceContext instance.")
        if dict(reference_context.genome_items) != genome:
            raise ValueError("Reference context does not match the genome FASTA contents.")
        if reference_context.cache_key[1][0] != str(Path(annotation).resolve(strict=True)):
            raise ValueError("Reference context does not match the annotation path.")
        declared_lengths, circular_contigs = reference_context.annotation_topology()
        annotation_features = reference_context.annotation_features
    for chrom, length in declared_lengths.items():
        if chrom not in genome:
            raise ValueError(
                f"Annotation sequence-region {chrom!r} is missing from the genome FASTA."
            )
        if length != len(genome[chrom]):
            raise ValueError(
                f"Annotation sequence-region {chrom!r} has length {length}, but the genome FASTA has length {len(genome[chrom])}."
            )
    for feature in annotation_features:
        chrom, start, stop, strand = (
            feature.chromosome, feature.start, feature.stop, feature.strand)
        if chrom not in genome:
            raise ValueError(f"Annotation sequence {chrom!r} is missing from the genome FASTA. Use files from the same assembly.")
        reference_length = len(genome[chrom])
        valid_coordinates = (
            1 <= start <= reference_length
            and start <= stop
            and (stop <= reference_length if chrom not in circular_contigs
                 else stop - start + 1 <= reference_length)
        )
        if not valid_coordinates or strand not in {"+", "-"}:
            raise ValueError(
                "Invalid CDS coordinates or strand in annotation: "
                f"{feature.feature_id} ({chrom}:{start}-{stop}:{strand})")
    for path in dict.fromkeys(alignments):
        with pysam.AlignmentFile(path) as alignment:
            sequence_records = alignment.header.to_dict().get("SQ", ())
            seen_references = set()
            for record in sequence_records:
                chrom, length = record.get("SN"), record.get("LN")
                if not isinstance(chrom, str) or not chrom:
                    raise ValueError(
                        f"Alignment header in {path} contains a reference without a name."
                    )
                if chrom in seen_references:
                    raise ValueError(
                        f"Alignment header in {path} contains duplicate reference name {chrom!r}."
                    )
                seen_references.add(chrom)
                if (isinstance(length, bool) or not isinstance(length, Integral)
                        or length < 1):
                    raise ValueError(
                        f"Alignment reference {chrom!r} in {path} has invalid length {length!r}."
                    )
                if chrom not in genome or len(genome[chrom]) != length:
                    raise ValueError(f"Alignment reference {chrom!r} (length {length}) in {path} does not match the genome FASTA. Use the same reference assembly.")


def generate_genome_dict(genome_file: Path) -> dict[str, str]:
    """
    read a genome fasta file into a dictionary
    """
    genome_dict = {}
    # Passing a path directly to SeqIO.parse leaves its internally opened
    # handle to generator finalization.  Own the handle explicitly so it is
    # closed promptly on both complete iteration and parse/validation errors.
    with Path(genome_file).open("r", encoding="utf-8") as handle:
        for entry in SeqIO.parse(handle, "fasta"):
            if not entry.id:
                raise ValueError("FASTA sequence headers must include a nonempty identifier after '>'.")
            if str(entry.id) in genome_dict:
                raise ValueError(f"Duplicate FASTA sequence identifier: {entry.id}")
            genome_dict[str(entry.id)] = str(entry.seq).upper()

    return genome_dict


def parse_read_lengths(read_length_json: Path) -> Optional[dict[str, list[str]]]:
    """
    Parse the read length input into a continuous list form.
    """

    if read_length_json is None or not Path(read_length_json).is_file():
        msg.warning(
            "Warning: Empty read-lengths parameter given, using all available read lengths."
        )
        return None

    with open(read_length_json, "r", encoding="utf-8") as json_file:
        tmp_dict = json.load(json_file, object_pairs_hook=_unique_json_object)

    if not isinstance(tmp_dict, dict) or not tmp_dict:
        raise ValueError("Read-length JSON must be a nonempty object mapping sample names (or default) to lengths.")

    read_length_dict = {}
    for file, read_length_string in tmp_dict.items():
        if isinstance(read_length_string, bool):
            raise ValueError("Read lengths must be positive integers, not true/false.")
        if isinstance(read_length_string, int):
            read_length_dict[file] = [str(read_length_string)]

        elif isinstance(read_length_string, str):
            parts = read_length_string.split(",")
            read_lengths = set()
            for part in parts:
                if "-" in part:
                    interval = part.split("-")
                    i1, i2 = int(interval[0]), int(interval[1])
                    if i1 > i2:
                        i1, i2 = i2, i1
                    if i2 - i1 + 1 > _MAX_EXPANDED_READ_LENGTHS:
                        raise ValueError(
                            f"Read-length range for {file} exceeds "
                            f"{_MAX_EXPANDED_READ_LENGTHS:,} lengths; select a narrower range."
                        )
                    read_lengths.update(range(i1, i2 + 1))
                else:
                    read_lengths.add(int(part))

                if len(read_lengths) > _MAX_EXPANDED_READ_LENGTHS:
                    raise ValueError(
                        f"Read-length selection for {file} exceeds "
                        f"{_MAX_EXPANDED_READ_LENGTHS:,} lengths; select fewer lengths."
                    )

            read_length_dict[file] = [str(x) for x in sorted(read_lengths)]
        else:
            raise ValueError("Error: Read-length JSON file is not in correct format!")

        if any(int(length) <= 0 for length in read_length_dict[file]):
            raise ValueError(f"Read lengths for {file} must be positive integers. Omit read_length_json to use all lengths.")

    return read_length_dict


def validate_read_length_samples(
    read_length_dict: dict[str, list[str]] | None,
    alignment_paths,
) -> None:
    """Require an explicit sample entry or default when filtering is enabled.

    Alignment readers historically fall back to all lengths for an unknown
    sample. At user-facing validation boundaries that would turn a typo in an
    explicitly supplied filter into silently unfiltered data.
    """
    if read_length_dict is None or "default" in read_length_dict:
        return
    sample_names = sorted({
        Path(path).stem.split("_", 1)[0]
        for path in alignment_paths
        if path is not None
    })
    missing = [sample for sample in sample_names if sample not in read_length_dict]
    if missing:
        raise ValueError(
            "Read-length JSON has no entry or default for: "
            f"{', '.join(missing)}. Add each sample name or a top-level 'default'."
        )


def validate_offset_samples(offset_dict: dict, alignment_paths) -> None:
    """Require an explicit offset entry or default for every calling sample.

    Checking the complete input set before opening any alignment prevents a
    direct multi-assay run from writing coverage for earlier samples and then
    failing only when a later sample is missing from the offset configuration.
    """
    if "default" in offset_dict:
        return
    sample_names = sorted({
        Path(path).stem.split("_", 1)[0]
        for path in alignment_paths
        if path is not None
    })
    missing = [sample for sample in sample_names if sample not in offset_dict]
    if missing:
        raise ValueError(
            "Offset JSON has no entry or default for: "
            f"{', '.join(missing)}. Add each sample name or a top-level 'default'."
        )


def parse_total_reads(
    mapped_counts_file_path: Path,
    normalization_method: str,
    normalization_scope: str = "contig",
    required_samples=None,
) -> Optional[dict[str, int]]:
    """
    Read tab-separated accepted counts and determine each normalization minimum.

    Format: sample, chromosome, total_reads. When supplied, ``required_samples``
    prevents an omitted or misspelled calling sample from silently changing the
    normalization reference; additional intentional reference samples remain
    allowed.
    """

    min_read_count_dict = {}
    if normalization_method == "min":
        error_msg = (
            "Error: min normalization method chosen but no mapped_counts_file_path given!\n"
            "Either use a different normalization method or provide a file containing total read counts for each sample and each chromosome.\n"
            "Consider using our helper script to create the required files."
        )


        if mapped_counts_file_path is None or not Path(mapped_counts_file_path).is_file():
            raise FileNotFoundError(error_msg)

        sample_counts = {}
        seen = set()
        with open(mapped_counts_file_path, "r", encoding="utf-8") as f:
            lines = list(filter(None, [line.strip() for line in f.readlines()]))

            for line in lines:
                fields = line.split("\t")
                if len(fields) != 3:
                    raise ValueError("Mapped counts must have three tab-separated fields: sample, chromosome, count (no header).")
                sample, chrom, cur_count = (field.strip() for field in fields)
                if not sample or not chrom:
                    raise ValueError(f"Mapped-count sample and chromosome names must not be empty: {line}")
                try:
                    count = int(cur_count)
                except ValueError as exc:
                    raise ValueError(f"Mapped counts must be integers: {line}") from exc
                if count <= 0:
                    raise ValueError(f"Mapped counts must be positive for min normalization: {line}")
                key = (sample, chrom)
                if key in seen:
                    raise ValueError(f"Duplicate mapped-count row for sample {sample!r}, chromosome {chrom!r}.")
                seen.add(key)
                sample_counts[sample] = sample_counts.get(sample, 0) + count
                min_read_count_dict[chrom] = min(min_read_count_dict.get(chrom, count), count)

    else:
        return None

    if not min_read_count_dict:
        raise ValueError("Mapped counts file is empty.")
    if required_samples is not None:
        expected = {str(sample) for sample in required_samples}
        observed = set(sample_counts)
        missing = sorted(expected - observed)
        if missing:
            raise ValueError(
                "Mapped counts must cover every calling sample "
                f"(missing samples: {', '.join(missing)})."
            )
    if normalization_scope == "library":
        return {LIBRARY_MINIMUM_KEY: min(sample_counts.values())}
    if normalization_scope != "contig":
        raise ValueError("normalization_scope must be contig or library.")
    return min_read_count_dict


def parse_alignment_input(
    alignment_file_tis: Path | None,
    alignment_file_tts: Path | None
) -> str:
    """
    Check whether the input alignment files are valid and determine the execution method for ORFBounder
    """

    if alignment_file_tis and alignment_file_tts:
        if not alignment_file_tis.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TIS does not exist: {alignment_file_tis}"
            )

        if not alignment_file_tts.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TTS does not exist: {alignment_file_tts}"
            )

        return "combined_methods"

    if alignment_file_tis:
        if not alignment_file_tis.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TIS does not exist: {alignment_file_tis}"
            )
        return "TIS"

    if alignment_file_tts:
        if not alignment_file_tts.is_file():
            raise FileNotFoundError(
                f"Error: Non-empty alignment file path given for TTS does not exist: {alignment_file_tts}"
            )
        return "TTS"

    raise FileNotFoundError("Error: Please ensure to either provide a TIS file, a TTS file or both!")


def check_alignment_path_input(
    alignment_file_path: Path | None,
    alignment_file_tis: Path | None,
    alignment_file_tts: Path | None
) -> Optional[set[Path]]:
    """
    Check alignment input path.
    Ensure that there is:
     - one sam/bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if not alignment_file_path or not alignment_file_path.is_dir():
        return None

    identities = set()
    for path in (alignment_file_tis, alignment_file_tts):
        if path:
            parts = Path(path).stem.split("_")[0].split("-")
            if len(parts) != 3:
                raise ValueError(f"Expected METHOD-condition-replicate filename: {path}")
            method, condition, replicate = parts
            identities.update({f"{method}-{condition}-{replicate}",
                               f"RNA{method}-{condition}-{replicate}",
                               f"RIBO-{condition}-{replicate}",
                               f"RNA-{condition}-{replicate}"})
    valid_bam = {}
    for path in sorted(alignment_file_path.iterdir()):
        if not path.is_file() or path.suffix.lower() not in {".bam", ".sam"}:
            continue
        identity = path.stem.split("_")[0]
        if identity in identities:
            if identity in valid_bam:
                raise ValueError(f"Duplicate expression alignments for {identity}: {valid_bam[identity]} and {path}")
            valid_bam[identity] = path
    return set(valid_bam.values()) or None


def parse_offset_json(offset_json: Path) -> dict:
    """
    Read offset JSON file into a dictionary
    """

    if offset_json is None or not Path(offset_json).is_file():
        raise FileNotFoundError(
            f"Error: Offset JSON file does not exist! {offset_json}"
        )

    if offset_json.stat().st_size == 0:
        raise ValueError(
            f"Error: Offset JSON file is empty or invalid! {offset_json}"
        )

    with open(offset_json, "r", encoding="utf-8") as json_file:
        offsets = json.load(json_file, object_pairs_hook=_unique_json_object)
    if not isinstance(offsets, dict) or not offsets:
        raise ValueError("Offset JSON must be a nonempty object, e.g. {\"default\": {\"default\": 0}}.")
    for sample, lengths in offsets.items():
        if not isinstance(lengths, dict) or not lengths:
            raise ValueError(f"Offsets for {sample} must map read lengths (or default) to integer offsets.")
        for length, offset in lengths.items():
            # Readers look up str(read.query_length), so accepting keys such as
            # "030" would silently use a default offset instead of the intended
            # calibration for a 30-base read.
            if length != "default" and re.fullmatch(r"[1-9][0-9]*", length) is None:
                raise ValueError(
                    f"Invalid offset read length for {sample}: {length!r}; "
                    "use a positive integer without leading zeros."
                )
            if isinstance(offset, bool) or not isinstance(offset, int):
                raise ValueError(f"Offset for {sample}/{length} must be an integer.")
    return offsets


def write_gff_file(
    dataframe_out: pd.DataFrame,
    output_path: Path,
    output_filename: Path,
    reference_lengths: dict[str, int] | None = None,
    circular_contigs: set[str] | None = None,
) -> None:
    """
    write a dataframe to a gff file
    """

    file_path = output_path / output_filename
    file_path.parent.mkdir(parents=True, exist_ok=True)

    msg.message(f"Writing: {file_path}")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
        f.writelines(gff3_topology_lines(reference_lengths, circular_contigs))
    with open(file_path, "a", encoding="utf-8") as f:
        dataframe_out.to_csv(
            f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
        )


def write_codon_interval_gff(
    output_path: Path,
    output_basename: str,
    codon_dict: dict,
    reference_lengths: dict[str, int] | None = None,
    circular_contigs: set[str] | None = None,
) -> None:
    """
    Create a gff3 file with all codon intervals.
    """

    if len(codon_dict) == 0:
        raise ValueError("Error: Codon dictionary is empty, cannot create codon interval gff file!")

    nTuple_gff = collections.namedtuple(
        "Pandas",
        [
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue
        chrom, start, stop, strand = parse_coordinate_id(key)

        identifier = quote(str(key), safe=":.-_")
        codon = quote(str(val[0]), safe=":.-_")
        attribute = f"ID={identifier};Peak_height={val[1]};Name={codon};Start_codon={codon}"

        if chrom in set(circular_contigs or ()):
            segments = gff3_feature_segments(
                start + 1, stop + 1, strand,
                dict(reference_lengths or {}).get(chrom), True,
            )
        else:
            segment_start, segment_stop = max(1, start + 1), stop + 1
            reference_length = dict(reference_lengths or {}).get(chrom)
            if reference_length is not None:
                reference_length = int(reference_length)
                if reference_length < 1:
                    raise ValueError(f"Invalid reference length for sequence {chrom!r}.")
                segment_stop = min(reference_length, segment_stop)
            if segment_stop < segment_start:
                raise ValueError(
                    f"Codon interval {key!r} does not overlap its physical reference."
                )
            segments = ((segment_start, segment_stop, "."),)
        rows.extend(
            nTuple_gff(
                _escape_gff3_seqid(chrom), "ORFBounder", "codon_interval",
                segment_start, segment_stop, ".", strand, ".", attribute,
            )
            for segment_start, segment_stop, _ in segments
        )

    df = pd.DataFrame.from_records(
        rows,
        columns=[
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    write_gff_file(
        df, output_path, output_basename,
        reference_lengths=reference_lengths, circular_contigs=circular_contigs,
    )


def excel_writer(out_file_name: Path, data_frames: dict[str, pd.DataFrame]) -> None:
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only = [
        "Nucleotide_Seq",
        "Amino_Acid_Seq",
        "Start_codon",
        "Stop_codon",
        "Strand",
        "Codon_count",
    ]
    # Annotation labels are input data, not spreadsheet programs or links.
    # Keeping XlsxWriter's automatic conversions enabled would let a value such
    # as ``=1+1`` become a formula when a user opens an otherwise trusted result.
    writer = pd.ExcelWriter(
        out_file_name,
        engine="xlsxwriter",
        engine_kwargs={
            "options": {
                "strings_to_formulas": False,
                "strings_to_urls": False,
            }
        },
    )
    for sheetname, df in data_frames.items():
        # XLSX numeric cells are binary64-backed and silently round integers
        # above 2^53. If a column contains such a value, store all of its
        # integer cells as decimal text so raw statistical counts remain exact
        # and missing cells remain missing. The caller's DataFrame is unchanged.
        excel_df = df
        unsafe_integer_columns = [
            column for column in df
            if any(
                not isinstance(value, bool)
                and isinstance(value, Integral)
                and abs(int(value)) > _MAX_EXACT_BINARY64_INTEGER
                for value in df[column]
            )
        ]
        if unsafe_integer_columns:
            excel_df = df.copy()
            for column in unsafe_integer_columns:
                excel_df[column] = pd.array([
                    str(int(value))
                    if not isinstance(value, bool) and isinstance(value, Integral)
                    else value
                    for value in df[column]
                ], dtype=object)

        excel_df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(excel_df):
            series = excel_df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(len(str(series.name)), series.astype(str).str.len().max() if not series.empty else 0) + 1

            worksheet.set_column(idx, idx, max_len)
    writer.close()


def write_results_to_gff(
    result_df: pd.DataFrame,
    output_path: Path,
    output_basename: str,
    split_gff: bool
) -> None:
    """
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    reference_lengths, circular_contigs = reference_topology_from_dataframe(result_df)

    nTuple_gff = collections.namedtuple(
        "Pandas",
        [
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    for row in result_df.itertuples(index=False, name=None):
        (
            gene_type,
            identifier,
            chrom,
            start,
            stop,
            strand,
            locus_tag,
            codon_count,
            rpm_tis,
            rpm_tts,
            rpm_ribo,
            start_codon,
            stop_codon,
        ) = row[0:13]

        attributes = {
            "ID": identifier, "Name": locus_tag, "Peak_height_TIS": rpm_tis,
            "Peak_height_TTS": rpm_tts, "Peak_height_RIBO": rpm_ribo,
            "Start_codon": start_codon, "Stop_codon": stop_codon,
            "Codon_count": codon_count, "Type": gene_type,
        }
        attribute = ";".join(f"{key}={quote(str(value), safe=':.-_')}" for key, value in attributes.items()) + ";"
        segments = gff3_feature_segments(
            int(start), int(stop), strand,
            reference_lengths.get(chrom), chrom in circular_contigs,
        )
        cur_tuples = [
            nTuple_gff(
                _escape_gff3_seqid(chrom), "ORFBounder", "CDS",
                segment_start, segment_stop, ".", strand, phase, attribute,
            )
            for segment_start, segment_stop, phase in segments
        ]

        gff_all.extend(cur_tuples)
        if split_gff:
            if gene_type in {"Annotated", "Annotated_Complex"}:
                gff_annotated.extend(cur_tuples)
            elif gene_type == "Unannotated":
                gff_unannotated.extend(cur_tuples)
            elif gene_type == "Near_Annotated":
                gff_near_annotated.extend(cur_tuples)
            elif gene_type == "Internal_Inframe":
                gff_internal_inframe.extend(cur_tuples)
            elif gene_type == "N-terminal_extension":
                gff_n_terminal.extend(cur_tuples)
            elif gene_type == "Internal_OutofFrame":
                gff_internal_out.extend(cur_tuples)

    msg.message("Generating gff files...")

    df_all = pd.DataFrame.from_records(
        gff_all,
        columns=[
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "score",
            "strand",
            "phase",
            "attribute",
        ],
    )
    write_gff_file(
        df_all, output_path, Path(f"{output_basename}.gff"),
        reference_lengths, circular_contigs,
    )

    if split_gff:
        df_annotated = pd.DataFrame.from_records(
            gff_annotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_annotated,
            output_path,
            Path(f"{output_basename}_annotated.gff"),
            reference_lengths,
            circular_contigs,
        )

        df_unannotated = pd.DataFrame.from_records(
            gff_unannotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_unannotated,
            output_path,
            Path(f"{output_basename}_unannotated.gff"),
            reference_lengths,
            circular_contigs,
        )

        df_near_annotated = pd.DataFrame.from_records(
            gff_near_annotated,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_near_annotated,
            output_path,
            Path(f"{output_basename}_near_annotated.gff"),
            reference_lengths,
            circular_contigs,
        )

        df_internal_inframe = pd.DataFrame.from_records(
            gff_internal_inframe,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_internal_inframe,
            output_path,
            Path(f"{output_basename}_internal_inframe.gff"),
            reference_lengths,
            circular_contigs,
        )

        df_n_terminal = pd.DataFrame.from_records(
            gff_n_terminal,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_n_terminal,
            output_path,
            Path(f"{output_basename}_n_terminal.gff"),
            reference_lengths,
            circular_contigs,
        )

        df_internal_out = pd.DataFrame.from_records(
            gff_internal_out,
            columns=[
                "chromosome",
                "source",
                "type",
                "start",
                "stop",
                "score",
                "strand",
                "phase",
                "attribute",
            ],
        )
        write_gff_file(
            df_internal_out,
            output_path,
            Path(f"{output_basename}_internal_out.gff"),
            reference_lengths,
            circular_contigs,
        )
    msg.success("Done")


def write_results_to_table(
    df_results: pd.DataFrame,
    output_path: Path,
    output_basename: str
) -> None:
    """
    write a csv and xlsx file containing all information,
    """

    msg.message("Generating output_tables...")

    out_csv = output_path / f"{output_basename}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Minimal quoting is invisible for ordinary rows and preserves decoded GFF3
    # labels containing tabs, newlines or quotes as one round-trippable field.
    df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    out_xlsx = output_path / f"{output_basename}.xlsx"
    df_dict = {"CDS": df_results}
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)

    excel_writer(out_xlsx, df_dict)
    msg.success("Done")
