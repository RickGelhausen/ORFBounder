"""
module for expression calculations
"""

from collections import Counter, OrderedDict
from pathlib import Path

import numpy as np
from interlap import InterLap

import lib.messaging as msg
from lib.alignment_reader import IntervalReader, validate_normalization_scope
from lib.coordinates import split_circular_interval

class OrderedCounter(Counter, OrderedDict):
    pass

RNAMAP = {"RIBO" : "RNA", "TIS" : "RNATIS", "TTS" : "RNATTS"}

def header_to_dictionary(
    method: str,
    condition: str,
    replicate: str,
    wildcards: list[str],
    cur_dict: OrderedDict
) -> None:
    """
    Add a header value to the correct sample in the given dictionary
    """

    if method in RNAMAP and f"{RNAMAP[method]}-{condition}-{replicate}" in wildcards:
        if (method, condition) in cur_dict:
            cur_dict[(method, condition)].append(replicate)
        else:
            cur_dict[(method, condition)] = [replicate]


def get_te_header(wildcards: list[str]) -> list[str]:
    """
    generate the correct te_header based on the available data
    """

    te_header = []
    te_header_dict = OrderedDict()
    for card in dict.fromkeys(wildcards):
        parts = card.split("-")
        if len(parts) != 3 or not all(parts) or parts[0] not in RNAMAP:
            continue
        method, condition, replicate = parts
        header_to_dictionary(method, condition, replicate, wildcards, te_header_dict)

    for key, val in te_header_dict.items():
        te_header.extend([f"{key[0]}-{key[1]}-{x}" for x in val])

    return te_header


def normalization_read_count(accepted_reads: dict[str, int], chrom: str, normalization_scope: str = "contig") -> int:
    """Return the accepted-read denominator across both strands in the requested scope."""
    validate_normalization_scope(normalization_scope)
    if normalization_scope == "library":
        return sum(accepted_reads.values())
    return accepted_reads.get(chrom, 0)


def calculate_rpkm(total_mapped: int, read_count: int, read_length: int) -> float:
    """
    calculate the rpkm
    """

    if read_length == 0:
        msg.warning("Warning: read_length: 0 detected! Setting RPKM to 0!")
        return 0.0
    if total_mapped == 0:
        msg.warning("Warning: total_mapped: 0 detected! Setting RPKM to 0!")
        return 0.0

    return float(f"{(read_count * 1000000000) / (total_mapped * read_length):.2f}")


def te_value(ribo_count: float, rna_count: float) -> float:
    """
    calculate the translational efficiency for one entry
    """

    if rna_count == 0:
        return np.nan

    return ribo_count / rna_count


def get_avg(t_eff: list[float]) -> list[float]:
    """
    get the final TE list
    """

    valid_count = 0
    total = 0.0

    for t in t_eff:
        if not np.isnan(t):
            valid_count += 1
            total += t

    if valid_count == 0:
        t_eff.extend([np.nan])
    else:
        t_eff.extend([total / valid_count])

    return t_eff

def te_value_to_dictionary(
    method: str,
    condition: str,
    replicate: str,
    read_dict: OrderedDict,
    cur_dict: OrderedDict
) -> None:
    """
    add Translational Efficiency value to the correct dictionary entry
    """

    if method in RNAMAP and (RNAMAP[method], condition, replicate) in read_dict:
        rpkm_ribo = read_dict[(method, condition, replicate)]
        rpkm_rna = read_dict[(RNAMAP[method], condition, replicate)]
        cur_te = te_value(rpkm_ribo, rpkm_rna)
        if (method, condition) in cur_dict:
            cur_dict[(method, condition)].append(cur_te)
        else:
            cur_dict[(method, condition)] = [cur_te]


def calculate_te(read_list: list[float], wildcards: list[str]) -> list[float]:
    """
    calculate the translational efficiency
    """

    read_dict = OrderedDict()
    te_dict = OrderedDict()
    if len(read_list) != len(wildcards):
        raise ValueError("Read counts and sample names must have the same length for TE calculation.")
    accepted_methods = set(RNAMAP) | set(RNAMAP.values())
    for idx, wildcard in enumerate(wildcards):
        parts = wildcard.split("-")
        if len(parts) != 3 or not all(parts) or parts[0] not in accepted_methods:
            continue
        method, condition, replicate = parts
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            msg.warning("Warning: multiple equal keys in calculate_te")

    for key, val in read_dict.items():
        method, condition, replicate = key
        if method not in RNAMAP:
            continue
        te_value_to_dictionary(method, condition, replicate, read_dict, te_dict)

    te_list = []
    for val in te_dict.values():
        te_list.extend(val)

    return te_list


def init_read_count_dict(
    read_count_dict: dict[tuple[str, int, int, str], list],
    result_dict: dict[tuple[str, str], dict[tuple[int, int], any]]
) -> dict[tuple[str, int, int, str], list]:
    """
    Collect intervals needed for read_counting.
    adds all intervals that are not yet in the read_count_dict from the predictions_dict.
    """
    for (chrom, strand) in result_dict.keys():
        for (start, stop) in result_dict[(chrom, strand)].keys():
            if (chrom, start, stop, strand) not in read_count_dict:
                read_count_dict[(chrom, start, stop, strand)] = []

    return read_count_dict


def create_interlap_dict(bam_file: Path, alignment_policy=None) -> tuple[dict[tuple[str, str], InterLap], dict[str, int]]:
    """
    create a dictionary with interlap objects for the current bam file.
    """
    return IntervalReader(Path(bam_file), None, True, alignment_policy).output()



def count_reads(
    chrom: str,
    start: int,
    stop: int,
    strand: str,
    read_interlap_dict: dict[tuple[str, str], InterLap],
    reference_lengths: dict[str, int] | None = None,
    circular_contigs: set[str] | None = None,
) -> int:
    """
    count the reads falling into a certain region
    """
    intervals = read_interlap_dict.get((chrom, strand))
    if intervals is None:
        return 0
    circular_contigs = set(circular_contigs or ())
    if chrom in circular_contigs:
        if reference_lengths is None or chrom not in reference_lengths:
            raise ValueError(f"Missing reference length for circular sequence {chrom!r}.")
        query_intervals = split_circular_interval(start, stop, reference_lengths[chrom])
    else:
        query_intervals = ((start, stop),)
    # New alignment intervals carry a per-record identity so CIGAR blocks count
    # once per read, including when a virtual circular query is split at the
    # origin. Retain support for callers using plain (start, stop) pairs.
    record_ids = set()
    plain_multiplicity = Counter()
    for query_start, query_stop in query_intervals:
        matches = list(intervals.find((query_start, query_stop)))
        record_ids.update(match[2] for match in matches if len(match) >= 3)
        current_plain = Counter(match[:2] for match in matches if len(match) < 3)
        for match, count in current_plain.items():
            plain_multiplicity[match] = max(plain_multiplicity[match], count)
    return sum(plain_multiplicity.values()) + len(record_ids)


def retrieve_read_counts(
    read_count_dict: dict[tuple[str, int, int, str], list],
    bam_files: list[Path],
    read_lengths: dict,
    all_reads_rpkm: bool,
    alignment_policy=None,
    diagnostics: dict | None = None,
    reference_lengths: dict[str, int] | None = None,
    circular_contigs: set[str] | None = None,
) -> tuple[dict[tuple[str, int, int, str], list], list[dict]]:
    """
    run over all available bam files and add read_counts for each interval in the interval dict.
    """
    accepted_read_list = []
    for bam_file in bam_files:
        reader = IntervalReader(
            bam_file, read_lengths, all_reads_rpkm, alignment_policy,
        )
        interlap_dict, accepted_read_dict = reader.output()
        if diagnostics is not None:
            diagnostics[Path(bam_file).stem] = {
                **reader.alignment_diagnostics,
                "accepted_reads_by_contig": accepted_read_dict,
            }

        accepted_read_list.append(accepted_read_dict)
        for (chrom, start, stop, strand) in read_count_dict.keys():
            read_count_dict[(chrom, start, stop, strand)].append(
                count_reads(
                    chrom, start, stop, strand, interlap_dict,
                    reference_lengths=reference_lengths,
                    circular_contigs=circular_contigs,
                )
            )

    return read_count_dict, accepted_read_list
