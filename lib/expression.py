"""
module for expression calculations
"""

from collections import Counter, OrderedDict
from pathlib import Path

import sys
import pysam
import numpy as np
from interlap import InterLap

import lib.messaging as msg
from lib.alignment_reader import IntervalReader

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

    if f"{RNAMAP[method]}-{condition}-{replicate}" in wildcards:
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
    for card in wildcards:
        if "-" not in card:
            continue
        method, condition, replicate = card.split("-")
        if "rna" in method.lower():
            continue
        header_to_dictionary(method, condition, replicate, wildcards, te_header_dict)

    for key, val in te_header_dict.items():
        te_header.extend([f"{key[0]}-{key[1]}-{x}" for x in val])

    return te_header


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

    if (RNAMAP[method], condition, replicate) in read_dict:
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
    for idx, wildcard in enumerate(wildcards):
        method, condition, replicate = wildcard.split("-")
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            msg.warning("Warning: multiple equal keys in calculate_te")

    for key, val in read_dict.items():
        method, condition, replicate = key
        if "rna" in method.lower():
            continue
        te_value_to_dictionary(method, condition, replicate, read_dict, te_dict)

    te_list = []
    for key, val in te_dict.items():
        if len(val) > 1:
            t_eff = get_avg(val)
        else:
            t_eff = val
        te_list.extend(t_eff)

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


def create_interlap_dict(bam_file: Path) -> tuple[dict[tuple[str, str], InterLap], dict[str, int]]:
    """
    create a dictionary with interlap objects for the current bam file.
    """
    bam_path = Path(bam_file)
    msg.message(f"Reading: {bam_path}")

    interlap_dict = {}
    total_mapped_reads = {}
    tmp_dict = {}

    with pysam.AlignmentFile(bam_file) as samfile:
        try:
            for read in samfile.fetch():
                chrom = read.reference_name
                if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                    continue

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                if chrom in total_mapped_reads:
                    total_mapped_reads[chrom] += 1
                else:
                    total_mapped_reads[chrom] = 1

                if not read.is_reverse:
                    if (chrom, "+") in tmp_dict:
                        tmp_dict[(chrom, "+")].append((start, stop))
                    else:
                        tmp_dict[(chrom, "+")] = [(start, stop)]
                else:
                    if (chrom, "-") in tmp_dict:
                        tmp_dict[(chrom, "-")].append((start, stop))
                    else:
                        tmp_dict[(chrom, "-")] = [(start, stop)]

        except ValueError as exc:
            raise ValueError("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.") from exc

    for key, val in tmp_dict.items():
        inter = InterLap()
        inter.update(val)
        interlap_dict[key] = inter

    msg.success("Done")
    return interlap_dict, total_mapped_reads


def count_reads(
    chrom: str,
    start: int,
    stop: int,
    strand: str,
    read_interlap_dict: dict[tuple[str, str], InterLap]
) -> int:
    """
    count the reads falling into a certain region
    """
    return len(list(read_interlap_dict[(chrom, strand)].find((start, stop))))


def retrieve_read_counts(
    read_count_dict: dict[tuple[str, int, int, str], list],
    bam_files: list[Path],
    read_lengths: dict,
    all_reads_rpkm: bool
) -> tuple[dict[tuple[str, int, int, str], list], list[dict]]:
    """
    run over all available bam files and add read_counts for each interval in the interval dict.
    """
    accepted_read_list = []
    for bam_file in bam_files:
        interlap_dict, accepted_read_dict = IntervalReader(bam_file, read_lengths, all_reads_rpkm).output()

        accepted_read_list.append(accepted_read_dict)
        for (chrom, start, stop, strand) in read_count_dict.keys():
            read_count_dict[(chrom, start, stop, strand)].append(
                count_reads(chrom, start, stop, strand, interlap_dict)
            )

    return read_count_dict, accepted_read_list