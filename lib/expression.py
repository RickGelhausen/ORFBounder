import pysam
import collections
import pandas as pd
import numpy as np
import itertools as iter

from interlap import InterLap
from collections import Counter, OrderedDict
import lib.messaging as msg

class OrderedCounter(Counter, OrderedDict):
    pass

def get_TE_header(wildcards):
    """
    generate the correct TE_header based on the available data
    """
    TE_header = []
    TE_header_dict = OrderedDict()
    for card in wildcards:
        if "-" not in card:
            continue
        method, condition, replicate = card.split("-")
        if method == "TIS":
            if "%s-%s-%s" %("RNATIS", condition, replicate) in wildcards:
                if ("TIS", condition) in  TE_header_dict:
                    TE_header_dict[("TIS", condition)].append(replicate)
                else:
                    TE_header_dict[("TIS", condition)] = [replicate]
        elif method == "TTS":
            if "%s-%s-%s" %("RNATTS", condition, replicate) in wildcards:
                if ("TTS", condition) in  TE_header_dict:
                    TE_header_dict[("TTS", condition)].append(replicate)
                else:
                    TE_header_dict[("TTS", condition)] = [replicate]

    for key, val in TE_header_dict.items():
        TE_header.extend(["%s-%s-%s" % (key[0], key[1], x) for x in val])

    return TE_header


def calculate_rpkm(total_mapped, read_count, read_length):
    """
    calculate the rpkm
    """
    if read_length == 0:
        msg.warning("Warning: read_length: 0 detected! Setting RPKM to 0!")
        return 0
    elif total_mapped == 0:
        msg.warning("Warning: total_mapped: 0 detected! Setting RPKM to 0!")
        return 0

    return float("%.2f" % ((read_count * 1000000000) / (total_mapped * read_length)))

def TE(ribo_count, rna_count):
    """
    calculate the translational efficiency for one entry
    """

    if ribo_count == 0 and rna_count == 0:
        return np.nan
    elif rna_count == 0:
        return np.nan
    else:
        return ribo_count / rna_count

def get_avg(t_eff):
    """
    get the final TE list
    """

    valid_count = 0
    sum = 0
    for t in t_eff:
        if t != np.nan:
            valid_count += 1
            sum += t

    if valid_count == 0:
        t_eff.extend([np.nan])

    else:
        t_eff.extend([sum / valid_count])

    return t_eff

def calculate_TE(read_list, wildcards):
    """
    calculate the translational efficiency
    """
    read_dict = OrderedDict()
    TE_dict = OrderedDict()
    for idx in range(len(wildcards)):
        method, condition, replicate = wildcards[idx].split("-")
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            msg.warning("Warning: multiple equal keys in calculate_TE")

    TE_list = []
    for key, val in read_dict.items():
        method, condition, replicate = key

        if method == "TIS":
            if ("RNATIS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATIS", condition, replicate)]
                cur_TE = TE(rpkm_ribo, rpkm_rna)
                if ("TIS", condition) in TE_dict:
                    TE_dict[("TIS", condition)].append(cur_TE)
                else:
                    TE_dict[("TIS", condition)] = [cur_TE]

        elif method == "TTS":
            if ("RNATTS", condition, replicate) in read_dict:
                rpkm_ribo = read_dict[key]
                rpkm_rna = read_dict[("RNATTS", condition, replicate)]
                cur_TE = TE(rpkm_ribo, rpkm_rna)
                if ("TTS", condition) in TE_dict:
                    TE_dict[("TTS", condition)].append(cur_TE)
                else:
                    TE_dict[("TTS", condition)] = [cur_TE]

    TE_list = []
    for key, val in TE_dict.items():
        if len(val) > 1:
            t_eff = get_avg(val)
        else:
            t_eff = val
        TE_list.extend(t_eff)

    return TE_list


def init_read_count_dict(read_count_dict, result_dict):
    """
    Collect intervals needed for read_counting.
    adds all intervals that are not yet in the read_count_dict from the predictions_dict.
    """

    for (chrom, strand) in result_dict.keys():
        for (start, stop) in result_dict[(chrom, strand)].keys():
            if (chrom, start, stop, strand) not in read_count_dict:
                read_count_dict[(chrom, start, stop, strand)] = []

    return read_count_dict

def create_interlap_dict(bam_file):
    """
    create a dictionary with interlap objects for the current bam file.
    """

    print("Reading: %s" % bam_file)
    interlap_dict = {}
    total_mapped_reads = {}
    tmp_dict = {}

    samfile = pysam.AlignmentFile(bam_file)
    try:
        for read in samfile.fetch():
            chrom = read.reference_name
            if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                continue

            #start, stop = read.reference_start, read.reference_start + read.query_length
            start, stop = read.reference_start, read.reference_end-1
            #start, stop = read.query_alignment_start, read.query_alignment_end

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

    except ValueError:
        msg.error("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

    for key, val in tmp_dict.items():
        inter = InterLap()
        inter.update(val)
        interlap_dict[key] = inter

    msg.success("Done")
    return interlap_dict, total_mapped_reads

def count_reads(chrom, start, stop, strand, read_interlap_dict):
    """
    count the reads falling into a certain region
    """
    return len(list(read_interlap_dict[(chrom, strand)].find((start, stop))))

def retrieve_read_counts(read_count_dict, bam_files):
    """
    run over all available bam files and add read_counts for each interval in the interval dict.
    """

    total_mapped_list = []
    for idx in range(len(bam_files)):
        interlap_dict, total_mapped = create_interlap_dict(bam_files[idx])
        total_mapped_list.append(total_mapped)
        for (chrom, start, stop, strand) in read_count_dict.keys():
            read_count_dict[(chrom,start,stop,strand)].append(count_reads(chrom, start, stop, strand, interlap_dict))

    return read_count_dict, total_mapped_list
