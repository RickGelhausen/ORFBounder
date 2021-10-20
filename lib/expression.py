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

RNAMAP = {"RIBO" : "RNA", "TIS" : "RNATIS", "TTS" : "RNATTS"}

def header_to_dictionary(method, condition, replicate, wildcards, cur_dict):
    """
    Add a header value to the correct sample in the given dictionary
    """
    if "%s-%s-%s" %(RNAMAP[method], condition, replicate) in wildcards:
        if (method, condition) in  cur_dict:
            cur_dict[(method, condition)].append(replicate)
        else:
            cur_dict[(method, condition)] = [replicate]

def get_te_header(wildcards):
    """
    generate the correct te_header based on the available data
    """
    te_header = []
    te_header_dict = OrderedDict()
    for card in wildcards:
        if "-" not in card or "rna" in card.lower():
            continue
        method, condition, replicate = card.split("-")
        header_to_dictionary(method, condition, replicate, wildcards, te_header_dict)

    for key, val in te_header_dict.items():
        te_header.extend(["%s-%s-%s" % (key[0], key[1], x) for x in val])

    return te_header


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

def te_value_to_dictionary(method, condition, replicate, read_dict, cur_dict):
    """
    add Translational Efficiency value to the correct dictionary entry
    """

    if (RNAMAP[method], condition, replicate) in read_dict:
        rpkm_ribo = read_dict[(method, condition, replicate)]
        rpkm_rna = read_dict[(RNAMAP[method], condition, replicate)]
        cur_TE = TE(rpkm_ribo, rpkm_rna)
        if (method, condition) in cur_dict:
            cur_dict[(method, condition)].append(cur_TE)
        else:
            cur_dict[(method, condition)] = [cur_TE]

def calculate_te(read_list, wildcards):
    """
    calculate the translational efficiency
    """
    read_dict = OrderedDict()
    te_dict = OrderedDict()
    for idx in range(len(wildcards)):
        method, condition, replicate = wildcards[idx].split("-")
        key = (method, condition, replicate)
        if key not in read_dict:
            read_dict[key] = read_list[idx]
        else:
            msg.warning("Warning: multiple equal keys in calculate_te")

    te_list = []
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

# def init_read_count_dict(read_count_dict, result_dict):
#     """
#     Collect intervals needed for read_counting.
#     adds all intervals that are not yet in the read_count_dict from the predictions_dict.
#     """
#     for offset in result_dict.keys():
#         for (chrom, strand) in result_dict[offset].keys():
#             for (start, stop) in result_dict[offset][(chrom, strand)].keys():
#                 if (chrom, start, stop, strand) not in read_count_dict:
#                     read_count_dict[(chrom, start, stop, strand)] = []

#     return read_count_dict

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
