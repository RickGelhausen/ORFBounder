#!/usr/bin/env python
# TTS_analysis/detect_longest_potential_ORFs.py
import pandas as pd
import csv
import sys, os
import collections
import argparse
import re

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import excel_utils as eu

def excel_writer(args, data_frames):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Aminoacid_seq", "Nucleotide_seq"]
    writer = pd.ExcelWriter(args.out_xlsx, engine='xlsxwriter')
    for sheetname, df in data_frames.items():
        df.to_excel(writer, sheet_name=sheetname, index=False)
        worksheet = writer.sheets[sheetname]
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.save()


def get_read_counts(args):
    """
    retrieve the reads and save them in a dictionary
    """

    read_df = pd.read_csv(args.read_counts, comment="#", header=None, sep="\t")

    read_count_dict = {}
    for row in read_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")

        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]
        read_count_dict["%s:%s-%s:%s" % (chromosome, start, stop, strand)] = read_list
    return read_count_dict


def calculate_expression_TIS(args, xlsx_df):
    header = list(xlsx_df.columns)

    dynamic_header_part1 = [x for x in header if "_peak_height" in x or "_log2FC" in x or "_avg" in x]
    dynamic_header_part2 = [x for x in header if "_relative_density" in x or "'distance" in x]

    total_mapped_dict = {}
    with open(args.total_mapped_reads, "r") as f:
        total = f.readlines()

    wildcards = []
    for line in total:
        wildcard, chromosome, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, chromosome)] = int(value)
        wildcards.append(wildcard)

    wildcards = eu.get_unique(wildcards)

    TE_header = eu.get_TE_header(wildcards)

    conditions = []
    for card in wildcards:
        conditions.append(card.split("-")[1])

    conditions = eu.get_unique(conditions)

    read_count_dict = get_read_counts(args)

    new_header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Gene_type",\
                  "Codon_count", "Start_codon", "Stop_codon"] + dynamic_header_part1 +\
                 ["15nt_upstream", "Nucleotide_seq", "Aminoacid_seq"] \
                  + dynamic_header_part2 + \
                 [cond + "_TE" for cond in TE_header] + [card + "_rpkm" for card in wildcards]

    name_list = ["s%s" % str(x) for x in range(len(new_header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    xlsx_df = xlsx_df.rename(columns={x:y for x,y in zip(xlsx_df.columns,range(0,len(xlsx_df.columns)))})
    rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        identifier = getattr(row, "_0")
        genome_id = getattr(row, "_1")
        start = int(getattr(row, "_2"))
        stop = int(getattr(row, "_3"))
        strand = getattr(row, "_4")
        locus_tag = getattr(row, "_5")
        gene_type = getattr(row, "_6")
        codon_count = int(getattr(row, "_7"))
        start_codon = getattr(row, "_8")
        stop_codon = getattr(row, "_9")
        #nucleotide_seq = getattr(row, "Nucleotide_seq")
        #aminoacid_seq = getattr(row, "Aminoacid_seq")
#        upstream_nt = getattr(row, "15nt_window")

        read_count = read_count_dict[identifier]

        length = codon_count * 3

        rpkm_list = []
        for idx, val in enumerate(read_count):
            rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], genome_id)], val, length))

        TE_list = eu.calculate_TE(rpkm_list, wildcards, conditions)

        result = [identifier, genome_id, int(start), int(stop), strand, locus_tag, gene_type, \
                  codon_count, start_codon, stop_codon] + \
                 [getattr(row, "_%s" % x) for x in range(10, 10+len(dynamic_header_part1))] + \
                 [getattr(row, "_%s" % (10+len(dynamic_header_part1))), getattr(row, "_%s" % (11+len(dynamic_header_part1))), getattr(row, "_%s" % (12+len(dynamic_header_part1)))] + \
                 [getattr(row, "_%s" % x) for x in range(13 + len(dynamic_header_part1), 13 + len(dynamic_header_part1)+len(dynamic_header_part2))] + \
                  TE_list + rpkm_list

        rows.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(rows, columns=new_header)
    dataframe_dict = { "all" : all_df }

    excel_writer(args, dataframe_dict)

def calculate_expression_TTS(args, xlsx_df):
    header = list(xlsx_df.columns)

    dynamic_header_part1 = [x for x in header if "_peak_height" in x or "_log2FC" in x or "_avg" in x]
    dynamic_header_part2 = [x for x in header if "_relative_density" in x or "'distance" in x]

    total_mapped_dict = {}
    with open(args.total_mapped_reads, "r") as f:
        total = f.readlines()

    wildcards = []
    for line in total:
        wildcard, chromosome, value = line.strip().split("\t")
        total_mapped_dict[(wildcard, chromosome)] = int(value)
        wildcards.append(wildcard)

    wildcards = eu.get_unique(wildcards)

    TE_header = eu.get_TE_header(wildcards)

    conditions = []
    for card in wildcards:
        conditions.append(card.split("-")[1])

    conditions = eu.get_unique(conditions)

    read_count_dict = get_read_counts(args)

    new_header = ["Identifier_short","Identifier_long", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Shortest_Gene_type",\
                  "Shortest_codon_count", "Shortest_start_codon", "position_alternativ_start", \
                  "Longest_Gene_type", "Longest_codon_count", "Longest_start_codon", "Stop_codon"] \
                  + dynamic_header_part1 + \
                 ["Shortest_15nt_upstream", "Shortest_Nucleotide_seq", "Shortest_Aminoacid_seq", \
                  "Longest_15nt_upstream", "Longest_Nucleotide_seq", "Longest_Aminoacid_seq", \
                  "Upstream_stop_codon", "Upstream_stop", "Stop_to_stop_nucleotide_seq"] + dynamic_header_part2 + \
                 ["Shortest_" + cond + "_TE" for cond in TE_header] + ["Shortest_" + card + "_rpkm" for card in wildcards] +\
                 ["Longest_" + cond + "_TE" for cond in TE_header] + ["Longest_" + card + "_rpkm" for card in wildcards]
    name_list = ["s%s" % str(x) for x in range(len(new_header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        genome_id = getattr(row, "Genome")
        start = int(getattr(row, "Start"))
        stop = int(getattr(row, "Stop"))
        strand = getattr(row, "Strand")
        locus_tag = getattr(row, "Locus_tag")
        shortest_gene_type = getattr(row, "Shortest_Gene_type")
        shortest_codon_count = int(getattr(row, "Shortest_codon_count"))
        longest_codon_count = int(getattr(row, "Longest_codon_count"))
        shortest_start_codon = getattr(row, "Shortest_start_codon")
        longest_start_codon = getattr(row, "Longest_start_codon")
        position_alternativ_start = int(getattr(row, "Position_alternativ_start"))
        longest_gene_type = getattr(row, "Longest_Gene_type")
        stop_codon = getattr(row, "Stop_codon")
        shortest_15nt_upstream = getattr(row, "Shortest_15nt_upstream")
        longest_15nt_upstream = getattr(row, "Longest_15nt_upstream")
        shortest_nucleotide_seq = getattr(row, "Shortest_Nucleotide_seq")
        shortest_aminoacid_seq = getattr(row, "Shortest_Aminoacid_seq")
        longest_nucleotide_seq = getattr(row, "Longest_Nucleotide_seq")
        longest_aminoacid_seq = getattr(row, "Longest_Aminoacid_seq")
        upstream_stop_codon = getattr(row, "Upstream_stop_codon")
        upstream_stop = getattr(row, "Upstream_stop")
        stop_to_stop_nucleotide_seq = getattr(row, "Stop_to_stop_nucleotide_seq")

        short_id = "%s:%s-%s:%s" % (genome_id, start, stop, strand)
        if strand == "+":
            long_id = "%s:%s-%s:%s" % (genome_id, position_alternativ_start, stop, strand)
        else:
            long_id = "%s:%s-%s:%s" % (genome_id, start, position_alternativ_start, strand)

        short_read_count = read_count_dict[short_id]
        long_read_count = read_count_dict[long_id]

        short_length = shortest_codon_count * 3
        long_length = longest_codon_count * 3

        short_rpkm_list = []
        for idx, val in enumerate(short_read_count):
            short_rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], genome_id)], val, short_length))

        short_TE_list = eu.calculate_TE(short_rpkm_list, wildcards, conditions)

        long_rpkm_list = []
        for idx, val in enumerate(long_read_count):
            long_rpkm_list.append(eu.calculate_rpkm(total_mapped_dict[(wildcards[idx], genome_id)], val, long_length))

        long_TE_list = eu.calculate_TE(long_rpkm_list, wildcards, conditions)

        result = [short_id, long_id, genome_id, int(start), int(stop), strand, locus_tag, shortest_gene_type, \
                  shortest_codon_count, shortest_start_codon, position_alternativ_start, \
                  longest_gene_type, longest_codon_count, longest_start_codon, stop_codon] + \
                 [getattr(row, "_%s" % x) for x in range(14, 14+len(dynamic_header_part1))] + \
                 [shortest_15nt_upstream, shortest_nucleotide_seq, shortest_aminoacid_seq, \
                  longest_15nt_upstream, longest_nucleotide_seq, longest_aminoacid_seq, \
                  upstream_stop_codon, upstream_stop, stop_to_stop_nucleotide_seq] + \
                 [getattr(row, "_%s" % x) for x in range(23 + len(dynamic_header_part1), 23 + len(dynamic_header_part1)+len(dynamic_header_part2))] + \
                  short_TE_list + short_rpkm_list + long_TE_list + long_rpkm_list

        rows.append(nTuple(*result))

    all_df = pd.DataFrame.from_records(rows, columns=new_header)
    all_df = all_df.sort_values(by=["Start", "Stop"])
    dataframe_dict = { "all" : all_df }

    excel_writer(args, dataframe_dict)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Post processing of the TTS xlsx file returned by the merge_TTS.py script. Check for the longest potential ORF using the upstream stop.')
    parser.add_argument("-i", "--input_xlsx", action="store", dest="in_xlsx", required=True, help="Input excel file.")
    parser.add_argument("-m", "--mapped_reads", action="store", dest="total_mapped_reads", required=True, help="file with total mapped reads for different samples.")
    parser.add_argument("-r", "--read_counts", action="store", dest="read_counts", required=True, help="gff file containing readcounts for each feature/sample")
    parser.add_argument("-o", "--output_xlsx", action="store", dest="out_xlsx", required=True, help= "Output excel file.")

    args = parser.parse_args()

    xlsx_df = pd.read_excel(args.in_xlsx, sheet_name=None)["all"]
    header = list(xlsx_df.columns)

    if "Position_alternativ_start" in header:
        calculate_expression_TTS(args, xlsx_df)
    else:
        calculate_expression_TIS(args, xlsx_df)




if __name__ == '__main__':
    main()
