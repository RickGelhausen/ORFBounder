#!/usr/bin/env python

import pandas as pd
import csv
import sys, os
import collections
import argparse
import re

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

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


def get_additional_ORF_information(genome_seq, cur_stop, strand, start_codons, stop_codons):
    """
    Search for stop upstream of current stop.
    Then search for start downstream from new stop.
    All done in-frame.
    """
    reverse_start_codons = [str(Seq(codon).reverse_complement()) for codon in start_codons]
    reverse_stop_codons = [str(Seq(codon).reverse_complement()) for codon in stop_codons]

    upstream_stop = 0
    longest_start = 0
    longest_nt = ""
    upstream_nt_stop = ""
    upstream_nt_start = ""

    if strand == "+":
        cur_position = cur_stop - 2
        nt=genome_seq[cur_position:cur_position+3]
        stop_nt_seq=nt
        while nt not in stop_codons:
            cur_position -= 3
            if cur_position < 0:
                upstream_stop = ""
                upstream_nt_stop = ""
                break;

            nt = genome_seq[cur_position:cur_position+3]
            stop_nt_seq += nt

        if upstream_stop == 0:
            upstream_stop = cur_position
            upstream_nt_stop = genome_seq[cur_position-15:cur_position]

        nt=genome_seq[cur_position:cur_position+3]
        while nt not in start_codons:
            cur_position += 3
            if cur_position == cur_stop - 2:
                return None

            nt = genome_seq[cur_position:cur_position+3]

        longest_start = cur_position
        longest_nt = genome_seq[longest_start:cur_stop+1]
        upstream_nt_start = genome_seq[longest_start-15:longest_start]

    else:
        cur_position = cur_stop
        nt=genome_seq[cur_position:cur_position+3]
        stop_nt_seq=nt
        while nt not in reverse_stop_codons:
            cur_position += 3
            if cur_position > len(genome_seq):
                upstream_stop = ""
                upstream_nt_stop = ""
                break;

            nt = genome_seq[cur_position:cur_position+3]
            stop_nt_seq += nt

        if upstream_stop == 0:
            upstream_stop = cur_position
            upstream_nt_stop = str(Seq(genome_seq[cur_position+2:cur_position+17]).reverse_complement())

        nt=genome_seq[cur_position:cur_position+3]
        while nt not in reverse_start_codons:
            cur_position -= 3
            if cur_position == cur_stop:
                return None

            nt = genome_seq[cur_position:cur_position+3]

        longest_start = cur_position + 2
        longest_nt = str(Seq(genome_seq[cur_stop:longest_start+1]).reverse_complement())
        upstream_nt_start = str(Seq(genome_seq[longest_start:longest_start+15]).reverse_complement())

    return longest_start, longest_nt, upstream_nt_start, upstream_stop, stop_nt_seq, upstream_nt_stop

def postprocess_excel_file(args, genome_dict):
    """
    Search for additional upstream stop codon and determine longest ORF
    """

    xlsx_df = pd.read_excel(args.in_xlsx, sheet_name=None)["all"]

    header = xlsx_df.columns
    rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        print(row)
        genome_id = getattr(row, "Genome")
        short_start = int(getattr(row, "Start")) - 1
        main_stop = int(getattr(row, "Stop")) - 1
        strand = getattr(row, "Strand")
        ORF_information = get_additional_ORF_information(genome_dict[genome_id][0], main_stop, strand, args.start_codons, args.stop_codons)
        if ORF_information != None:
            longest_start, longest_nt, upstream_nt_start, upstream_stop, stop_nt_seq, upstream_nt_stop = ORF_information


    # all_df = pd.DataFrame.from_records(rows, columns=header)
    # dataframe_dict = { "all" : all_df }
    #
    # excel_writer(args, dataframe_dict)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Post processing of the TTS xlsx file returned by the merge_TTS.py script.')
    parser.add_argument("-i", "--input_xlsx", action="store", dest="in_xlsx", required=True, help= "Output excel file.")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", required=True, help= "Genome file.")
    parser.add_argument("--target_site", action="store", dest="target_site", default="TTS", help="TTS")
    parser.add_argument("-o", "--output_xlsx", action="store", dest="out_xlsx", required=True, help= "Output excel file.")
    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])

    args = parser.parse_args()

    print("Fetching genome...")
    # read the genome file
    genome_file = SeqIO.parse(args.genome_file, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    print("Done.")

    postprocess_excel_file(args, genome_dict)



if __name__ == '__main__':
    main()
