#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import collections
import csv
from collections import Counter, OrderedDict
import itertools as iter
import numpy as np

class OrderedCounter(Counter, OrderedDict):
    pass

def df_to_dictionary(args, filepath, meta_dict):
    tmp_dict = {}
    table_df = pd.read_csv(filepath, sep="\t", comment="#")

    for row in table_df.itertuples(index=False, name='Pandas'):
        gene_type = getattr(row, "Type")
        unique_id = getattr(row, "Identifier")
        genome = getattr(row, "Genome")
        start = int(getattr(row, "Start"))
        stop = int(getattr(row, "Stop"))
        strand = getattr(row, "Strand")
        gene_name = getattr(row, "locus_tag")
        aa_count = int(getattr(row, "codon_count"))
        peak_height = float(getattr(row, "peak_height"))
        start_codon = getattr(row, "start_codon")
        stop_codon = getattr(row, "stop_codon")
        nt_upstream = getattr(row, "_11")
        nt_seq = getattr(row, "nt_seq")
        aa_seq = getattr(row, "aa_seq")
        relative_density = float(getattr(row, "relative_density"))
        fiveprime = getattr(row, "_15")
        threeprime = getattr(row, "_16")

        tmp_dict[unique_id] = (peak_height, relative_density)
        if unique_id not in meta_dict.keys():
            meta_dict[unique_id] = (gene_name, gene_type, aa_count, start_codon, stop_codon, nt_upstream, nt_seq, aa_seq, fiveprime, threeprime)

    return tmp_dict, meta_dict

#
def parse_input(args):
    """
    Check if each RIBO file has a TIS counterpart and create two lists of dictionaries as well as a meta dictionary
    """

    constrasts = [constrast.split("_") for constrast in args.contrast_list]

    table_path = os.path.dirname(os.path.abspath(args.table_list[0]))
    suffix = "." + ".".join(os.path.splitext(os.path.basename(args.table_list[0]))[0].split(".")[1:]) + ".csv"
    wildcards = [os.path.splitext(os.path.basename(x))[0].split(".")[0] for x in sorted(args.table_list)]

    meta_dict = {}
    wildcard_dicts = {}
    for card in wildcards:
        tmp_dict, meta_dict = df_to_dictionary(args, os.path.join(table_path, card + suffix), meta_dict)
        wildcard_dicts[card] = tmp_dict

    return meta_dict, wildcard_dicts, wildcards, constrasts

def calculate_fold_changes(unique_id, wildcard_dicts, contrasts):
    """
    calculate log2FC of the different samples
    """
    data = []

    for key, wild_dict in wildcard_dicts.items():
        if unique_id in wild_dict:
            peak_height = wild_dict[unique_id][0]
        else:
            peak_height = ""

        data.append(peak_height)

    for contrast in contrasts:
        con1, con2 = contrast
        if unique_id in wildcard_dicts[con1]:
            con1_height = wildcard_dicts[con1][unique_id][0]
        else:
            con1_height = ""

        if unique_id in wildcard_dicts[con2]:
            con2_height = wildcard_dicts[con2][unique_id][0]
        else:
            con2_height = ""


        if con1_height != "" and con2_height != "":
            fc = np.log2(con2_height / con1_height)
        else:
            fc = ""

        data.append(fc)

    return data


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


def create_wildcard_header(wildcards, contrasts):
    """
    create the header based on the given wildcards
    """
    header = ["%s_peak_height" % wildcard for wildcard in wildcards]
    header.extend(["%s_%s_log2FC" % (contrast[0], contrast[1]) for contrast in contrasts])
    return header


def get_relative_density(unique_id, wildcard_dicts):
    """
    extract relative_density
    """
    data = []
    for key, wild_dict in wildcard_dicts.items():
        if unique_id in wild_dict:
            rel_density = wild_dict[unique_id][1]
        else:
            rel_density = ""

        data.extend([rel_density])

    return data

def create_excel_file(args):
    """
    generate an overview excel file with the combined results
    """
    meta_dict, wildcard_dicts, wildcards, contrasts = parse_input(args)

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Gene_type", "Codon_count", "Start_codon", "Stop_codon"] +\
             create_wildcard_header(wildcards, contrasts) + ["15nt upstream", "Nucleotide_seq", "Aminoacid_seq"] +\
             ["%s_relative_density" % (sample) for sample in wildcards] + ["5'distance", "3'distance"]

    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    all_sheet = []
    for unique_id in meta_dict.keys():
        meta_info = meta_dict[unique_id]

        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        result = [unique_id, chrom, int(start), int(stop), strand, meta_info[0], meta_info[1], meta_info[2], meta_info[3], meta_info[4]] +\
                 calculate_fold_changes(unique_id, wildcard_dicts, contrasts) + [meta_info[5], meta_info[6], meta_info[7]] +\
                 get_relative_density(unique_id, wildcard_dicts) + [meta_info[8], meta_info[9]]

        all_sheet.append(nTuple(*result))
    all_df = pd.DataFrame.from_records(all_sheet, columns=[header[x] for x in range(len(header))])
    all_df = all_df.astype({"Start" : "int32", "Stop" : "int32"})
    all_df = all_df.sort_values(by=["Genome", "Start", "Stop"])
    dataframe_dict = { "all" : all_df }
    excel_writer(args, dataframe_dict)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='create an overview table containing merged RET results.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", required=True, help= "list of input tables.")
    parser.add_argument("--contrasts", nargs="+", dest="contrast_list", required=True, help= "list of contrasts used for log2FC.")
    parser.add_argument("-x", "--xlsx", action="store", dest="out_xlsx", required=True, help= "Output excel file.")

    args = parser.parse_args()

    create_excel_file(args)


if __name__ == '__main__':
    main()
