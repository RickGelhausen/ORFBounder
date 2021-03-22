#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd
import numpy as np

import collections
import csv
from collections import Counter, OrderedDict
from pathlib import Path

class OrderedCounter(Counter, OrderedDict):
    pass

def extend_combined_dictionary(xlsx_df, meta_dict, dynamic_dict):
    """
    collect data from current file and add it to the existing dictionary
    { chrom:start-stop:strand : { wildcard : peak_height, relative_density, RPKM, TE} }
    { chrom:start-stop:strand : metadata }
    """
    peak_height_map = {}
    relative_density_map = {}
    rpkm_map = {}
    TE_map = {}
    for i, val in enumerate(xlsx_df.columns):
        if val.endswith("_peak_height"):
            peak_height_map[val[:-12]] = i
        elif val.endswith("_relative_density"):
            relative_density_map[val[:-17]] = i
        elif val.endswith("_rpkm"):
            rpkm_map[val[:-5]] = i
        elif val.endswith("_TE"):
            TE_map[val[:-3]] = i

    for row in xlsx_df.itertuples(index=False, name=None):
        gene_type, unique_id = row[0], row[1]
        gene_name, codon_count = row[6], row[7]
        start_codon, stop_codon, nt_upstream, nt_seq, aa_seq = row[10:15]
        fiveprime, threeprime = row[17], row[18]

        meta_dict[unique_id] = (gene_type, gene_name, codon_count, start_codon, stop_codon, nt_upstream, nt_seq, aa_seq, fiveprime, threeprime)
        if unique_id not in dynamic_dict:
            dynamic_dict[unique_id] = {}

        for key, val in peak_height_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][0] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [float(row[val]), np.nan, np.nan, np.nan]

        for key, val in relative_density_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][1] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [np.nan, float(row[val]), np.nan, np.nan]

        for key, val in rpkm_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][2] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [np.nan, np.nan, float(row[val]), np.nan]

        for key, val in TE_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][3] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [np.nan, np.nan, np.nan, float(row[val])]

    return meta_dict, dynamic_dict

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

def build_merged_dataframe(meta_dict, dynamic_dict):
    """
    Given the input data of all tables build a new dataframe with sorted wildcards
    """
    wildcards = set()
    for unique_id, val in dynamic_dict.items():
         wildcards.update(val.keys())

    wildcards = sorted(list(wildcards))
    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"] \
           + [card + "_peak_height" for card in wildcards if ("TIS" or "TTS") and not "RNA" in card] \
           + ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"] \
           + [card + "_relative_density" for card in wildcards if ("TIS" or "TTS") and not "RNA" in card] \
           + [card + "_rpkm" for card in wildcards] \
           + [card + "_TE" for card in get_TE_header(wildcards)]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for unique_id, val in meta_dict.items():
        result = []
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        wild_dict = dynamic_dict[unique_id]
        result.extend([val[0], unique_id, chrom, int(start), int(stop), strand, val[1], val[2]])
        for card in wildcards:
            if ("TIS" or "TTS") and not "RNA" in card:
                if card in wild_dict:
                    result.append(wild_dict[card][0])
                else:
                    result.append(np.nan)

        result.extend(val[3:])
        for card in wildcards:
            if ("TIS" or "TTS") and not "RNA" in card:
                if card in wild_dict:
                    result.append(wild_dict[card][1])
                else:
                    result.append(np.nan)

        for card in wildcards:
            if card in wild_dict:
                result.append(wild_dict[card][2])
            else:
                result.append(np.nan)

        for card in get_TE_header(wildcards):
            if card in wild_dict:
                result.append(wild_dict[card][3])
            else:
                result.append(np.nan)
        result_rows.append(nTuple(*result))

    return pd.DataFrame.from_records(result_rows, columns=header), wildcards

def screen_input_tables(table_list):
    """
    screen over the input tables
    """

    meta_dict, dynamic_dict = {}, {}
    for table in table_list:
        xlsx_df = pd.read_excel(table, sheet_name=None)["CDS"]

        meta_dict, dynamic_dict = extend_combined_dictionary(xlsx_df, meta_dict, dynamic_dict)

    return meta_dict, dynamic_dict

def excel_writer(out_file_name, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Nucleotide_Seq", "Amino_Acid_Seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"]
    writer = pd.ExcelWriter(out_file_name, engine='xlsxwriter')
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

def merge_tables(table_list, output_path):
    """
    collect information from all input tables and merge them into one final output table
    """

    meta_dict, dynamic_dict = screen_input_tables(table_list)

    df_results, wildcards = build_merged_dataframe(meta_dict, dynamic_dict)

    #-empty_cols = [col for col in df_results.columns if list(df_results[col].unique()) == (["", nan])]
    #df_results.drop(empty_cols, axis=1, inplace=True)
    df_results.dropna(how="all", axis=1, inplace=True)
    df_results = df_results.sort_values(by=["Genome", "Start", "Stop", "Strand"])
    Path(os.path.dirname(output_path)).mkdir(parents=True, exist_ok=True)
    df_results.to_csv(output_path[:-4]+"csv", sep="\t", index=False, quoting=csv.QUOTE_NONE)

    excel_writer(output_path, {"CDS" : df_results}, wildcards)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Merge all results for different files into one xlsx and csv file.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", required=True, help= "list of input tables.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    merge_tables(args.table_list, args.output_path)


if __name__ == '__main__':
    main()
