#!/usr/bin/env python
import argparse
import os
import pandas as pd
import numpy as np

import collections
import csv
from collections import Counter, OrderedDict
from pathlib import Path

import lib.io as io
import lib.expression as expr

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
    te_map = {}
    for i, val in enumerate(xlsx_df.columns):
        if val.endswith("_peak_height"):
            peak_height_map[val[:-12]] = i
        elif val.endswith("_relative_density"):
            relative_density_map[val[:-17]] = i
        elif val.endswith("_rpkm"):
            rpkm_map[val[:-5]] = i
        elif val.endswith("_TE"):
            te_map[val[:-3]] = i

    for row in xlsx_df.itertuples(index=False, name=None):
        gene_type, unique_id = row[0], row[1]
        gene_name, codon_count = row[6], row[7]
        start_codon, stop_codon, nt_upstream, nt_seq, aa_seq = row[11:16]
        fiveprime, threeprime = row[19], row[20]

        meta_dict[unique_id] = (gene_type, gene_name, codon_count, start_codon, stop_codon, nt_upstream, nt_seq, aa_seq, fiveprime, threeprime)
        if unique_id not in dynamic_dict:
            dynamic_dict[unique_id] = {}

        for key, val in peak_height_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][0] = row[val]
            else:
                dynamic_dict[unique_id][key] = [row[val], np.nan, np.nan, np.nan, np.nan, np.nan]

        for key, val in relative_density_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][1] = row[val]
            else:
                dynamic_dict[unique_id][key] = [np.nan, np.nan, np.nan, row[val], np.nan, np.nan]

        for key, val in rpkm_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][2] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [np.nan, np.nan, np.nan, np.nan, float(row[val]), np.nan]

        for key, val in te_map.items():
            if key in dynamic_dict[unique_id]:
                dynamic_dict[unique_id][key][3] = float(row[val])
            else:
                dynamic_dict[unique_id][key] = [np.nan, np.nan, np.nan,  np.nan, np.nan, float(row[val])]

    return meta_dict, dynamic_dict

def get_log_fc_contrast(wildcards):
    """
    Create log fold change contrast list
    """

    contrasts = []

    wildcard_dict = {}
    for card in wildcards:
        if "-" not in card:
            continue

        method, condition, replicate = card.split("-")
        if "ribo" in method.lower():
            wildcard_dict[f"{condition}-{replicate}"] = []

    for card in wildcards:
        if "-" not in card:
            continue

        method, condition, replicate = card.split("-")
        if method.lower() in ["tts","tis"] and f"{condition}-{replicate}" in wildcard_dict:
            wildcard_dict[f"{condition}-{replicate}"].append(card)

    for ribo, tt_list in wildcard_dict.items():
        for val in tt_list:
            contrasts.append((f"RIBO-{ribo}", val))

    return contrasts

def calculate_fold_changes(row, con1, con2):
    """
    Calculate the contrast between two columns
    """

    con1_heights = float(row[con1 + "_peak_height"])
    con2_heights = float(row[con2 + "_peak_height"])

    if np.nan in [con1_heights, con2_heights]:
        return np.nan

    return np.log2(con2_heights / con1_heights)

def build_merged_dataframe(meta_dict, dynamic_dict):
    """
    Given the input data of all tables build a new dataframe with sorted wildcards
    """
    wildcards = set()
    for unique_id, val in dynamic_dict.items():
         wildcards.update(val.keys())

    wildcards = sorted(list(wildcards))
    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"] \
           + [card + "_peak_height" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]] \
           + ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"] \
           + [card + "_relative_density" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]] \
           + [card + "_rpkm" for card in wildcards] \
           + [card + "_TE" for card in expr.get_te_header(wildcards)]
    name_list = [f"s{x}" for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for unique_id, val in meta_dict.items():
        result = []
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        wild_dict = dynamic_dict[unique_id]
        result.extend([val[0], unique_id, chrom, int(start), int(stop), strand, val[1], val[2]])
        for card in wildcards:
            if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]:
                if card in wild_dict:
                    result.append(wild_dict[card][0])
                else:
                    result.append(np.nan)

        result.extend(val[3:])
        for card in wildcards:
            if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]:
                if card in wild_dict:
                    result.append(wild_dict[card][1])
                else:
                    result.append(np.nan)

        for card in wildcards:
            if card in wild_dict:
                result.append(wild_dict[card][2])
            else:
                result.append(np.nan)

        for card in expr.get_te_header(wildcards):
            if card in wild_dict:
                result.append(wild_dict[card][3])
            else:
                result.append(np.nan)
        result_rows.append(nTuple(*result))

    result_df = pd.DataFrame.from_records(result_rows, columns=header)
    contrasts = get_log_fc_contrast(wildcards)

    contrast_header = []
    if contrasts != []:
        tis_columns = [x for x in result_df.columns if ("TIS" in x.split("_")[0] and not "RNA" in x.split("_")[0]) ]
        result_df = result_df[result_df[tis_columns].any(axis="columns")]

        for contrast in contrasts:
            con1, con2 = contrast
            result_df[f"{con1}_{con2}_log2FC"] = result_df.apply(lambda row: calculate_fold_changes(row, con1, con2), axis=1)
            contrast_header.append(f"{con1}_{con2}_log2FC")


        new_header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"] \
                   + [card + "_peak_height" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]] \
                   + contrast_header \
                   + ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"] \
                   + [card + "_relative_density" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and not "RNA" in card.split("-")[0]] \
                   + [card + "_rpkm" for card in wildcards] \
                   + [card + "_TE" for card in expr.get_te_header(wildcards)]

        result_df = result_df[new_header]

    return result_df, wildcards


def screen_input_tables(table_list):
    """
    screen over the input tables
    """

    meta_dict, dynamic_dict = {}, {}
    for table in table_list:
        xlsx_df = pd.read_excel(table, sheet_name=None)["CDS"]

        meta_dict, dynamic_dict = extend_combined_dictionary(xlsx_df, meta_dict, dynamic_dict)

    return meta_dict, dynamic_dict

def write_merged_table(meta_dict, dynamic_dict, output_path):
    """
    create final merged table and write it to xlsx/csv file
    """
    df_res, _ = build_merged_dataframe(meta_dict, dynamic_dict)

    df_res.dropna(how="all", axis=1, inplace=True)
    df_res = df_res.sort_values(by=["Genome", "Start", "Stop", "Strand"])
    Path(os.path.dirname(output_path)).mkdir(parents=True, exist_ok=True)
    df_res.to_csv(output_path[:-4]+"csv", sep="\t", index=False, quoting=csv.QUOTE_NONE)

    io.excel_writer(output_path, {"CDS" : df_res})

def write_merged_gff(meta_dict, output_path):
    """
    create final merged annotation file in gff3 file
    """
    nTuple = collections.namedtuple('Pandas', ["chromosome", "source", "type", "start", "stop", "score", "strand", "phase", "attribute"])

    result_rows = []
    for unique_id, val in meta_dict.items():
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        attribute = "ID=%s;Name=%s" % (unique_id, val[1])
        result_rows.append(nTuple(chrom, "ORFBounder", "CDS", int(start), int(stop), ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(result_rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    Path(os.path.dirname(output_path)).mkdir(parents=True, exist_ok=True)

    with open(output_path.replace(".xlsx", ".gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(output_path.replace(".xlsx", ".gff"), "a") as f:
        df.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)



def merge_tables(table_list, output_path):
    """
    collect information from all input tables and merge them into one final output table
    """

    meta_dict, dynamic_dict = screen_input_tables(table_list)
    write_merged_table(meta_dict, dynamic_dict, output_path)
    write_merged_gff(meta_dict, output_path)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Merge all results for different files into one xlsx and csv file.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", required=True, help= "list of input tables.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    merge_tables(args.table_list, args.output_path)


if __name__ == '__main__':
    main()
