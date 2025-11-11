#!/usr/bin/env python

"""
Module to merge multiple output tables into one final table
"""

from collections import Counter, OrderedDict
from pathlib import Path
import collections
import csv

import argparse
import pandas as pd
import numpy as np

from lib import io
from lib import messaging as msg
import lib.expression as expr

class OrderedCounter(Counter, OrderedDict):
    pass

def extend_combined_dictionary(
    xlsx_df: pd.DataFrame,
    meta_dict: dict,
    dynamic_dict: dict
) -> tuple[dict, dict]:
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


def get_log_fc_contrast(wildcards: list[str]) -> list[tuple[str, str]]:
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
        if method.lower() in ["tts", "tis"] and f"{condition}-{replicate}" in wildcard_dict:
            wildcard_dict[f"{condition}-{replicate}"].append(card)

    for ribo, tt_list in wildcard_dict.items():
        for val in tt_list:
            contrasts.append((f"RIBO-{ribo}", val))

    return contrasts

def calculate_fold_changes(row: pd.Series, tts: str, ribo: str, min_val: float) -> float:
    """
    Calculate the contrast between two columns
    """
    tts_heights = np.float64(row[f"{tts}_peak_height"])
    ribo_heights = np.float64(row[f"{ribo}_peak_height"])

    if np.isnan(tts_heights):
        return np.nan
    if (np.isnan(ribo_heights) or ribo_heights == 0) and tts_heights > 0:
        return np.log2(tts_heights / min_val)

    return np.log2(tts_heights / ribo_heights)

def build_merged_dataframe(
    meta_dict: dict,
    dynamic_dict: dict
) -> tuple[pd.DataFrame, list[str]]:
    """
    Given the input data of all tables build a new dataframe with sorted wildcards

    ### MAJOR CHANGE 1: Fixed operator precedence bug
    # Original: not "RNA" in card.split("-")[0]
    # This evaluates as: (not "RNA") in card.split("-")[0] -> True in [...] -> always True
    # Fixed: "RNA" not in card.split("-")[0]
    # This correctly checks if "RNA" is not in the split result
    """
    wildcards = set()
    for unique_id, val in dynamic_dict.items():
         wildcards.update(val.keys())

    wildcards = sorted(list(wildcards))
    contrasts = get_log_fc_contrast(wildcards)

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"] \
           + [f"{card}_peak_height" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]] \
           + ["Evidence"] \
           + ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"] \
           + [f"{card}_relative_density" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]] \
           + [f"{card}_rpkm" for card in wildcards] \
           + [f"{card}_TE" for card in expr.get_te_header(wildcards)]

    name_list = [f"s{x}" for x in range(len(header))]
    n_tuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for unique_id, val in meta_dict.items():
        result = []
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        wild_dict = dynamic_dict[unique_id]
        result.extend([val[0], unique_id, chrom, int(start), int(stop), strand, val[1], val[2]])

        for card in wildcards:
            if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]:
                if card in wild_dict:
                    result.append(wild_dict[card][0])
                else:
                    result.append(np.nan)

        evidence = []
        for card in wildcards:
            if ("TIS" in card or "TTS" in card) and "RNA" not in card.split("-")[0]:
                if card in wild_dict:
                    if wild_dict[card][0] > 0:
                        evidence.append(card)

        result.append(",".join(evidence))
        result.extend(val[3:])

        for card in wildcards:
            if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]:
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
        result_rows.append(n_tuple(*result))

    result_df = pd.DataFrame.from_records(result_rows, columns=header)

    heights_list = []
    for wild in wildcards:
        if "RNA" not in wild:
            heights_list.append(f"{wild}_peak_height")

    min_val = result_df[heights_list].min().min() * 0.9
    msg.message(f"Minimum value for log2FC calculation: {min_val}")

    contrast_header = []
    if contrasts:
        ts_columns = [x for x in result_df.columns if (("TIS" in x.split("_")[0] or "TTS" in x.split("_")[0]) and "RNA" not in x.split("_")[0])]
        result_df = result_df[result_df[ts_columns].any(axis="columns")]

        for contrast in contrasts:
            con1, con2 = contrast
            result_df[f"{con2}_{con1}_log2FC"] = result_df.apply(lambda row: calculate_fold_changes(row, con2, con1, min_val), axis=1)
            contrast_header.append(f"{con2}_{con1}_log2FC")

        new_header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"] \
                   + [f"{card}_peak_height" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]] \
                   + contrast_header \
                   + ["Evidence"] \
                   + ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"] \
                   + [f"{card}_relative_density" for card in wildcards if ("TIS" in card or "TTS" in card or "RIBO" in card) and "RNA" not in card.split("-")[0]] \
                   + [f"{card}_rpkm" for card in wildcards] \
                   + [f"{card}_TE" for card in expr.get_te_header(wildcards)]

        result_df = result_df[new_header]

    return result_df, wildcards


def screen_input_tables(table_list: list[Path]) -> tuple[dict, dict]:
    """
    screen over the input tables
    """
    meta_dict, dynamic_dict = {}, {}
    for table in table_list:
        xlsx_df = pd.read_excel(table, sheet_name=None)["CDS"]
        meta_dict, dynamic_dict = extend_combined_dictionary(xlsx_df, meta_dict, dynamic_dict)

    return meta_dict, dynamic_dict


def write_merged_table(
    meta_dict: dict,
    dynamic_dict: dict,
    output_path: Path
) -> pd.DataFrame:
    """
    create final merged table and write it to xlsx/csv file
    """

    df_res, _ = build_merged_dataframe(meta_dict, dynamic_dict)

    df_res.dropna(how="all", axis=1, inplace=True)
    df_res = df_res.sort_values(by=["Genome", "Start", "Stop", "Strand"])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    csv_path = output_path.with_suffix(".csv")
    df_res.to_csv(csv_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    io.excel_writer(output_path, {"CDS": df_res})
    return df_res


def write_merged_gff(res_df: pd.DataFrame, output_path: Path) -> None:
    """
    create final merged annotation file in gff3 file format
    """

    nTuple = collections.namedtuple('Pandas', ["chromosome", "source", "type", "start", "stop", "score", "strand", "phase", "attribute"])

    log2fc_cols = [(x, f"_{list(res_df.columns).index(x)}") for x in res_df.columns if "log2FC" in x]

    result_rows = []
    for row in res_df.itertuples(index=False):
        gene_type = getattr(row, "Type")
        identifier = getattr(row, "Identifier")
        chrom = getattr(row, "Genome")
        start = getattr(row, "Start")
        stop = getattr(row, "Stop")
        strand = getattr(row, "Strand")
        locus_tag = getattr(row, "Locus_tag")
        evidence = getattr(row, "Evidence")

        attribute = f"ID={identifier};Name={locus_tag};type={gene_type};evidence={evidence}"
        for col in log2fc_cols:
            log2fc = getattr(row, col[1])
            if not pd.isna(log2fc):
                attribute += f";{str(col[0]).lower()}={log2fc}"

        result_rows.append(nTuple(chrom, "ORFBounder", "CDS", int(start), int(stop), ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(result_rows, columns=["chromosome", "source", "type", "start", "stop", "score", "strand", "phase", "attribute"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    gff_path = output_path.with_suffix(".gff")

    with open(gff_path, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
    with open(gff_path, "a", encoding="utf-8") as f:
        df.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


def merge_tables(table_list: list[Path], output_path: Path) -> None:
    """
    collect information from all input tables and merge them into one final output table
    """
    meta_dict, dynamic_dict = screen_input_tables(table_list)
    res_df = write_merged_table(meta_dict, dynamic_dict, output_path)
    write_merged_gff(res_df, output_path)


def main() -> None:
    """
    Command line interface for merging multiple result tables into one final table
    """
    parser = argparse.ArgumentParser(description='Merge all results for different files into one xlsx and csv file.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", required=True, help="list of input tables.")
    parser.add_argument("-o", "--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    merge_tables(args.table_list, args.output_path)


if __name__ == '__main__':
    main()