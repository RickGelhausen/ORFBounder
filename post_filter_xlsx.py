#!/usr/bin/env python
import os,sys
import re
import argparse
import collections
import csv

import pandas as pd
import numpy as np

from pathlib import Path

import lib.misc as misc
import lib.expression as expr
import lib.io as io


def read_input_table(input_table):
    """
    Read .xlsx table and return a dataframe
    """

    return pd.read_excel(input_table, sheet_name=None)["CDS"]

# def create_idr_format_data_frame(xlsx_df, replicate):
#     """
#     create a dataframe for a given replicate, dest="genome"
#     """
#
#     replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
#     non_header = ["chromosome", "start", "stop", "Name", "score", "strand", "signalValue", "pvalue", "qvalue", "summit"]
#     nTuple = collections.namedtuple('Pandas', non_header)
#
#     rows = []
#     for row in xlsx_df.itertuples(index=False, name=None):
#         genome, start, stop, strand = row[2], int(row[3]), int(row[4]), row[5]
#         peak_height = row[replicate_index]
#         if np.isnan(peak_height) or peak_height == 0:
#             continue
#
#         # TODO MISSING TIS/TTS SPECIFICATIONS
#         if strand == "+":
#             out_start, out_stop = start-1, start+1
#         else:
#             out_start, out_stop = stop-3, stop-1
#
#         result = [genome, out_start, out_stop, ".", 0, strand, float(peak_height), -1, -1, -1]
#
#         rows.append(nTuple(*result))
#
#     return pd.DataFrame.from_records(rows, columns=non_header)
#
# def create_gff_format_data_frame(xlsx_df, replicate):
#     """
#     create a dataframe for a given replicate
#     """
#
#     replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
#     non_header = ["chrom", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"]
#     nTuple = collections.namedtuple('Pandas', non_header)
#
#     rows = []
#     for row in xlsx_df.itertuples(index=False, name=None):
#         unique_id, genome, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
#         peak_height = row[replicate_index]
#         if np.isnan(peak_height) or peak_height == 0:
#             continue
#         # TODO MISSING TIS/TTS SPECIFICATIONS
#         if strand == "+":
#             out_start, out_stop = start-1, start+1
#         else:
#             out_start, out_stop = stop-3, stop-1
#
#         result = [genome, "ORFBounder", "TSS", out_start, out_stop, peak_height, strand, ".", "ID=%s;" % (unique_id)]
#
#         rows.append(nTuple(*result))
#
#     return pd.DataFrame.from_records(rows, columns=non_header)
#
# def retrieve_replicates(xlsx_df):
#     """
#     Return a list of the replicates in the input file
#     """
#
#     return [entry[:-12] for entry in xlsx_df.columns if "_peak_height" in entry]
#
# # def write_dataframe_to_file(cur_df, tmp_folder, replicate):
# #     """
# #     write peak_height dataframe to file
# #     """
# #
# #     Path(os.path.dirname(tmp_folder)).mkdir(parents=True, exist_ok=True)
# #     cur_df.to_csv(os.path.join(tmp_folder, replicate + "_peak_height.gff"), sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
#
# def write_dataframe_to_file(cur_df, tmp_folder, replicate):
#     """
#     write peak_height dataframe to file
#     """
#
#     Path(os.path.dirname(tmp_folder)).mkdir(parents=True, exist_ok=True)
#     cur_df.to_csv(os.path.join(tmp_folder, replicate + "_peaks"), sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)
#
# def generate_peak_height_output_files(xlsx_df, tmp_folder):
#     """
#     Create an idr input file for each replicate in the input file
#     """
#
#     for replicate in retrieve_replicates(xlsx_df):
#         cur_df = create_idr_format_data_frame(xlsx_df, replicate)
#
#         write_dataframe_to_file(cur_df, tmp_folder, replicate)

def get_wildcards(xlsx_df):
    """
    retrieve wildcards from _peak_height columns
    """

    return [entry[:-12] for entry in xlsx_df.columns if "_peak_height" in entry]


def check_bamfile_input(bam_file_path, wildcards):
    """
    Check bam input path.
    Ensure that there is:
     - one bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if bam_file_path == "" or not isinstance(bam_file_path, str):
        return -1

    valid_bam = set()

    _, _, bam_file_list = next(os.walk(bam_file_path))
    bam_file_list = [ file for file in bam_file_list if file.endswith(".bam")]

    for card in wildcards:

        if "TIS" in card:
            RNATIS_prefix = "RNATIS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file or RNATIS_prefix in file:
                    valid_bam.add(os.path.join(bam_file_path,file))

        if "TTS" in card:
            RNATTS_prefix = "RNATTS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file or RNATTS_prefix in file:
                    valid_bam.add(os.path.join(bam_file_path,file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def retrieve_original_offset_position(start, stop, strand, offset, method):
    """
    retrieve the original position of the peaks
    """
    original_position = 0
    if method == "TIS":
        if strand == "+":
            original_position = start + offset
        else:
            original_position = stop - 2 - offset

    else:
        if strand == "+":
            original_position = stop - 2 + offset
        else:
            original_position = start - offset

    return original_position

def create_replicate_peak_dict(xlsx_df, replicate, offset_dict):
    """
    create dictionary for a given replicate with the original positions of the peak.
    """

    replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
    if "TIS" in card:
        method = "TIS"
    else:
        method = "TTS"

    replicate_dict = {}
    for row in xlsx_df.itertuples(index=False, name=None):
        unique_id, chrom, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
        peak_height = row[replicate_index]

        if np.isnan(peak_height) or peak_height == 0:
            continue

        offset = offset_dict[method][replicate]

        original_position = retrieve_original_offset_position(chrom, start-1, stop-1, strand, offset)

        replicate_dict[(chrom, original_position-2, original_position+2, strand)] = (unique_id, peak_height, 0)

    return replicate_dict

def write_test_file(read_count_dict, output_path, replicate):
    """
    write a test file
    """

    header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Original_peak_pos"] \
           + [replicate + "_peak_height", replicate + "_peak_rpkm"]

    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for (chrom, interval_start, interval_stop, strand), (unique_id, peak_height, rpkm) in read_count_dict.items():
        chrom, mid, strand = unique_id.split(":")
        start, stop = mid.split("-")

        result_rows.append(nTuple(unique_id, chrom, start, stop, strand, interval_start+2, peak_height, rpkm))


    Path(os.path.dirname(output_path)).mkdir(parents=True, exist_ok=True)

    io.excel_writer(output_path, {"CDS" : pd.DataFrame.from_records(result_rows, columns=header)})




def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Filter the final output file.')
    parser.add_argument("-i", "--input_table", action="store", dest="input_table", required=True, help= "Table created by ORFBounder merge or call_ORFBounder.py (Requires atleast 2 replicates).")
    parser.add_argument("--offset_json", action="store", dest="offset_json", required=True, help="A JSON file containing all information about the offsets.")
    parser.add_argument("--mapping_tis", action="store", dest="mapping_tis", default="", help="The mapping method used for the TIS data.")
    parser.add_argument("--mapping_tts", action="store", dest="mapping_tts", default="", help="The mapping method used for the TTS data.")
    parser.add_argument("--genome", action="store", dest="genome", required=True, help="The genome file.")
    parser.add_argument("--bamfiles", action="store", dest="bam_files", required=True)
    # parser.add_argument("-t","--tmp_folder", action="store", dest="tmp_folder", required=True, help="folder for storing temporary files.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    xlsx_df = read_input_table(args.input_table)
    offset_dict = misc.build_offset_dictionary(args.offset_json, args.genome, args.mapping_tis, args.mapping_tts)

    wildcards = get_wildcards(xlsx_df)
    bam_files = check_bamfile_input(args.bam_files, wildcards)

    for idx, card in enumerate(wildcards[0]):
        replicate_dict = get_original_peak_positions(xlsx_df, card, offset_dict)
        interlap_dict, total_mapped = misc.create_interlap_dict(bam_file[idx])
        for (chrom, start, stop, strand), val in read_count_dict.keys():
            read_count = expr.count_reads(chrom, start, stop, strand, interlap_dict)
            rpkm = expr.calculate_rpkm(total_mapped, read_count, stop-start+1)

            read_count_dict[(chrom,start,stop,strand)] = (val[0], val[1], rpkm)


if __name__ == '__main__':
    main()
