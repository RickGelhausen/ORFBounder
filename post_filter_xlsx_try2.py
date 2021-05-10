#!/usr/bin/env python
import os,sys
import re
import argparse
import collections
import csv

import pandas as pd
import numpy as np
import pysam
import interlap

from pathlib import Path

import lib.misc as misc
import lib.expression as expr
import lib.io as io
import lib.messaging as msg

import matplotlib.pyplot as plt

import scipy.stats

def read_data(data_file):
    with open(data_file, "r") as f:
        data = [float(x.rstrip("\n")) for x in f.readlines()]

    return data


def read_input_table(input_table):
    """
    Read .xlsx table and return a dataframe
    """

    return pd.read_excel(input_table, sheet_name=None)["CDS"]

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

    valid_bam = []

    _, _, bam_file_list = next(os.walk(bam_file_path))
    bam_file_list = [ file for file in bam_file_list if file.endswith(".bam")]

    for card in wildcards:

        if "TIS" in card:
            #RNATIS_prefix = "RNATIS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file and not "RNA" in file:# or RNATIS_prefix in file:
                    valid_bam.append(os.path.join(bam_file_path,file))

        if "TTS" in card:
            #RNATTS_prefix = "RNATTS-" + "-".join(card.split("-")[1:])
            for file in bam_file_list:
                if card in file and not "RNA" in file:# or RNATTS_prefix in file:
                    valid_bam.append(os.path.join(bam_file_path,file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def retrieve_original_offset_position(start, stop, strand, offset, method):
    """
    retrieve the original position of the peaks
    """

    peak_position = 0
    if method == "TIS":
        if strand == "+":
            peak_position = start + offset
        else:
            peak_position = stop - offset

    else:
        if strand == "+":
            peak_position = stop + offset
        else:
            peak_position = start - offset

    return peak_position

def create_replicate_peak_dict(xlsx_df, replicate, offset_dict):
    """
    create dictionary for a given replicate with the original positions of the peak.
    """

    replicate_index = list(xlsx_df.columns).index(replicate+"_peak_height")
    if "TIS" in replicate:
        method = "TIS"
    else:
        method = "TTS"

    replicate_dict = {}
    for row in xlsx_df.itertuples(index=False, name=None):
        unique_id, chrom, start, stop, strand = row[1], row[2], int(row[3]), int(row[4]), row[5]
        peak_height = row[replicate_index]

        try:
            offset = offset_dict[method][replicate]
        except KeyError:
            offset = offset_dict[method]["default"]

        peak_position = retrieve_original_offset_position(int(start)-1, int(stop)-1, strand, int(offset), method)

        if np.isnan(peak_height):
            peak_height = 0

        replicate_dict[(chrom, peak_position-2, peak_position+2, strand)] = (unique_id, peak_height, 0, 0)

    return replicate_dict

def prepare_dataframe(xlsx_df):
    """
    create dataframe with results
    """

    header = ["Type","Identifier", "Genome", "Start", "Stop", "Strand", "Codon_count"]

    return xlsx_df[header].copy()

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
            currentreadlength = read.query_length
            if currentreadlength not in [31,32]:
                continue

            if chrom in total_mapped_reads:
                total_mapped_reads[chrom] += 1
            else:
                total_mapped_reads[chrom] = 1

            if not read.is_reverse:
                # if start > 1619730 and start < 1619874:
                #     print(start, stop)
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
        inter = interlap.InterLap()
        inter.update(val)
        interlap_dict[key] = inter

    msg.success("Done")
    return interlap_dict, total_mapped_reads


def calculate_tpm(read_count_list, length_list):
    """
    """

    tmp = [x/y for x,y in zip(read_count_list, length_list)]

    return [1000000*(i / sum(tmp)) for i in tmp]

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Filter the final output file.')
    parser.add_argument("-i", "--input_table", action="store", dest="input_table", required=True, help= "Table created by ORFBounder merge or call_ORFBounder.py (Requires atleast 2 replicates).")
    parser.add_argument("--offset_json", action="store", dest="offset_json", required=True, help="A JSON file containing all information about the offsets.")
    parser.add_argument("--mapping_tis", action="store", dest="mapping_tis", default="", help="The mapping method used for the TIS data.")
    parser.add_argument("--mapping_tts", action="store", dest="mapping_tts", default="", help="The mapping method used for the TTS data.")
    parser.add_argument("--genome", action="store", dest="genome", required=True, help="The genome file.")
    parser.add_argument("--bamfiles", action="store", dest="bam_files", required=True)
    parser.add_argument("--validated", action="store", dest="validated_orfs", required=True)

    # parser.add_argument("-t","--tmp_folder", action="store", dest="tmp_folder", required=True, help="folder for storing temporary files.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    xlsx_df = read_input_table(args.input_table)
    offset_dict = misc.build_offset_dictionary(args.offset_json, args.genome, args.mapping_tis, args.mapping_tts)

    wildcards = get_wildcards(xlsx_df)
    bam_files = list(check_bamfile_input(args.bam_files, wildcards))

    new_df = prepare_dataframe(xlsx_df)

    with open(args.validated_orfs, "r") as f:
        validated_orfs = [x.strip() for x in f.readlines() if x != ""]

    print(validated_orfs)

    for idx, card in enumerate(wildcards):
        read_count_dict = create_replicate_peak_dict(xlsx_df, card, offset_dict)
        interlap_dict, total_mapped = create_interlap_dict(bam_files[idx])
        peak_height_list = []
        rpkm_list = []
        read_count_list = []
        length_list = []
        for (chrom, start, stop, strand), val in read_count_dict.items():
            read_count = expr.count_reads(chrom, start, stop, strand, interlap_dict)
            rpkm = expr.calculate_rpkm(total_mapped[chrom], read_count, int(stop)-int(start)+1)

            read_count_dict[(chrom,start,stop,strand)] = (val[0], val[1], rpkm, read_count)

            length_list.append(stop-start+1)
            peak_height_list.append(val[1])
            read_count_list.append(read_count)
            rpkm_list.append(rpkm)

        tpm_list = calculate_tpm(read_count_list, length_list)

        new_df["%s_peak_height" % card] = peak_height_list
        new_df["%s_read_count" % card] = read_count_list
        new_df["%s_peak_rpkm" % card] = rpkm_list
        new_df["%s_peak_tpm" % card] = tpm_list

    l_thresh = 30
    u_thresh = 5000
    new_df = new_df.loc[(new_df["TIS-dcmeB-1_peak_rpkm"] > l_thresh) & (new_df["TIS-dcmeB-2_peak_rpkm"] > l_thresh) & (new_df["TIS-dcmeB-3_peak_rpkm"] > l_thresh)]
    new_df = new_df.loc[(new_df["TIS-dcmeB-1_peak_rpkm"] < u_thresh) & (new_df["TIS-dcmeB-2_peak_rpkm"] < u_thresh) & (new_df["TIS-dcmeB-3_peak_rpkm"] < u_thresh)]

    rep1_rpkm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-1_peak_rpkm"].to_list())
    rep2_rpkm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-2_peak_rpkm"].to_list())
    rep3_rpkm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-3_peak_rpkm"].to_list())

    rep1_read_count = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-1_read_count"].to_list())
    rep2_read_count = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-2_read_count"].to_list())
    rep3_read_count = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-3_read_count"].to_list())

    rep1_peak_height = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-1_peak_height"].to_list())
    rep2_peak_height = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-2_peak_height"].to_list())
    rep3_peak_height = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-3_peak_height"].to_list())

    rep1_peak_tpm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-1_peak_tpm"].to_list())
    rep2_peak_tpm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-2_peak_tpm"].to_list())
    rep3_peak_tpm = np.array(new_df.loc[new_df["Strand"].isin(["+","-"]), "TIS-dcmeB-3_peak_tpm"].to_list())
    #
    print("--------------------------------------------------------------------------------------------------------")
    print("RPKM_peak r1 v r2: ",scipy.stats.wilcoxon(rep1_rpkm,rep2_rpkm))
    print("RPKM_peak r1 v r3: ",scipy.stats.wilcoxon(rep1_rpkm,rep3_rpkm))
    print("RPKM_peak r2 v r3: ",scipy.stats.wilcoxon(rep2_rpkm,rep3_rpkm))

    print("Read_count r1 v r2: ",scipy.stats.wilcoxon(rep1_read_count,rep2_read_count))
    print("Read_count r1 v r3: ",scipy.stats.wilcoxon(rep1_read_count,rep3_read_count))
    print("Read_count r2 v r3: ",scipy.stats.wilcoxon(rep2_read_count,rep3_read_count))

    print("Peak_height r1 v r2: ",scipy.stats.wilcoxon(rep1_peak_height,rep2_peak_height))
    print("Peak_height r1 v r3: ",scipy.stats.wilcoxon(rep1_peak_height,rep3_peak_height))
    print("Peak_height r2 v r3: ",scipy.stats.wilcoxon(rep2_peak_height,rep3_peak_height))

    print("TPM_peak r1 v r2: ",scipy.stats.wilcoxon(rep1_peak_tpm,rep2_peak_tpm))
    print("TPM_peak r1 v r3: ",scipy.stats.wilcoxon(rep1_peak_tpm,rep3_peak_tpm))
    print("TPM_peak r2 v r3: ",scipy.stats.wilcoxon(rep2_peak_tpm,rep3_peak_tpm))

    print("--------------------------------------------------------------------------------------------------------")
    print("RPKM_peak r1 v r2: ",scipy.stats.ks_2samp(rep1_rpkm,rep2_rpkm))
    print("RPKM_peak r1 v r3: ",scipy.stats.ks_2samp(rep1_rpkm,rep3_rpkm))
    print("RPKM_peak r2 v r3: ",scipy.stats.ks_2samp(rep2_rpkm,rep3_rpkm))

    print("Read_count r1 v r2: ",scipy.stats.ks_2samp(rep1_read_count,rep2_read_count))
    print("Read_count r1 v r3: ",scipy.stats.ks_2samp(rep1_read_count,rep3_read_count))
    print("Read_count r2 v r3: ",scipy.stats.ks_2samp(rep2_read_count,rep3_read_count))

    print("Peak_height r1 v r2: ",scipy.stats.ks_2samp(rep1_peak_height,rep2_peak_height))
    print("Peak_height r1 v r3: ",scipy.stats.ks_2samp(rep1_peak_height,rep3_peak_height))
    print("Peak_height r2 v r3: ",scipy.stats.ks_2samp(rep2_peak_height,rep3_peak_height))

    print("TPM_peak r1 v r2: ",scipy.stats.ks_2samp(rep1_peak_tpm,rep2_peak_tpm))
    print("TPM_peak r1 v r3: ",scipy.stats.ks_2samp(rep1_peak_tpm,rep3_peak_tpm))
    print("TPM_peak r2 v r3: ",scipy.stats.ks_2samp(rep2_peak_tpm,rep3_peak_tpm))

    print("--------------------------------------------------------------------------------------------------------")
    print("RPKM_peak r1 v r2: ",scipy.stats.kruskal(rep1_rpkm,rep2_rpkm))
    print("RPKM_peak r1 v r3: ",scipy.stats.kruskal(rep1_rpkm,rep3_rpkm))
    print("RPKM_peak r2 v r3: ",scipy.stats.kruskal(rep2_rpkm,rep3_rpkm))

    print("Read_count r1 v r2: ",scipy.stats.kruskal(rep1_read_count,rep2_read_count))
    print("Read_count r1 v r3: ",scipy.stats.kruskal(rep1_read_count,rep3_read_count))
    print("Read_count r2 v r3: ",scipy.stats.kruskal(rep2_read_count,rep3_read_count))

    print("Peak_height r1 v r2: ",scipy.stats.kruskal(rep1_peak_height,rep2_peak_height))
    print("Peak_height r1 v r3: ",scipy.stats.kruskal(rep1_peak_height,rep3_peak_height))
    print("Peak_height r2 v r3: ",scipy.stats.kruskal(rep2_peak_height,rep3_peak_height))

    print("TPM_peak r1 v r2: ",scipy.stats.kruskal(rep1_peak_tpm,rep2_peak_tpm))
    print("TPM_peak r1 v r3: ",scipy.stats.kruskal(rep1_peak_tpm,rep3_peak_tpm))
    print("TPM_peak r2 v r3: ",scipy.stats.kruskal(rep2_peak_tpm,rep3_peak_tpm))

    print("--------------------------------------------------------------------------------------------------------")
    with open(os.path.join(args.output_path, "test_table.xlsx"), "w") as f:
        new_df.to_csv(f, sep="\t", index=None)
    with open(os.path.join(args.output_path, "exp11and16_table_validated.xlsx"), "w") as f:
        new_df.loc[new_df["Identifier"].isin(validated_orfs)].to_csv(f, sep="\t", index=None)

if __name__ == '__main__':
    main()
