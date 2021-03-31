#!/usr/bin/env python

import os, sys
import re
import csv
import collections
import pandas as pd
import numpy as np

from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import lib.misc as misc
import lib.messaging as msg

def generate_genome_dict(genome_file):
    """
    read a genome fasta file into a dictionary
    """
    genome_file = SeqIO.parse(genome_file, "fasta")
    genome_dict = {}
    for entry in genome_file:
        genome_dict[str(entry.id)] = str(entry.seq)

    return genome_dict

def handle_input(fwd_wig_file_TIS, rev_wig_file_TIS, fwd_wig_file_TTS, rev_wig_file_TTS):
    """
    Check if input is valid.
    """

    if fwd_wig_file_TIS != "" and rev_wig_file_TIS != "" and fwd_wig_file_TTS != "" and rev_wig_file_TTS != "":
        return "combined_methods"

    if fwd_wig_file_TIS != "" and rev_wig_file_TIS != "":
        return "TIS"

    if fwd_wig_file_TTS != "" and rev_wig_file_TTS != "":
        return "TTS"

    msg.error("Error: Please ensure to either provide 2 TIS files, 2 TTS files OR both!")
    sys.exit()

def check_bamfile_input(bam_file_path, fwd_wig_file_TIS, fwd_wig_file_TTS):
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
    if fwd_wig_file_TIS != "":
        TIS_prefix = os.path.basename(fwd_wig_file_TIS).split(".")[0]
        RNATIS_prefix = "RNATIS-" + "-".join(TIS_prefix.split("-")[1:])
        for file in bam_file_list:
            if TIS_prefix in file or RNATIS_prefix in file:
                valid_bam.add(os.path.join(bam_file_path,file))

    if fwd_wig_file_TTS != "":
        TTS_prefix = os.path.basename(fwd_wig_file_TTS).split(".")[0]
        RNATTS_prefix = "RNATTS-" + "-".join(TTS_prefix.split("-")[1:])
        for file in bam_file_list:
            if TTS_prefix in file or RNATTS_prefix in file:
                valid_bam.add(os.path.join(bam_file_path,file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def load_wig(wig_path):
    """
    load wig file into a dictionary
    """
    with open(wig_path, 'r') as wig_file:
        chromosome = ""
        wig_data_dict = {}
        for line in wig_file.readlines():
            line = line.rstrip()

            if line[0].isdigit() and line[0] != "0":
                if chromosome not in wig_data_dict.keys():
                    msg.error("Error: Incomplete header in wig file! Missing chrom= field!")
                    sys.exit()

                wig_data_dict[chromosome].append(line)

            elif "chrom=" in line:
                tmp = re.split('[ =]', line)
                chromosome = tmp[tmp.index("chrom")+1]
                if chromosome not in wig_data_dict:
                    wig_data_dict[chromosome] = []

    return wig_data_dict

def write_gff_file(dataframe_out, output_path, output_filename):
    """
    write a dataframe to a gff file
    """
    filename = os.path.join(output_path, output_filename)
    Path(os.path.dirname(filename)).mkdir(parents=True, exist_ok=True)

    with open(filename, "w") as f:
        f.write("##gff-version 3\n")
    with open(filename, "a") as f:
        dataframe_out.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

def write_codon_interval_gff(output_path, output_basename, codon_dict, p_offset, method):
    """
    Create a gff3 file with all codon intervals.
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue
        chrom, mid, strand = key.split(":")
        start, stop = mid.split("-")

        if method == "TIS":
            if strand == "+":
                cur_position = int(start) - p_offset + 2
            elif strand == "-":
                cur_position = int(start) + p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Start_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)
        else:# change here if interval changes
            if strand == "+":
                cur_position = int(start) - p_offset + 2
            elif strand == "-":
                cur_position = int(start) + p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Stop_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)

        rows.append(nTuple_gff(chrom, "ORFBounder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    write_gff_file(df, output_path, output_basename)

def excel_writer(out_file_name, data_frames):
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
            #print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.save()

def write_results_to_gff(result_df, output_path, output_basename, split_gff):
    """
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    for row in result_df.itertuples(index=False, name=None):
        gene_type, identifier, chrom, start, stop, strand, gene_name, codon_count, rpm_start, rpm_stop, start_codon, stop_codon = row[0:12]

        attribute = "ID=%s;Name=%s;Peak_height_TIS=%s;Peak_height_TTS=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;Type=%s" \
                    % (identifier, gene_name, rpm_start, rpm_stop, start_codon, stop_codon, codon_count, gene_type)
        cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start)+1, int(stop)+1, ".", strand, ".", attribute)

        gff_all.append(cur_tuple)
        if split_gff:
            if gene_type == "Annotated":
                gff_annotated.append(cur_tuple)
            elif gene_type == "Unannotated":
                gff_unannotated.append(cur_tuple)
            elif gene_type == "Near_Annotated":
                gff_near_annotated.append(cur_tuple)
            elif gene_type == "Internal_Inframe":
                gff_internal_inframe.append(cur_tuple)
            elif gene_type == "N-terminal_extension":
                gff_n_terminal.append(cur_tuple)
            elif gene_type == "Internal_OutofFrame":
                gff_internal_out.append(cur_tuple)

    msg.message("Generating gff files...")
    df_all = pd.DataFrame.from_records(gff_all, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    write_gff_file(df_all, output_path, os.path.join("result_gffs","%s.gff" % output_basename))

    if split_gff:
        df_annotated = pd.DataFrame.from_records(gff_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_annotated, output_path, os.path.join("result_gffs","%s_annotated.gff" % output_basename))

        df_unannotated = pd.DataFrame.from_records(gff_unannotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_unannotated, output_path, os.path.join("result_gffs","%s_unannotated.gff" % output_basename))

        df_near_annotated = pd.DataFrame.from_records(gff_near_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_near_annotated, output_path, os.path.join("result_gffs","%s_near_annotated.gff" % output_basename))

        df_internal_inframe = pd.DataFrame.from_records(gff_internal_inframe, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_internal_inframe, output_path, os.path.join("result_gffs","%s_internal_inframe.gff" % output_basename))

        df_n_terminal = pd.DataFrame.from_records(gff_n_terminal, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_n_terminal, output_path, os.path.join("result_gffs","%s_n_terminal.gff" % output_basename))

        df_internal_out = pd.DataFrame.from_records(gff_internal_out, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_internal_out, output_path, os.path.join("result_gffs","%s_internal_out.gff" % output_basename))
    msg.success("Done")

def write_results_to_table(df_results, output_path, output_basename):
    """
    write a csv and xlsx file containing all information,
    """

    msg.message("Generating output_tables...")
    out_csv = os.path.join(output_path, "%s.csv" % output_basename)
    Path(os.path.dirname(out_csv)).mkdir(parents=True, exist_ok=True)

    df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    out_xlsx = os.path.join(output_path, "%s.xlsx" % output_basename)
    df_dict = {"CDS" : df_results}
    Path(os.path.dirname(out_xlsx)).mkdir(parents=True, exist_ok=True)

    excel_writer(out_xlsx, df_dict)
    msg.success("Done")
