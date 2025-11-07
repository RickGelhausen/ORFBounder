#!/usr/bin/env python

import os
import csv
import json
import collections
import pandas as pd

from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO


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

def parse_read_lengths(read_length_json):
    """
    Parse the read length input into a continuous list form.
    """

    if read_length_json == "" or read_length_json == -1:
        msg.warning("Warning: Empty read-lengths parameter given, using all available read lengths.")
        return -1

    if os.path.isfile(read_length_json):
        with open(read_length_json, 'r') as json_file:
            tmp_dict = json.load(json_file)
    else:
        msg.error("Error: Read-length JSON file does not exist! %s" % read_length_json)

    read_length_dict = {}
    for file, read_length_string in tmp_dict.items():
        if type(read_length_string) is int:
            read_length_dict[file] = [str(read_length_string)]

        elif type(read_length_string) is str:
            parts = read_length_string.split(",")
            read_lengths = set()
            for part in parts:
                if "-" in part:
                    interval = part.split("-")
                    i1, i2 = int(interval[0]), int(interval[1])
                    if i1 > i2:
                        i1, i2 = i2, i1
                    print(i1, i2)
                    for i in range(i1, i2+1):
                        read_lengths.add(i)
                else:
                    read_lengths.add(part)

            read_length_dict[file] = [str(x) for x in sorted(read_lengths)]
        else:
            msg.error("Error: Read-length JSON file is not in correct format!")

    print(read_length_dict)
    return read_length_dict

def parse_total_reads(mapped_counts_file_path, normalization_method):
    """
    Takes a tab seperated file of total read counts and determines the minimum for each chromosome
    Format:  sample chromosome total_reads
    """

    min_read_count_dict = {}
    if normalization_method == "min":
        error_msg = "Error: min normalization method chosen but no mapped_counts_file_path given!\n"\
                    "Either use a different normalization method or provide a file containing total read counts for each sample and each chromosome.\n"\
                    "Consider using our helper script to create the required files."

        if type(mapped_counts_file_path) is not str:
            msg.error(error_msg)
        elif not Path(mapped_counts_file_path).is_file():
            msg.error(error_msg)
        else:
            with open(mapped_counts_file_path, "r") as f:
                lines = list(filter(None, [line.strip() for line in f.readlines()]))

            for line in lines:
                sample, chrom, cur_count = line.split("\t")

                if chrom in min_read_count_dict:
                    if min_read_count_dict[chrom] > int(cur_count):
                        min_read_count_dict[chrom] = int(cur_count)
                else:
                    min_read_count_dict[chrom] = int(cur_count)

    else:
        return None

    return min_read_count_dict

def parse_alignment_input(alignment_file_tis, alignment_file_tts):
    """
    Check whether the input alignment files are valid and determine the execution method for ORFBounder
    """

    if alignment_file_tis != "" and alignment_file_tts != "":
        if not os.path.isfile(alignment_file_tis):
            msg.error("Error: Non-empty alignment file path given for TIS does not exist: %s " % alignment_file_tis)
        if not os.path.isfile(alignment_file_tts):
            msg.error("Error: Non-empty alignment file path given for TTS does not exist: %s " % alignment_file_tts)

        return "combined_methods"

    if alignment_file_tis != "":
        if not os.path.isfile(alignment_file_tis):
            msg.error("Error: Non-empty alignment file path given for TIS does not exist: %s " % alignment_file_tis)
        return "TIS"

    if alignment_file_tts != "":
        if not os.path.isfile(alignment_file_tts):
            msg.error("Error: Non-empty alignment file path given for TTS does not exist: %s " % alignment_file_tts)
        return "TTS"

    msg.error("Error: Please ensure to either provide a TIS file, a TTS file or both!")


def check_alignment_path_input(alignment_file_path, alignment_file_tis, alignment_file_tts):
    """
    Check alignment input path.
    Ensure that there is:
     - one sam/bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if alignment_file_path == "" or not isinstance(alignment_file_path, str):
        return -1

    valid_bam = set()

    _, _, alignment_file_list = next(os.walk(alignment_file_path))
    alignment_file_list = [ file for file in alignment_file_list if file.endswith(".bam") or file.endswith(".sam")]

    condition, replicate = "", ""
    if alignment_file_tis != "":
        tis_prefix = os.path.basename(alignment_file_tis).split(".")[0]
        condition, replicate = tis_prefix.split("-")[1:]
        rnatis_prefix = "RNATIS-%s-%s" % (condition, replicate)
        for file in alignment_file_list:
            if tis_prefix in file or rnatis_prefix in file:
                valid_bam.add(os.path.join(alignment_file_path, file))

    if alignment_file_tts != "":
        tts_prefix = os.path.basename(alignment_file_tts).split(".")[0]
        condition, replicate = tts_prefix.split("-")[1:]
        rnatts_prefix = "RNATTS-%s-%s" % (condition, replicate)
        for file in alignment_file_list:
            if tts_prefix in file or rnatts_prefix in file:
                valid_bam.add(os.path.join(alignment_file_path, file))

    if condition != "" and replicate != "":
        for file in alignment_file_list:
            if "RIBO-%s-%s" % (condition, replicate) in file or "RNA-%s-%s" % (condition, replicate) in file:
                valid_bam.add(os.path.join(alignment_file_path, file))

    if len(valid_bam) == 0:
        return -1

    return valid_bam

def parse_offset_json(offset_json):
    """
    Read offset JSON file into a dictionary
    """

    if os.path.isfile(offset_json):
        with open(offset_json, 'r') as json_file:
            offset_dict = json.load(json_file)
        #except:
        #    msg.error("Error: Provided Offset JSON file is not in correct JSON format!")
    else:
        msg.error("Error: Offset JSON file does not exist! %s" % offset_json)

    return offset_dict

def write_gff_file(dataframe_out, output_path, output_filename):
    """
    write a dataframe to a gff file
    """
    file_name = os.path.join(output_path, output_filename)
    Path(os.path.dirname(file_name)).mkdir(parents=True, exist_ok=True)

    msg.message("Writing: %s" % file_name)
    with open(file_name, "w") as f:
        f.write("##gff-version 3\n")
    with open(file_name, "a") as f:
        dataframe_out.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

def write_codon_interval_gff(output_path, output_basename, codon_dict):
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

        attribute = "ID=%s;Peak_height=%s;Name=%s;Start_codon=%s" % (key, val[1], val[0], val[0])

        rows.append(nTuple_gff(chrom, "ORFBounder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    write_gff_file(df, output_path, output_basename)

# def write_area_interval_gff(output_path, output_basename, area_dict, offset, method):
#     """
#     Create a gff3 file with all codon intervals.
#     """

#     nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

#     rows = []
#     for key, val in area_dict.items():
#         if val[1] <= 0:
#             continue
#         chrom, mid, strand = key.split(":")
#         start, stop = mid.split("-")

#         if method == "TIS":
#             if strand == "+":
#                 cur_position = int(start) - offset + 24
#             elif strand == "-":
#                 cur_position = int(start) + offset + 24

#             attribute = "ID=%s;Area_coverage=%s;Name=%s;Start_codon=%s;Original_position=%s" % ("%s:%s-%s:%s" % (chrom,int(start)+1, int(stop)+1, strand), val[1], val[0], val[0], cur_position)
#         else:# change here if interval changes
#             if strand == "+":
#                 cur_position = int(start) - offset + 24
#             elif strand == "-":
#                 cur_position = int(start) + offset + 24

#             attribute = "ID=%s;Area_coverage=%s;Name=%s;Stop_codon=%s;Original_position=%s" % ("%s:%s-%s:%s" % (chrom, int(start)+1, int(stop)+1, strand), val[1], val[0], val[0], cur_position)

#         rows.append(nTuple_gff(chrom, "ORFBounder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

#     df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

#     write_gff_file(df, output_path, output_basename)

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
        worksheet.freeze_panes(1, 0)
        for idx, col in enumerate(df):
            series = df[col]
            if col in header_only:
                max_len = len(str(series.name)) + 2
            else:
                max_len = max(( series.astype(str).str.len().max(), len(str(series.name)) )) + 1
            #print("Sheet: %s | col: %s | max_len: %s" % (sheetname, col, max_len))
            worksheet.set_column(idx, idx, max_len)
    writer.close()

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
        gene_type, identifier, chrom, start, stop, strand, locus_tag, codon_count, rpm_tis, rpm_tts, rpm_ribo, start_codon, stop_codon = row[0:13]

        attribute = f"ID={identifier};Name={locus_tag};Peak_height_TIS={rpm_tis};Peak_height_TTS={rpm_tts};Peak_height_RIBO={rpm_ribo};Start_codon={start_codon};Stop_codon={stop_codon};Codon_count={codon_count};Type={gene_type};"
        cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start), int(stop), ".", strand, ".", attribute)

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
