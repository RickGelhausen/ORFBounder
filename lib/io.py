#!/usr/bin/env python

import os, sys
import re
import csv
import collections
import pandas as pd
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import lib.misc as misc
import lib.expression as expr
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

def handle_input(args):
    """
    Check if input is valid.
    """

    if args.fwd_wig_file_TIS != "" and args.rev_wig_file_TIS != "" and args.fwd_wig_file_TTS != "" and args.rev_wig_file_TTS != "":
        return "combined_methods"

    if args.fwd_wig_file_TIS != "" and args.rev_wig_file_TIS != "":
        return "TIS"

    if args.fwd_wig_file_TTS != "" and args.rev_wig_file_TTS != "":
        return "TTS"

    sys.exit("Please ensure to either provide 2 TIS files, 2 TTS files OR both")

def check_bamfile_input(args):
    """
    Check bam input path.
    Ensure that there is:
     - one bam file corresponding to each input method (TIS, TTS)
     - (optional) one RNA bam file corresponding to each input method (TIS, TTS)
    """

    if args.bam_file_path == "":
        return -1

    valid_bam = set()

    _, _, bam_file_list = next(os.walk(args.bam_file_path))
    bam_file_list = [ file for file in bam_file_list if file.endswith(".bam")]
    if args.fwd_wig_file_TIS != "":
        TIS_prefix = os.path.basename(args.fwd_wig_file_TIS).split(".")[0]
        RNATIS_prefix = "RNATIS-" + "-".join(TIS_prefix.split("-")[1:])
        for file in bam_file_list:
            if TIS_prefix in file or RNATIS_prefix in file:
                valid_bam.add(os.path.join(args.bam_file_path,file))

    if args.fwd_wig_file_TTS != "":
        TTS_prefix = os.path.basename(args.fwd_wig_file_TTS).split(".")[0]
        RNATTS_prefix = "RNATTS-" + "-".join(TTS_prefix.split("-")[1:])
        for file in bam_file_list:
            if TTS_prefix in file or RNATTS_prefix in file:
                valid_bam.add(os.path.join(args.bam_file_path,file))

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
                    sys.exit("Incomplete header in wig file! Missing chrom= field!")

                wig_data_dict[chromosome].append(line)

            elif "chrom=" in line:
                tmp = re.split('[ =]', line)
                chromosome = tmp[tmp.index("chrom")+1]
                if chromosome not in wig_data_dict:
                    wig_data_dict[chromosome] = []

    return wig_data_dict

def write_gff_file(dataframe_out, output_path, output_filename, method):
    """
    write a dataframe to a gff file
    """
    filename = os.path.join(output_path, method, output_filename)
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

    write_gff_file(df, output_path, output_basename, method)

def excel_writer(out_file_name, data_frames, wildcards):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Aminoacid_seq", "Nucleotide_seq", "Start_codon", "Stop_codon", "Strand", "Codon_count"] + [card + "_rpkm" for card in wildcards]
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

def write_results_to_output_files(detected_ORFs_dict, gene_dict_TIS, gene_dict_TTS, genome, output_path, output_basename, \
                                    split_gff, read_count_dict, total_mapped_list, wildcards, method):
    """
    write a csv file containing all information,
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    TIS_header = "TIS"
    TTS_header = "TTS"
    for card in wildcards:
        if "TIS" in card and not "RNA" in card:
            TIS_header = card

        if "TTS" in card and not "RNA" in card:
            TTS_header = card

    TE_header = expr.get_TE_header(wildcards)
    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count", \
              TIS_header + "_peak_height", TTS_header + "_peak_height", "Start_codon", "Stop_codon", "15nt_window",\
              "Nucleotide_Seq", "Amino_Acid_Seq", TIS_header + "_relative_density", TTS_header + "_relative_density", \
              "5'-distance", "3'-distance"] + [card + "_rpkm" for card in wildcards] +\
              [cond + "_TE" for cond in TE_header]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    result_rows = []
    for (chrom, strand) in detected_ORFs_dict.keys():
        for (start, stop) in detected_ORFs_dict[(chrom, strand)].keys():
            rpm_start, rpm_stop = detected_ORFs_dict[(chrom, strand)][(start, stop)]
            if gene_dict_TIS != {}:
                gene_type, gene_name = misc.get_gene_information(chrom, start, stop, strand, gene_dict_TIS)
            else:
                gene_type, gene_name = misc.get_gene_information(chrom, start, stop, strand, gene_dict_TTS)
            nt_seq, aa_seq, nt_window, start_codon, stop_codon = misc.get_genome_information(start, stop, strand, genome[chrom], method)

            if aa_seq.count("*") > 1:
                continue
            rpkm_list = []
            TE_list = []
            if read_count_dict != {}:
                for idx, val in enumerate(read_count_dict[(chrom, start, stop, strand)]):
                    rpkm_list.append(expr.calculate_rpkm(total_mapped_list[idx][chrom], val, len(nt_seq)))

                TE_list = expr.calculate_TE(rpkm_list, wildcards)

            identifier = "%s:%s-%s:%s" % (chrom, start+1, stop+1, strand)
            codon_count = int(len(nt_seq)/3)

            if gene_dict_TIS != {}:
                fiveprime_dist, threeprime_dist = misc.calculate_utr_distance(start, stop, gene_name, gene_dict_TIS, method)
                relative_density_start = misc.calculate_relative_density(rpm_start, gene_name, gene_type, gene_dict_TIS)
                relative_density_stop = misc.calculate_relative_density(rpm_stop, gene_name, gene_type, gene_dict_TTS)
            else:
                fiveprime_dist, threeprime_dist = misc.calculate_utr_distance(start, stop, gene_name, gene_dict_TTS, method)
                relative_density_start = misc.calculate_relative_density(rpm_start, gene_name, gene_type, gene_dict_TIS)
                relative_density_stop = misc.calculate_relative_density(rpm_stop, gene_name, gene_type, gene_dict_TTS)

            rpm_start = rpm_start if rpm_start != -1 else "NaN"
            rpm_stop = rpm_stop if rpm_stop != -1 else "NaN"

            attribute = "ID=%s;Name=%s;Peak_height_TIS=%s;Peak_height_TTS=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;Type=%s" \
                        % (identifier, gene_name, rpm_start, rpm_stop, start_codon, stop_codon, codon_count, gene_type)
            cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start)+1, int(stop)+1, ".", strand, ".", attribute)

            result = [gene_type, identifier, chrom, start+1, stop+1, strand, gene_name, codon_count, rpm_start, rpm_stop, \
                      start_codon, stop_codon, nt_window, nt_seq, aa_seq, relative_density_start, relative_density_stop, \
                      fiveprime_dist, threeprime_dist] + rpkm_list + TE_list

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

            result_rows.append(nTuple(*result))

    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])

    print("Generating gff files...")
    df_all = pd.DataFrame.from_records(gff_all, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    write_gff_file(df_all, output_path, os.path.join("results_gff","%s.gff" % output_basename), method)

    if split_gff:
        df_annotated = pd.DataFrame.from_records(gff_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_annotated, output_path, os.path.join("results_gff","%s_annotated.gff" % output_basename), method)

        df_unannotated = pd.DataFrame.from_records(gff_unannotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_unannotated, output_path, os.path.join("results_gff","%s_unannotated.gff" % output_basename), method)

        df_near_annotated = pd.DataFrame.from_records(gff_near_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_near_annotated, output_path, os.path.join("results_gff","%s_near_annotated.gff" % output_basename), method)

        df_internal_inframe = pd.DataFrame.from_records(gff_internal_inframe, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_internal_inframe, output_path, os.path.join("results_gff","%s_internal_inframe.gff" % output_basename), method)

        df_n_terminal = pd.DataFrame.from_records(gff_n_terminal, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_n_terminal, output_path, os.path.join("results_gff","%s_n_terminal.gff" % output_basename), method)

        df_internal_out = pd.DataFrame.from_records(gff_internal_out, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        write_gff_file(df_internal_out, output_path, os.path.join("results_gff","%s_internal_out.gff" % output_basename), method)

    print("Done.")

    print("Generating output_table...")
    out_csv = os.path.join(output_path, method, "result_tables", "%s.csv" % output_basename)
    Path(os.path.dirname(out_csv)).mkdir(parents=True, exist_ok=True)

    df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    out_xlsx = os.path.join(output_path, method, "result_tables", "%s.xlsx" % output_basename)
    df_dict = {"CDS" : df_results}
    Path(os.path.dirname(out_xlsx)).mkdir(parents=True, exist_ok=True)

    excel_writer(out_xlsx, df_dict, wildcards)
    print("Done.")
