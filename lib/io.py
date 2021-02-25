#!/usr/bin/env python

import os
import re
import csv
import collections
import pandas as pd
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

import lib.misc as misc

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

def load_wig(wig_path):
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

def write_results_to_output_files(detected_ORFs_dict, gene_dict_TIS, gene_dict_TTS, genome, output_path, output_basename, split_gff, method):
    """
    write a csv file containing all information,
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count", \
              "Peak_height_TIS", "Peak_height_TTS", "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", \
              "Relative_density_start", "Relative_density_stop", "5'-distance", "3'-distance"]
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
        for (start, stop, rpm_start, rpm_stop) in detected_ORFs_dict[(chrom, strand)]:
            gene_type, gene_name = misc.get_gene_information(chrom, start, stop, strand, gene_dict_TIS)
            nt_seq, aa_seq, nt_window, start_codon, stop_codon = misc.get_genome_information(start, stop, strand, genome[chrom], method)

            if aa_seq.count("*") > 1:
                continue

            identifier = "%s:%s-%s:%s" % (chrom, start+1, stop+1, strand)
            codon_count = int(len(nt_seq)/3)

            fiveprime_dist, threeprime_dist = misc.calculate_utr_distance(start, stop, gene_name, gene_dict_TIS, method)
            relative_density_start = misc.calculate_relative_density(rpm_start, gene_name, gene_type, gene_dict_TIS)
            relative_density_stop = misc.calculate_relative_density(rpm_stop, gene_name, gene_type, gene_dict_TTS)

            rpm_start = rpm_start if rpm_start != -1 else "NaN"
            rpm_stop = rpm_stop if rpm_stop != -1 else "NaN"

            attribute = "ID=%s;Name=%s;Peak_height_TIS=%s;Peak_height_TTS=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;Type=%s" \
                        % (identifier, gene_name, rpm_start, rpm_stop, start_codon, stop_codon, codon_count, gene_type)
            cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start)+1, int(stop)+1, ".", strand, ".", attribute)

            result = [gene_type, identifier, chrom, start+1, stop+1, strand, gene_name, codon_count, rpm_start, rpm_stop, \
                      start_codon, stop_codon, nt_window, nt_seq, aa_seq, relative_density_start, relative_density_stop, \
                      fiveprime_dist, threeprime_dist]

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

    if not os.path.isfile(out_csv):
        df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    else:
        df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE, header=False, mode="a")
    print("Done.")
