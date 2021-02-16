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

def generate_genome_dict(genome_file):
    """
    read a genome fasta file into a dictionary
    """
    genome_file = SeqIO.parse(genome_file, "fasta")
    genome_dict = {}
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

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

def write_results_to_output_files(detected_ORFs_dict, gene_dict, , output_path, output_basename, split_gff, method):
    """
    write a csv file containing all information,
    write a gff file comtaining the ORFs from the csv,

    if split_gff == True then write one additional gff file for each gene_type
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "locus_tag", "codon_count", "peak_height", "start_codon", "stop_codon", "15nt window", "nt_seq", "aa_seq", "relative_density", "5'-distance", "3'-distance"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    for (chrom, strand), (start, stop, rpm_start, rpm_stop) in detected_ORFs_dict.items():
        gene_type, gene_name = get_gene_information(chrom, start, stop, strand, gene_dict)


        identifier = getattr(row, "Identifier")
        start_codon = getattr(row, "start_codon")
        stop_codon = getattr(row, "stop_codon")
        peak_height = float(getattr(row, "peak_height"))
        aa_length = len(getattr(row, "aa_seq"))
        locus_tag = getattr(row, "locus_tag")

        chrom, mid, strand = identifier.split(":")
        start, stop = mid.split("-")

        attribute = "ID=%s;Name=%s;Peak_height=%s;Start_codon=%s;Stop_codon=%s;AA_length=%s;Type=%s" % (identifier, locus_tag, peak_height, start_codon, stop_codon, aa_length, gene_type)
        cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start)-1, int(stop)-1, ".", strand, ".", attribute)

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

def write_output_files():
    """

    """



            if (chrom, cur_position+1, strand) in a_codon_pos:
                cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                gene_type = "Annotated"
                nt_seq = genome_seq[cur_start-1:cur_stop]
                aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
            else:
            aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))


        if aa_seq.count("*") > 1:
            continue

        if gene_type != "Annotated":

        else:
            out_start, out_stop = cur_start, cur_stop

    if gene_name in gene_dict:
        gene_rpm = gene_dict[gene_name][4]
        if gene_rpm != 0:
            relative_density = rpm / gene_rpm
        else:
            relative_density = "NaN"

        if method == "TIS":
            fiveprime_dist = out_start - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - out_start
        else:
            fiveprime_dist = out_stop - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - out_stop
    else:
        fiveprime_dist, threeprime_dist, relative_density = "NaN", "NaN", "NaN"


        if gene_type=="N-terminal_extension":
            relative_density = "NaN"

        if method == "TIS":
            if strand == "+":
                nt_window = genome_seq[out_start-16:out_start-1]
            else:
                nt_window = str(Seq(genome_seq[out_stop:out_stop+15]).reverse_complement())
        else:
            if strand == "+":
                nt_window = genome_seq[out_stop:out_stop+15]
            else:
                nt_window = str(Seq(genome_seq[out_start-16:out_start-1]).reverse_complement())

        start_codon, stop_codon = nt_seq[:3], nt_seq[-3:]

            unique_id="%s:%s-%s:%s" % (chrom, out_start, out_stop, strand)
            result = [gene_type, unique_id, chrom, out_start, out_stop, strand, gene_name, int(len(nt_seq)/3), rpm, start_codon, stop_codon, nt_window, nt_seq, aa_seq, relative_density, fiveprime_dist, threeprime_dist]
    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])

                                attributes = "ID=%s:%s-%s:%s;Name=%s:%s-%s:%s;Start_codon=%s;Stop_codon=%s;Start_peak=%s;Stop_peak=%s" \
                                                % (chrom, start, stop, strand, stop, start, stop, strand, start_codon, stop_codon, start_rpm, stop_rpm)
                        combined_ORFs_gff.append((chrom, "ORFBounder", "CDS", start, stop, ".", strand, ".", attributes))
nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    return pd.DataFrame.from_records(combined_ORFs_gff, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
