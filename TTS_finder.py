#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd

import collections
import csv
from interlap import InterLap
import operator

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def generate_annotation_dict(annotation_path):
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """

    annotation_df = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)
    annotation_dict = {}
    gene_dict = {}
    cds_dict = {}

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name")+1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list, old_locus_tag)
        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)

    for key in cds_dict.keys():
        gene_name = ""
        gene_id, locus_tag, name, read_list, old_locus_tag = cds_dict[key]

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_dict[key] = (gene_id, locus_tag, name, read_list, gene_name, old_locus_tag)

    return annotation_dict

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

def calculate_density(wig_file_data, annotation_interlap, gene_dict):
    """
    calculate density for each annotated gene
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))

        matching_genes = list(annotation_interlap.find((position, position)))
        for gene in matching_genes:
            gene_dict[gene[2]][4] += read_count

    return gene_dict

def annotation_interlap(args):
    """
    create an interlap object for the annotation
    """
    annotation_dict = generate_annotation_dict(args.annotation_file)

    annotation_fwd_interlap = InterLap()
    annotation_rev_interlap = InterLap()
    gene_dict = {}
    a_codon_pos = {}

    for key, val in annotation_dict.items():
        genome, mid, strand = key.split(":")
        start, stop = mid.split("-")
        start, stop = int(start), int(stop)

        if val[1] != "":
            locus_tag = val[1]
        else:
            locus_tag = val[4]

        if strand == "+":
            annotation_fwd_interlap.add((start, stop, locus_tag))
        else:
            annotation_rev_interlap.add((start, stop, locus_tag))

        if locus_tag not in gene_dict:
            gene_dict[locus_tag] = [genome, start, stop, strand, 0]

            if args.target_site == "TIS":
                if strand == "+":
                    a_codon_pos[(genome, int(start), strand)] = (start, stop, locus_tag)
                else:
                    a_codon_pos[(genome, int(stop)-2, strand)] = (start, stop, locus_tag)
            else:
                if strand == "+":
                    a_codon_pos[(genome, int(stop)-2, strand)] = (start, stop, locus_tag)
                else:
                    a_codon_pos[(genome, int(start), strand)] = (start, stop, locus_tag)

    return annotation_fwd_interlap, annotation_rev_interlap, gene_dict, a_codon_pos

def create_codon_interlaps(args, chrom, genome_seq, codons):
    """
    create interlaps around each codon, incorporating the offset
    """
    reverse_codons = [str(Seq(codon).reverse_complement()) for codon in codons]

    fwd_codon_interlap = InterLap()
    rev_codon_interlap = InterLap()
    codon_dict = {}
    for pos in range(len(genome_seq)-2):
        codon = genome_seq[pos:pos+3]
        # TODO MERGE TOGETHER ONCE IT IS CLEAR THAT IT STAYS LIKE THIS
        if codon in codons:
            if args.target_site == "TIS":
                interval_start = pos + args.p_offset - 2
                interval_stop = pos + args.p_offset + 2
            else:
                interval_start = pos + args.p_offset - 2
                interval_stop = pos + args.p_offset + 2
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "+")
            fwd_codon_interlap.add((interval_start, interval_stop, key))
            codon_dict[key] = [codon, 0]
        elif codon in reverse_codons:
            if args.target_site == "TIS":
                interval_start = pos - args.p_offset+2 - 2
                interval_stop = pos - args.p_offset+2 + 2
            else:
                interval_start = pos - args.p_offset+2 - 2
                interval_stop = pos - args.p_offset+2 + 2
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "-")
            rev_codon_interlap.add((interval_start, interval_stop, key))
            codon_dict[key] = [str(Seq(codon).reverse_complement()), 0]

    return fwd_codon_interlap, rev_codon_interlap, codon_dict

def screen_wig_for_tts(wig_file_data, codon_interlap, codon_dict, read_count_threshold):
    """
    screen over wig file and update the according codon entries
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))
        # if read_count > x here could be a readcount restriction
        if read_count <= read_count_threshold:# change here if interval changes
            continue
        matching_codons = list(codon_interlap.find((position, position)))
        for match in matching_codons:
            codon_dict[match[2]][1] += read_count
    return codon_dict

def get_gene_information(chrom, start_position, stop_position, strand, gene_dict):
    # val : genome, start, stop, strand, 0
    type = "Unannotated"
    for key, val in gene_dict.items():

        if val[0] != chrom:
            continue

        if strand == "+":
            gene_start, gene_stop = val[1], val[2]

            if abs(start_position-gene_start)<10 and val[3]==strand:
                type="Near_Annotated"
                break
            elif stop_position==gene_stop and val[3]==strand:
                if start_position > gene_start:
                    type="Internal_Inframe"
                    break
                else:
                    type="N-terminal_extension"
                    break
            else:
                if start_position >= gene_start and start_position <= gene_stop and val[3]==strand:
                    type="Internal_OutofFrame"
                    break
        else:
            gene_start, gene_stop = val[2], val[1]

            if abs(start_position-gene_start)<10 and val[3]==strand:
                type="Near_Annotated"
                break
            elif stop_position==gene_stop and val[3]==strand:
                if start_position < gene_start:
                    type="Internal_Inframe"
                    break
                else:
                    type="N-terminal_extension"
                    break
            else:
                if start_position <= gene_start and start_position >= gene_stop and val[3]==strand:
                    type="Internal_OutofFrame"
                    break

    if type == "Unannotated":
        key = "%s:%s-%s:%s" % (chrom, start_position, stop_position, strand)

    return type, key

def write_codon_interval_gff(args, codon_gff_path, codon_dict):
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

        if args.target_site == "TIS":
            if strand == "+":
                cur_position = int(start) - args.p_offset + 2
            elif strand == "-":
                cur_position = int(start) + args.p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Start_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)
        else:# change here if interval changes
            if strand == "+":
                cur_position = int(start) - args.p_offset + 2
            elif strand == "-":
                cur_position = int(start) + args.p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Stop_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)

        rows.append(nTuple_gff(chrom, "TTS_finder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    with open(codon_gff_path, "w") as f:
        f.write("##gff-version 3\n")
    with open(codon_gff_path, "a") as f:
        df.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

def prepare_output_file(args, codon_dict, gene_dict, genome_seq, match_codons, a_codon_pos):
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]

    rows_all = []
    rows_annotated = []
    rows_unannotated = []
    rows_near_annotated = []
    rows_internal_inframe = []
    rows_n_terminal = []
    rows_internal_out = []

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "locus_tag", "codon_count", "peak_height", "start_codon", "stop_codon", "15nt window", "nt_seq", "aa_seq", "relative_density", "5'-distance", "3'-distance"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, mid, strand = key.split(":")
        interval_start, interval_stop = mid.split("-")

        if args.target_site == "TIS":
            if strand == "+":
                cur_start = int(interval_start) - args.p_offset + 2 # change here if interval changes
                cur_position = cur_start

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = genome_seq[cur_start-1:cur_stop]
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in match_codons:
                        cur_position+=3
                        if cur_position > len(genome_seq)-2:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq += nt

                    if nt not in match_codons:
                        continue

                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_stop = cur_position + 2

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

            elif strand == "-":
                cur_start = int(interval_start) + 4 + args.p_offset - 2 # change here if interval changes
                cur_position = cur_start - 2

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = str(Seq(genome_seq[cur_start-1:cur_stop]).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in reverse_match_codons:
                        cur_position-=3
                        if cur_position < 0:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq = nt + nt_seq

                    if nt not in reverse_match_codons:
                        continue

                    nt_seq = str(Seq(nt_seq).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_stop = cur_position

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

        elif args.target_site == "TTS":
            if strand == "+":
                cur_stop = int(interval_start) - args.p_offset + 4 # change here if interval changes
                cur_position = cur_stop - 2

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = genome_seq[cur_start-1:cur_stop]
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in match_codons:
                        cur_position-=3
                        if cur_position < 0:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq = nt + nt_seq

                    if nt not in match_codons:
                        continue

                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_start = cur_position

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

            elif strand == "-":
                cur_stop = int(interval_start) + 2 + args.p_offset - 2 # change here if interval changes
                cur_position = cur_stop

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = str(Seq(genome_seq[cur_start-1:cur_stop]).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt# change here if interval changes
                    while nt not in reverse_match_codons:
                        cur_position+=3
                        if cur_position > len(genome_seq)-2:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq += nt

                    if nt not in reverse_match_codons:
                        continue

                    nt_seq = str(Seq(nt_seq).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_start = cur_position + 2

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)


        if aa_seq.count("*") > 1:
            continue

        if gene_type != "Annotated":
            if strand == "+":
                out_start, out_stop = cur_start+1, cur_stop+1
            else:
                out_start, out_stop = cur_stop+1, cur_start+1
        else:
            out_start, out_stop = cur_start, cur_stop

        rpm = val[1]
        if gene_name in gene_dict:
            gene_rpm = gene_dict[gene_name][4]
            if gene_rpm != 0:
                relative_density = rpm / gene_rpm
            else:
                relative_density = "NaN"

            if args.target_site == "TIS":
                fiveprime_dist = out_start - gene_dict[gene_name][1]
                threeprime_dist = gene_dict[gene_name][2] - out_start
            else:
                fiveprime_dist = out_stop - gene_dict[gene_name][1]
                threeprime_dist = gene_dict[gene_name][2] - out_stop
        else:
            fiveprime_dist, threeprime_dist, relative_density = "NaN", "NaN", "NaN"

        start_codon, stop_codon = nt_seq[:3], nt_seq[-3:]

        if gene_type=="N-terminal_extension":
            relative_density = "NaN"

        if args.target_site == "TIS":
            if strand == "+":
                nt_window = genome_seq[out_start-16:out_start-1]
            else:
                nt_window = str(Seq(genome_seq[out_stop:out_stop+15]).reverse_complement())
        else:
            if strand == "+":
                nt_window = genome_seq[out_stop:out_stop+15]
            else:
                nt_window = str(Seq(genome_seq[out_start-16:out_start-1]).reverse_complement())

        unique_id="%s:%s-%s:%s" % (chrom, out_start, out_stop, strand)
        result = [gene_type, unique_id, chrom, out_start, out_stop, strand, gene_name, int(len(nt_seq)/3), rpm, start_codon, stop_codon, nt_window, nt_seq, aa_seq, relative_density, fiveprime_dist, threeprime_dist]
        result_rows.append(nTuple(*result))

        attribute = "ID=%s;Name=%s;Peak_height=%s;Start_codon=%s;Stop_codon=%s;AA_length=%s;Type=%s" % (unique_id, gene_name, rpm, start_codon, stop_codon, int(len(nt_seq)/3), gene_type)
        cur_tuple = nTuple_gff(chrom, "TTS_finder", "CDS", out_start, out_stop, ".", strand, ".", attribute)
        rows_all.append(cur_tuple)
        if gene_type == "Annotated":
            rows_annotated.append(cur_tuple)
        elif gene_type == "Unannotated":
            rows_unannotated.append(cur_tuple)
        elif gene_type == "Near_Annotated":
            rows_near_annotated.append(cur_tuple)
        elif gene_type == "Internal_Inframe":
            rows_internal_inframe.append(cur_tuple)
        elif gene_type == "N-terminal_extension":
            rows_n_terminal.append(cur_tuple)
        elif gene_type == "Internal_OutofFrame":
            rows_internal_out.append(cur_tuple)

    df_all = pd.DataFrame.from_records(rows_all, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_annotated = pd.DataFrame.from_records(rows_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_unannotated = pd.DataFrame.from_records(rows_unannotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_near_annotated = pd.DataFrame.from_records(rows_near_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_internal_inframe = pd.DataFrame.from_records(rows_internal_inframe, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_n_terminal = pd.DataFrame.from_records(rows_n_terminal, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_internal_out = pd.DataFrame.from_records(rows_internal_out, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    print("Generating gff files...")
    with open(args.output_gff, "a") as f:
        df_all.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_annotated.gff"), "a") as f:
        df_annotated.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_unannotated.gff"), "a") as f:
        df_unannotated.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_near_annotated.gff"), "a") as f:
        df_near_annotated.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_internal_inframe.gff"), "a") as f:
        df_internal_inframe.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_internal_inframe.gff"), "a") as f:
        df_internal_inframe.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_n_terminal.gff"), "a") as f:
        df_n_terminal.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)
    with open(args.output_gff.replace(".gff", "_internal_out.gff"), "a") as f:
        df_internal_out.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

    print("Done.")
    print("Generating output_table...")
    output_df = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])
    if not os.path.isfile(args.output_file):
        output_df.to_csv(args.output_file, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    else:
        output_df.to_csv(args.output_file, sep="\t", index=False, quoting=csv.QUOTE_NONE, header=False, mode="a")
    print("Done.")

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("-f", "--fwd_file", action="store", dest="fwd_wig_file", help="input forward wig file.", required=True)
    parser.add_argument("-r", "--rev_file", action="store", dest="rev_wig_file", help="input reverse wig file.", required=True)
    parser.add_argument("-a", "--annotation_file", action="store", dest="annotation_file", help="input annotation file.", required=True)
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", help="input sequence file.", required=True)
    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])
    parser.add_argument("--p_offset", action="store", dest="p_offset", type=int, default=15)
    parser.add_argument("--output_gff", action="store", dest="output_gff", required=True, help="The gff output path.")
    parser.add_argument("--target_site", action="store", dest="target_site", default="TTS", help="TTS / TIS")
    parser.add_argument("-c", "--read_count_threshold", action="store", dest="read_count_threshold", default=5, type=int, help="skip reads lower than this threshold.")
    parser.add_argument("--codon_interval_out", action="store", dest="codon_interval_out", default="", help="If desired report all codon_intervals as gff3 file.")
    parser.add_argument("-o","--output_file", action="store", dest="output_file", required=True)
    args = parser.parse_args()

    if args.target_site == "TTS":
        search_codons = args.stop_codons
        match_codons = args.start_codons
    else:
        search_codons = args.start_codons
        match_codons = args.stop_codons

    print("Fetching genome...")
    # read the genome file
    genome_file = SeqIO.parse(args.genome_file, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    print("Done.")
    #print(gene_density_dict["Cj0004c"])
    #chrom = max(genome_dict.items(), key=operator.itemgetter(1))[0]
    #print("Chromosome: %s" % chrom)

    print("Preparing empty output files")
    with open(args.output_gff, "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_annotated.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_unannotated.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_near_annotated.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_internal_inframe.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_internal_inframe.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_n_terminal.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.output_gff.replace(".gff", "_internal_out.gff"), "w") as f:
        f.write("##gff-version 3\n")

    if os.path.isfile(args.output_file):
        sys.exit("File already found! Please ensure that prior output tables with the same name are deleted.")

    fwd_wig_dict = load_wig(args.fwd_wig_file)
    rev_wig_dict = load_wig(args.rev_wig_file)
    print("Computing predictions...")
    for key, val in genome_dict.items():
        print(key)
        if key not in fwd_wig_dict:
            print("No forward wig entry found for chrom: %s" % key)
            print("Skipping...")
            continue
        if key not in rev_wig_dict:
            print("No reverse wig entry found for chrom: %s" % key)
            print("Skipping...")
            continue

        annotation_fwd_interlap, annotation_rev_interlap, gene_density_dict, a_codon_pos = annotation_interlap(args)
        gene_density_dict = calculate_density(fwd_wig_dict[key], annotation_fwd_interlap, gene_density_dict)
        gene_density_dict = calculate_density(rev_wig_dict[key], annotation_rev_interlap, gene_density_dict)

        fwd_codon_interlap, rev_codon_interlap, codon_dict = create_codon_interlaps(args, key, val[0], search_codons)
        codon_dict = screen_wig_for_tts(fwd_wig_dict[key], fwd_codon_interlap, codon_dict, args.read_count_threshold)
        codon_dict = screen_wig_for_tts(rev_wig_dict[key], rev_codon_interlap, codon_dict, args.read_count_threshold)

        if args.codon_interval_out != "":
            write_codon_interval_gff(args, args.codon_interval_out.replace(".gff", "_%s.gff" % key), codon_dict)

        prepare_output_file(args, codon_dict, gene_density_dict, val[0], match_codons, a_codon_pos)

    print("Terminating...")

if __name__ == '__main__':
    main()
