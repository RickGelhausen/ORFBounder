#!/usr/bin/env python
# TTS_analysis/detect_longest_potential_ORFs.py
import pandas as pd
import csv
import sys, os
import collections
import argparse
import re

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def excel_writer(args, data_frames):
    """
    create an excel sheet out of a dictionary of data_frames
    correct the width of each column
    """
    header_only =  ["Aminoacid_seq", "Nucleotide_seq"]
    writer = pd.ExcelWriter(args.out_xlsx, engine='xlsxwriter')
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

def generate_gene_dict(args):
    """
    create an interlap object for the annotation
    """
    annotation_dict = generate_annotation_dict(args.annotation_file)

    gene_dict = {}

    for key, val in annotation_dict.items():
        genome, mid, strand = key.split(":")
        start, stop = mid.split("-")
        start, stop = int(start), int(stop)

        if val[1] != "":
            locus_tag = val[1]
        else:
            locus_tag = val[4]

        if locus_tag not in gene_dict:
            gene_dict[locus_tag] = [genome, start, stop, strand]

    return gene_dict

def get_gene_information(chrom, start_position, stop_position, strand, gene_dict):
    # val : genome, start, stop, strand
    start_position, stop_position = start_position+1, stop_position+1
    # if strand == "-":
    #     start_position, stop_position = stop_position, start_position

    type = "Unannotated"
    for key, val in gene_dict.items():

        if val[0] != chrom:
            continue

        if strand == "+":
            gene_start, gene_stop = val[1], val[2]

            if start_position == gene_start and stop_position == gene_stop:
                type="Annotated"
                break
            elif abs(start_position-gene_start)<10 and val[3]==strand:
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

            if start_position == gene_start and stop_position == gene_stop:
                type="Annotated"
                break
            elif abs(start_position-gene_start)<10 and val[3]==strand:
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

    return type

def get_additional_ORF_information(chrom, genome_seq, cur_stop, strand, start_codons, stop_codons, gene_dict):
    """
    Search for stop upstream of current stop.
    Then search for start downstream from new stop.
    All done in-frame.
    """
    reverse_start_codons = [str(Seq(codon).reverse_complement()) for codon in start_codons]
    reverse_stop_codons = [str(Seq(codon).reverse_complement()) for codon in stop_codons]

    upstream_stop = 0
    longest_start = 0
    longest_nt = ""
    upstream_nt_stop = ""
    upstream_nt_start = ""

    loop_counter = 0
    if strand == "+":
        cur_position = cur_stop - 2
        nt=genome_seq[cur_position:cur_position+3]
        stop_nt_seq=nt
        while nt not in stop_codons or loop_counter == 0:
            loop_counter+=1
            cur_position -= 3
            if cur_position < 0:
                upstream_stop = ""
                upstream_nt_stop = ""
                break;

            nt = genome_seq[cur_position:cur_position+3]
            stop_nt_seq = nt + stop_nt_seq

        if upstream_stop == 0:
            upstream_stop = cur_position
            upstream_nt_stop = genome_seq[cur_position-15:cur_position]

        nt=genome_seq[cur_position:cur_position+3]
        while nt not in start_codons:
            cur_position += 3
            if cur_position == cur_stop - 2:
                return None

            nt = genome_seq[cur_position:cur_position+3]

        longest_start = cur_position
        longest_nt = genome_seq[longest_start:cur_stop+1]
        upstream_nt_start = genome_seq[longest_start-15:longest_start]

    else:
        cur_position = cur_stop
        nt=genome_seq[cur_position:cur_position+3]
        stop_nt_seq=nt

        while nt not in reverse_stop_codons or loop_counter == 0:
            cur_position += 3
            loop_counter+=1
            if cur_position > len(genome_seq):
                upstream_stop = ""
                upstream_nt_stop = ""
                break;

            nt = genome_seq[cur_position:cur_position+3]
            stop_nt_seq += nt

        if upstream_stop == 0:
            upstream_stop = cur_position
            upstream_nt_stop = str(Seq(genome_seq[cur_position+2:cur_position+17]).reverse_complement())

        nt=genome_seq[cur_position:cur_position+3]
        while nt not in reverse_start_codons:
            cur_position -= 3
            if cur_position == cur_stop:
                return None

            nt = genome_seq[cur_position:cur_position+3]

        stop_nt_seq = str(Seq(stop_nt_seq).reverse_complement())
        longest_start = cur_position + 2
        longest_nt = str(Seq(genome_seq[cur_stop:longest_start+1]).reverse_complement())
        upstream_nt_start = str(Seq(genome_seq[longest_start:longest_start+15]).reverse_complement())

    longest_gene_type = get_gene_information(chrom, longest_start, cur_stop, strand, gene_dict)

    return longest_start, longest_nt, upstream_nt_start, upstream_stop, stop_nt_seq, upstream_nt_stop, longest_gene_type

def postprocess_excel_file(args, genome_dict):
    """
    Search for additional upstream stop codon and determine longest ORF
    """

    gene_dict = generate_gene_dict(args)

    xlsx_df = pd.read_excel(args.in_xlsx, sheet_name=None)["all"]

    header = list(xlsx_df.columns)
    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    dynamic_header_part1 = [x for x in header if "_peak_height" in x or "_log2FC" in x]
    dynamic_header_part2 = [x for x in header if "_relative_density" in x or "'distance" in x]

    new_header = ["Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Shortest_Gene_type",\
                  "Shortest_codon_count", "Shortest_start_codon", "Position_alternativ_start", \
                  "Longest_Gene_type", "Longest_codon_count", "Longest_start_codon", "Stop_codon"] + dynamic_header_part1 +\
                 ["Shortest_15nt_upstream", "Shortest_Nucleotide_seq", "Shortest_Aminoacid_seq", \
                  "Longest_15nt_upstream", "Longest_Nucleotide_seq", "Longest_Aminoacid_seq", \
                  "Upstream_stop_codon", "Upstream_stop", "Stop_to_stop_nucleotide_seq"] + dynamic_header_part2

    rows_short = []
    rows_long = []
    name_list = ["s%s" % str(x) for x in range(len(new_header))]
    nTuple = collections.namedtuple('Pandas', name_list)
    result_rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        identifier = getattr(row, "Identifier")
        genome_id = getattr(row, "Genome")
        short_start = int(getattr(row, "Start")) - 1
        main_stop = int(getattr(row, "Stop")) - 1
        strand = getattr(row, "Strand")
        locus_tag = getattr(row, "Locus_tag")
        shortest_gene_type = getattr(row, "Gene_type")
        codon_count = getattr(row, "Codon_count")
        start_codon = getattr(row, "Start_codon")
        stop_codon = getattr(row, "Stop_codon")
        nucleotide_seq = getattr(row, "Nucleotide_seq")
        aminoacid_seq = getattr(row, "Aminoacid_seq")

        if strand == "-":
            main_stop, short_start = short_start, main_stop

        ORF_information = get_additional_ORF_information(genome_id, genome_dict[genome_id][0], main_stop, strand, args.start_codons, args.stop_codons, gene_dict)
        if ORF_information == None:
            continue

        longest_start, longest_nt, upstream_nt_start, upstream_stop, stop_nt_seq, upstream_nt_stop, longest_gene_type = ORF_information

        if strand == "-":
            main_stop, short_start = short_start, main_stop


        result = [identifier, genome_id, short_start+1, main_stop+1, strand, locus_tag, shortest_gene_type, codon_count, start_codon, longest_start+1, longest_gene_type] + \
                 [int(len(longest_nt) / 3), longest_nt[0:3], stop_codon] + [getattr(row, "_%s" % x) for x in range(10, 10+len(dynamic_header_part1))] + \
                 [getattr(row, "_%s" % (10+len(dynamic_header_part1)))] + [nucleotide_seq, aminoacid_seq] + \
                 [upstream_nt_start, longest_nt, str(Seq(longest_nt, generic_dna).translate(table=11, to_stop=False))] + \
                 [upstream_nt_stop, upstream_stop, stop_nt_seq] + [getattr(row, "_%s" % x) for x in range(13 + len(dynamic_header_part1), 13 + len(dynamic_header_part1)+len(dynamic_header_part2))]

        result_rows.append(nTuple(*result))

        attribute_short = "ID=%s;Name=%s;Gene_type=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;" % ("%s:%s-%s:%s" % (genome_id, short_start+1, main_stop+1, strand), locus_tag, shortest_gene_type, start_codon, stop_codon, codon_count)
        cur_tuple_short = nTuple_gff(genome_id, "TTS_finder", "CDS", short_start+1, main_stop+1, ".", strand, ".", attribute_short)
        if strand == "+":
            attribute_long = "ID=%s;Name=%s;Gene_type=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;" % ("%s:%s-%s:%s" % (genome_id, longest_start+1, main_stop+1, strand), locus_tag, longest_gene_type, longest_nt[0:3], stop_codon, int(len(longest_nt) / 3))
            cur_tuple_long = nTuple_gff(genome_id, "TTS_finder", "CDS", longest_start+1, main_stop+1, ".", strand, ".", attribute_long)
        else:
            attribute_long = "ID=%s;Name=%s;Gene_type=%s;Start_codon=%s;Stop_codon=%s;Codon_count=%s;" % ("%s:%s-%s:%s" % (genome_id, main_stop+1, longest_start+1, strand), locus_tag, longest_gene_type, longest_nt[0:3], stop_codon, int(len(longest_nt) / 3))
            cur_tuple_long = nTuple_gff(genome_id, "TTS_finder", "CDS", main_stop+1, longest_start+1, ".", strand, ".", attribute_long)

        rows_short.append(cur_tuple_short)
        rows_long.append(cur_tuple_long)
    df_short = pd.DataFrame.from_records(rows_short, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    df_long = pd.DataFrame.from_records(rows_long, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    all_df = pd.DataFrame.from_records(result_rows, columns=new_header)
    dataframe_dict = { "all" : all_df }


    with open(args.out_xlsx.replace(".xlsx", "_short.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.out_xlsx.replace(".xlsx", "_short.gff"), "a") as f:
        df_short.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

    with open(args.out_xlsx.replace(".xlsx", "_long.gff"), "w") as f:
        f.write("##gff-version 3\n")
    with open(args.out_xlsx.replace(".xlsx", "_long.gff"), "a") as f:
        df_long.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

    excel_writer(args, dataframe_dict)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Post processing of the TTS xlsx file returned by the merge_TTS.py script. Check for the longest potential ORF using the upstream stop.')
    parser.add_argument("-i", "--input_xlsx", action="store", dest="in_xlsx", required=True, help="Input excel file.")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", required=True, help="Genome file.")
    parser.add_argument("-a", "--annotation_file", action="store", dest="annotation_file", required=True, help="Annotation file used to determine gene_type.")
    parser.add_argument("--target_site", action="store", dest="target_site", default="TTS", help="TTS")
    parser.add_argument("-o", "--output_xlsx", action="store", dest="out_xlsx", required=True, help= "Output excel file.")
    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])

    args = parser.parse_args()

    print("Fetching genome...")
    # read the genome file
    genome_file = SeqIO.parse(args.genome_file, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    print("Done.")

    postprocess_excel_file(args, genome_dict)



if __name__ == '__main__':
    main()
