#!/usr/bin/env python

import pandas as pd
import csv
import sys, os
import collections
import argparse
import re

from interlap import InterLap

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


def annotation_interlap(args):
    """
    create an interlap object for the annotation
    """
    annotation_dict = generate_annotation_dict(args.annotation_gff)

    annotation_fwd_interlap = InterLap()
    annotation_rev_interlap = InterLap()
    gene_dict = {}

    for key, val in annotation_dict.items():
        genome, mid, strand = key.split(":")
        start, stop = mid.split("-")
        start, stop = int(start), int(stop)

        if val[1] != "":
            locus_tag = val[1]
        else:
            locus_tag = val[4]

        if args.target_site == "TIS":
            if args.direction == "both":
                if strand == "+":
                    annotation_fwd_interlap.add((start-args.size, start+args.size+2, locus_tag))
                else:
                    annotation_rev_interlap.add((stop-args.size-2, stop+args.size, locus_tag))
            elif args.direction == "up":
                if strand == "+":
                    annotation_fwd_interlap.add((start-args.size, start+2, locus_tag))
                else:
                    annotation_rev_interlap.add((stop-2, stop+args.size, locus_tag))
            elif args.direction == "down":
                if strand == "+":
                    annotation_fwd_interlap.add((start, start+args.size+2, locus_tag))
                else:
                    annotation_rev_interlap.add((stop-args.size-2, stop, locus_tag))

        else:
            if args.direction == "both":
                if strand == "+":
                    annotation_fwd_interlap.add((stop-args.size-2, stop+args.size, locus_tag))
                else:
                    annotation_rev_interlap.add((start-args.size, start+args.size+2, locus_tag))
            elif args.direction == "up":
                if strand == "+":
                    annotation_fwd_interlap.add((stop-args.size-2, stop, locus_tag))
                else:
                    annotation_rev_interlap.add((start, start+args.size+2, locus_tag))
            elif args.direction == "down":
                if strand == "+":
                    annotation_fwd_interlap.add((stop-2, stop+args.size, locus_tag))
                else:
                    annotation_rev_interlap.add((start-args.size, start+2, locus_tag))

    return annotation_fwd_interlap, annotation_rev_interlap

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

def postprocess_excel_file(args):
    annotation_fwd_interlap, annotation_rev_interlap = annotation_interlap(args)

    xlsx_df = pd.read_excel(args.in_xlsx, sheet_name=None)["all"]

    header = xlsx_df.columns
    rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        start = getattr(row, "Start")
        stop = getattr(row, "Stop")
        strand = getattr(row, "Strand")

        if args.target_site == "TIS":
            if strand == "+":
                matching_genes = list(annotation_fwd_interlap.find((start, start)))
            else:
                matching_genes = list(annotation_rev_interlap.find((stop, stop)))
        else:
            if strand == "+":
                matching_genes = list(annotation_fwd_interlap.find((stop, stop)))
            else:
                matching_genes = list(annotation_rev_interlap.find((start, start)))

        if len(matching_genes) == 0:
            rows.append(row)

    all_df = pd.DataFrame.from_records(rows, columns=header)
    dataframe_dict = { "all" : all_df }

    excel_writer(args, dataframe_dict)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Post processing of the TTS/TIS xlsx file returned by the merge_TTS.py script.')
    parser.add_argument("-i", "--input_xlsx", action="store", dest="in_xlsx", required=True, help= "Output excel file.")
    parser.add_argument("-a", "--annotation_gff", action="store", dest="annotation_gff", required=True, help= "Input annotation used for processing.")
    parser.add_argument("--target_site", action="store", dest="target_site", default="TTS", help="TTS / TIS")
    parser.add_argument("--direction", action="store", dest="direction", default="both", help="up/down/both (up/down of annotated stop)")
    parser.add_argument("--size", action="store", dest="size", default=25, type=int, help="size of the interval checked (if 'both' then 2*size)")
    parser.add_argument("-o", "--output_xlsx", action="store", dest="out_xlsx", required=True, help= "Output excel file.")

    args = parser.parse_args()

    postprocess_excel_file(args)



if __name__ == '__main__':
    main()
