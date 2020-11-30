#!/usr/bin/env python
# TTS_analysis/detect_longest_potential_ORFs.py
import pandas as pd
import csv
import sys, os
import collections
import argparse
import re



def xlsx_to_gff3(args):
    """
    Read excel file and create either one or two gff3 files depending on the target_site
    """

    target_site = ""
    xlsx_df = pd.read_excel(args.in_xlsx, sheet_name=None)["all"]


    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    header = list(xlsx_df.columns)
    if "Position_alternativ_start" in header:
        target_site = "TTS"
    else:
        target_site = "TIS"

    rows = []
    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        if target_site == "TTS":
            genome_id = getattr(row, "Genome")
            start = int(getattr(row, "Start"))
            stop = int(getattr(row, "Stop"))
            upstream_start = int(getattr(row, "Position_alternativ_start"))
            strand = getattr(row, "Strand")

            attribute_short = "ID=%s" % ("%s:%s-%s:%s" % (genome_id, start, stop, strand))
            cur_tuple_short = nTuple_gff(genome_id, "TTS_finder", "CDS", start, stop, ".", strand, ".", attribute_short)

            if strand == "+":
                attribute_long = "ID=%s" % ("%s:%s-%s:%s" % (genome_id, start, stop, strand))
                cur_tuple_long = nTuple_gff(genome_id, "TTS_finder", "CDS", upstream_start, stop, ".", strand, ".", attribute_long)
            else:
                attribute_long = "ID=%s" % ("%s:%s-%s:%s" % (genome_id, start, upstream_start, strand))
                cur_tuple_long = nTuple_gff(genome_id, "TTS_finder", "CDS", start, upstream_start, ".", strand, ".", attribute_long)

            rows.append(cur_tuple_short)
            if (strand == "+" and start != upstream_start) or (strand == "-" and stop != upstream_start):
                rows.append(cur_tuple_long)

        elif target_site == "TIS":
            genome_id = getattr(row, "Genome")
            start = int(getattr(row, "Start"))
            stop = int(getattr(row, "Stop"))
            strand = getattr(row, "Strand")

            attribute = "ID=%s" % ("%s:%s-%s:%s" % (genome_id, start, stop, strand))
            cur_tuple = nTuple_gff(genome_id, "TTS_finder", "CDS", start, stop, ".", strand, ".", attribute)

            rows.append(cur_tuple)

    df_gff = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    with open(args.out_gff, "w") as f:
        f.write("##gff-version 3\n")
    with open(args.out_gff, "a") as f:
        df_gff.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Transfer TTS/TIS xlsx file to gff3 for further anaylsis.')
    parser.add_argument("-i", "--input_xlsx", action="store", dest="in_xlsx", required=True, help="Input excel file.")
    parser.add_argument("-o", "--output_gff", action="store", dest="out_gff", required=True, help= "Output gff file.")

    args = parser.parse_args()

    xlsx_to_gff3(args)



if __name__ == '__main__':
    main()
