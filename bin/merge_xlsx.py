#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd

import collections
import csv

def extend_combined_dictionary(xlsx_df, combined_dict):
    """
    collect data from current file and add it to the existing dictionary
    { chrom:start-stop:strand : { wildcard : peak_height, relative_density, RPKM, TE} }
    """

    for row in xlsx_df.itertuples(index=False, name='Pandas'):
        gene_type = getattr(row, "Type")
        unique_id = getattr(row, "Identifier")
        genome = getattr(row, "Genome")
        start = int(getattr(row, "Start"))
        stop = int(getattr(row, "Stop"))
        strand = getattr(row, "Strand")
        gene_name = getattr(row, "locus_tag")
        aa_count = int(getattr(row, "codon_count"))
        peak_height = float(getattr(row, "peak_height"))
        start_codon = getattr(row, "start_codon")
        stop_codon = getattr(row, "stop_codon")
        nt_upstream = getattr(row, "_11")
        nt_seq = getattr(row, "nt_seq")
        aa_seq = getattr(row, "aa_seq")
        relative_density = float(getattr(row, "relative_density"))
        fiveprime = getattr(row, "_15")
        threeprime = getattr(row, "_16")

def screen_input_tables(table_list):
    """
    screen over the input tables
    """

    combined_dict = {}
    for table in table_list:
        xlsx_df = pd.read_excel(table, sheet_name=None)["CDS"]

        combined_dict = extend_combined_dictionary(xlsx_df, combined_dict)



def merge_tables(table_list, output_path, output_basename):
    """
    collect information from all input tables and merge them into one final output table
    """

    screen_input_tables(table_list)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='Merge all results for different files into one xlsx file.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", required=True, help= "list of input tables.")
    parser.add_argument("--output_basename", action="store", dest="output_basename", required=True, help="the basename for all output files." )
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output path to the result folder.")
    args = parser.parse_args()

    merge_tables(table_list, output_path, output_basename)


if __name__ == '__main__':
    main()
