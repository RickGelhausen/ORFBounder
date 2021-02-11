#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd

import collections
import csv

import operator

import lib.io
import lib.misc
import lib.predictions


def prediction_call(annotation_file, start_codons, stop_codons, fwd_wig_file, rev_wig_file, output_path, output_basename, p_offset, method, read_count_threshold):
    """
    execute the script for either TIS or TTS
    """
    if method == "TIS":
        search_codons = start_codons
        match_codons = stop_codons
    else:
        search_codons = stop_codons
        match_codons = start_codons

    fwd_wig_dict = io.load_wig(fwd_wig_file)
    rev_wig_dict = io.load_wig(rev_wig_file)

    print("Checking output folder...")
    if os.path.isdir(os.path.join(output_path, method)):
        sys.exit("Result directory already found! Please ensure that prior output folders with the same name are deleted.")

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

        annotation_fwd_interlap, annotation_rev_interlap, gene_density_dict, a_codon_pos = misc.annotation_interlap(annotation_file, method)
        gene_density_dict = predictions.calculate_density(fwd_wig_dict[key], annotation_fwd_interlap, gene_density_dict)
        gene_density_dict = predictions.calculate_density(rev_wig_dict[key], annotation_rev_interlap, gene_density_dict)

        fwd_codon_interlap, rev_codon_interlap, codon_dict = misc.create_codon_interlaps(args, key, val[0], search_codons)
        codon_dict = predictions.screen_wig_for_tts(fwd_wig_dict[key], fwd_codon_interlap, codon_dict, read_count_threshold)
        codon_dict = predictions.screen_wig_for_tts(rev_wig_dict[key], rev_codon_interlap, codon_dict, read_count_threshold)


        io.write_codon_interval_gff(output_path, "codon_intervals/%s_%s_intervals.gff" % (output_basename, key), codon_dict, p_offset, method)

        df_results, detected_codons_list = predictions.detect_potential_CDS(args, codon_dict, gene_density_dict, val[0], match_codons, p_offset, a_codon_pos, method)
        io.write_results_to_output_files(df_results, output_path, output_basename, split_gff, method)

    return detected_codons_list

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("--fwd_file_TIS", action="store", dest="fwd_wig_file_TIS", default="", help="input forward wig file for TIS.")
    parser.add_argument("--rev_file_TIS", action="store", dest="rev_wig_file_TIS", default="", help="input reverse wig file for TIS.")
    parser.add_argument("--fwd_file_TTS", action="store", dest="fwd_wig_file_TTS", default="", help="input forward wig file for TTS.")
    parser.add_argument("--rev_file_TTS", action="store", dest="rev_wig_file_TTS", default="", help="input reverse wig file for TTS.")

    parser.add_argument("-a", "--annotation_file", action="store", dest="annotation_file", help="input annotation file.", required=True)
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", help="input sequence file.", required=True)

    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])

    parser.add_argument("--offset_TIS", action="store", dest="p_offset_TIS", type=int, default=15)
    parser.add_argument("--offset_TTS", action="store", dest="p_offset_TTS", type=int, default=15)

    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", help="The max length to search for when using TIS and TTS.")
    parser.add_argument("--output_basename", action="store", dest="output_basename", default="result.gff", help="the basename for all output files." )
    parser.add_argument("-c", "--read_count_threshold", action="store", dest="read_count_threshold", default=5, type=int, help="skip reads lower than this threshold.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output path to the result folder.")
    args = parser.parse_args()

    method = handle_input(args)

    print("Fetching genome...")
    io.generate_genome_dict(genome_file)


    print("Done.")

    tis_predictions = []
    tts_predictions = []

    if method == "TIS":
        tis_predictions = prediction_call(args.annotation_file, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TIS, args.rev_wig_file_TIS, args.output_path, \
                                          args.output_basename, args.p_offset_TIS, "TIS", args.read_count_threshold)
    elif method == "TTS":
        tts_predictions = prediction_call(args.annotation_file, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TTS, args.rev_wig_file_TTS, args.output_path, \
                                          args.output_basename, args.p_offset_TTS, "TTS", args.read_count_threshold)

    else:
        tis_predictions = prediction_call(args.annotation_file, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TIS, args.rev_wig_file_TIS, args.output_path, \
                                          args.output_basename, args.p_offset_TIS, "TIS", args.read_count_threshold)

        tts_predictions = prediction_call(args.annotation_file, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TTS, args.rev_wig_file_TTS, args.output_path, \
                                          args.output_basename, args.p_offset_TTS, "TTS", args.read_count_threshold)


        detect_

    print("Terminating...")

if __name__ == '__main__':
    main()
