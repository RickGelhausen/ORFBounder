#!/usr/bin/env python
import argparse
import re
import os, sys
import pandas as pd

import collections
import csv

import operator

import lib.io as io
import lib.misc as misc
import lib.predictions as predictions
import lib.expression as expr
import lib.messaging as msg


def prediction_call(annotation_file, genome_dict, start_codons, stop_codons, fwd_wig_file, rev_wig_file, output_path, \
                    output_basename, p_offset, method, longest_potential_ORF, read_count_threshold, detected_ORFs_dict={}):
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

    msg.message("Checking output folder...")
    if os.path.isfile(os.path.join(output_path, method, "result_tables", output_basename + ".csv")):
        msg.warning("Result table already found! Please ensure that prior output files with the same name are deleted.")
        sys.exit()

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

        annotation_fwd_interlap, annotation_rev_interlap, gene_density_dict = misc.annotation_interlap(annotation_file, method)
        gene_density_dict = misc.calculate_density(fwd_wig_dict[key], annotation_fwd_interlap, gene_density_dict)
        gene_density_dict = misc.calculate_density(rev_wig_dict[key], annotation_rev_interlap, gene_density_dict)

        fwd_codon_interlap, rev_codon_interlap, codon_dict = misc.create_codon_interlaps(key, val, search_codons, p_offset)
        codon_dict = predictions.screen_wig_for_tss(fwd_wig_dict[key], fwd_codon_interlap, codon_dict, read_count_threshold)
        codon_dict = predictions.screen_wig_for_tss(rev_wig_dict[key], rev_codon_interlap, codon_dict, read_count_threshold)

        io.write_codon_interval_gff(output_path, os.path.join("codon_intervals","%s_%s_intervals.gff" % (output_basename, key)), codon_dict, p_offset, method)

        detected_ORFs_dict = predictions.detect_potential_ORFs(codon_dict, val, search_codons, match_codons, p_offset, method, longest_potential_ORF, detected_ORFs_dict)

    return detected_ORFs_dict, gene_density_dict

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
    parser.add_argument("--use_longest_TTS_ORF", action="store_true", dest="use_longest_TTS_ORF", help="Use the furthest possible inframe start codon for each detected stop codon to form the longest possible ORF that contains only one inframe stop codon. \
                                                                                                        Default uses the first detected start codon and may result in very short ORFs.")
    parser.add_argument("--bam_file_path", action="store", dest="bam_file_path", default="", help="(optional) bam file to calculate RPKM and TE values for the final results.")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=100, help="The max length to search for when using TIS and TTS combined.")
    parser.add_argument("--output_basename", action="store", dest="output_basename", required=True, help="the basename for all output files." )
    parser.add_argument("-c", "--read_count_threshold", action="store", dest="read_count_threshold", default=5, type=int, help="skip reads lower than this threshold.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True, help="Output path to the result folder.")
    args = parser.parse_args()

    method = io.handle_input(args)
    bam_files = io.check_bamfile_input(args)
    wildcards = []
    if bam_files == -1:
        print("No valid bam files detected, skipping readcount calculation")
    else:
        for file in bam_files:
            wildcards.append(re.split('_|\.', os.path.basename(file))[0])

        wildcards, bam_files = (list(t) for t in zip(*sorted(zip(wildcards, bam_files))))

    print("Fetching genome...")
    genome_dict = io.generate_genome_dict(args.genome_file)
    print("Done.")

    read_count_dict, total_mapped_list = {}, []
    if method == "TIS":
        predictions, gene_density_dict_TIS \
                        = prediction_call(args.annotation_file, genome_dict, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TIS, args.rev_wig_file_TIS, args.output_path, \
                                          args.output_basename, args.p_offset_TIS, "TIS", args.use_longest_TTS_ORF, \
                                          args.read_count_threshold)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, wildcards, bam_files)

        io.write_results_to_output_files(predictions, gene_density_dict_TIS, {}, genome_dict, args.output_path, \
                                         args.output_basename, args.split_gff, read_count_dict, total_mapped_list, \
                                         wildcards, method)

    elif method == "TTS":
        predictions, gene_density_dict_TTS \
                        = prediction_call(args.annotation_file, genome_dict, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TTS, args.rev_wig_file_TTS, args.output_path, \
                                          args.output_basename, args.p_offset_TTS, "TTS", args.use_longest_TTS_ORF, \
                                          args.read_count_threshold)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, wildcards, bam_files)

        io.write_results_to_output_files(predictions, {}, gene_density_dict_TTS, genome_dict, args.output_path, \
                                         args.output_basename, args.split_gff, read_count_dict, total_mapped_list, \
                                         wildcards, method)


    else:
        predictions, gene_density_dict_TIS \
                        = prediction_call(args.annotation_file, genome_dict, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TIS, args.rev_wig_file_TIS, args.output_path, \
                                          args.output_basename, args.p_offset_TIS, "TIS", args.use_longest_TTS_ORF, \
                                          args.read_count_threshold)

        predictions, gene_density_dict_TTS \
                        = prediction_call(args.annotation_file, genome_dict, args.start_codons, args.stop_codons, \
                                          args.fwd_wig_file_TTS, args.rev_wig_file_TTS, args.output_path, \
                                          args.output_basename, args.p_offset_TTS, "TTS", args.use_longest_TTS_ORF, \
                                          args.read_count_threshold, predictions)

        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, wildcards, bam_files)

        io.write_results_to_output_files(predictions, gene_density_dict_TIS, gene_density_dict_TTS, \
                                         genome_dict, args.output_path, args.output_basename, args.split_gff, \
                                         read_count_dict, total_mapped_list, wildcards, method)

        #combined_predictions_dict = predictions.combined_data_detection(tis_predictions, tts_predictions, args.max_ORF_length)


    print("Terminating...")

if __name__ == '__main__':
    main()
