#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd

import collections
import csv

import operator

import lib.io as io
import lib.misc as misc
import lib.predictions as pred
import lib.expression as expr
import lib.messaging as msg


def prediction_call(annotation_file, genome_dict, start_codons, stop_codons, fwd_wig_file, rev_wig_file, output_path, \
                    output_basename, offset, method, longest_potential_ORF, detected_ORFs_dict, min_peak_height,
                    peak_height_calculation):
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
    if os.path.isfile(os.path.join(output_path, "result_tables", output_basename + ".csv")):
        msg.error("Error: Result table already found! Please ensure that prior output files with the same name are deleted.")

    msg.success("Done.")

    for key, val in genome_dict.items():
        msg.message("Current chromosome: %s" % key)
        if key not in fwd_wig_dict:
            msg.warning("Warning: No forward wig entry found for chrom: %s" % key)
            msg.warning("Skipping...")
            continue
        if key not in rev_wig_dict:
            msg.warning("Warning: No reverse wig entry found for chrom: %s" % key)
            msg.warning("Skipping...")
            continue

        codon_dict, gene_density_dict = {}, {}
        annotation_fwd_interlap, annotation_rev_interlap, gene_density_dict = misc.annotation_interlap(annotation_file, method)
        gene_density_dict = misc.calculate_density(fwd_wig_dict[key], annotation_fwd_interlap, gene_density_dict)
        gene_density_dict = misc.calculate_density(rev_wig_dict[key], annotation_rev_interlap, gene_density_dict)

        fwd_codon_interlap, rev_codon_interlap, codon_dict = misc.create_codon_interlaps(key, val, search_codons, offset)
        codon_dict = pred.screen_wig_for_tss(fwd_wig_dict[key], fwd_codon_interlap, codon_dict, min_peak_height, peak_height_calculation)
        codon_dict = pred.screen_wig_for_tss(rev_wig_dict[key], rev_codon_interlap, codon_dict, min_peak_height, peak_height_calculation)

        io.write_codon_interval_gff(output_path, os.path.join("codon_intervals","%s-%s_%s_intervals.gff" % (method, output_basename, key)), codon_dict, offset, method)

        fwd_area_interlap, rev_area_interlap, area_dict = misc.create_area_interlaps(key, val, search_codons, offset)
        area_dict = pred.screen_area_for_tss(fwd_wig_dict[key], fwd_area_interlap, area_dict)
        area_dict = pred.screen_area_for_tss(rev_wig_dict[key], rev_area_interlap, area_dict)

        #io.write_area_interval_gff(output_path, os.path.join("area_intervals","%s-%s_%s_intervals.gff" % (method, output_basename, key)), area_dict, offset, method)

        detected_ORFs_dict = pred.detect_potential_ORFs(codon_dict, val, search_codons, match_codons, offset, method, detected_ORFs_dict, longest_potential_ORF)

    return detected_ORFs_dict, gene_density_dict, codon_dict

def run_ORFBounder(fwd_wig_file_tis, rev_wig_file_tis, fwd_wig_file_TTS, rev_wig_file_TTS, bam_file_path, \
                annotation_file, genome_file, start_codons, stop_codons, output_path, output_basename, \
                offset_tis, offset_TTS, TTS_start_selection, min_peak_height, split_gff, max_ORF_length, \
                peak_height_calculation):
    """
    run functions necessary to generate the final output of ORFBounder
    """

    method = io.handle_input(fwd_wig_file_tis, rev_wig_file_tis, fwd_wig_file_TTS, rev_wig_file_TTS)
    bam_files = io.check_bamfile_input(bam_file_path, fwd_wig_file_tis, fwd_wig_file_TTS)

    headers = ["TIS","TTS"]
    if fwd_wig_file_tis != "":
        headers[0] = re.split('_|\.', os.path.basename(fwd_wig_file_tis))[0]
    if fwd_wig_file_TTS != "":
        headers[1] = re.split('_|\.', os.path.basename(fwd_wig_file_TTS))[0]

    wildcards = []
    if bam_files == -1:
        msg.message("No valid bam files detected, skipping readcount calculation")
    else:
        for file in bam_files:
            wildcards.append(re.split('_|\.', os.path.basename(file))[0])

        wildcards, bam_files = (list(t) for t in zip(*sorted(zip(wildcards, bam_files))))

    msg.message("Fetching genome...")
    genome_dict = io.generate_genome_dict(genome_file)
    msg.success("Done.")

    read_count_dict, total_mapped_list = {}, []
    predictions, gene_density_dict_tis, gene_density_dict_TTS = {}, {}, {}
    combined_result_df = pd.DataFrame()
    if method == "TIS":
        predictions, gene_density_dict_tis, _ \
                        = prediction_call(annotation_file, genome_dict, start_codons, stop_codons, \
                                          fwd_wig_file_tis, rev_wig_file_tis, output_path, \
                                          output_basename, offset_tis, "TIS", TTS_start_selection, \
                                          predictions, min_peak_height, peak_height_calculation)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, bam_files)

        result_df = misc.generate_result_dataframe(predictions, gene_density_dict_tis, {}, genome_dict, read_count_dict, \
                                            total_mapped_list, wildcards, method, headers)
        msg.success("Potential ORFs detected: %s" % len(result_df))
    elif method == "TTS":
        predictions, gene_density_dict_TTS, _ \
                        = prediction_call(annotation_file, genome_dict, start_codons, stop_codons, \
                                          fwd_wig_file_TTS, rev_wig_file_TTS, output_path, \
                                          output_basename, offset_TTS, "TTS", TTS_start_selection, \
                                          predictions, min_peak_height, peak_height_calculation)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, bam_files)

        result_df = misc.generate_result_dataframe(predictions, {}, gene_density_dict_TTS, genome_dict, read_count_dict, \
                                            total_mapped_list, wildcards, method, headers)

        msg.success("Potential ORFs detected: %s" % len(result_df))
    else:
        predictions, gene_density_dict_tis, codon_dict_tis \
                        = prediction_call(annotation_file, genome_dict, start_codons, stop_codons, \
                                          fwd_wig_file_tis, rev_wig_file_tis, output_path, \
                                          output_basename, offset_tis, "TIS", TTS_start_selection, \
                                          predictions, min_peak_height, peak_height_calculation)

        predictions, gene_density_dict_TTS, codon_dict_tts \
                        = prediction_call(annotation_file, genome_dict, start_codons, stop_codons, \
                                          fwd_wig_file_TTS, rev_wig_file_TTS, output_path, \
                                          output_basename, offset_TTS, "TTS", TTS_start_selection, \
                                          predictions, min_peak_height, peak_height_calculation)

        #combined_predictions = pred.combined_data_detection(codon_dict_tis, codon_dict_tts, offset_tis, offset_TTS, max_ORF_length)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict = expr.init_read_count_dict(read_count_dict, combined_predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, bam_files)

        result_df = misc.generate_result_dataframe(predictions, gene_density_dict_tis, gene_density_dict_TTS, genome_dict, \
                                            read_count_dict, total_mapped_list, wildcards, method, headers)

        #combined_result_df = misc.generate_result_dataframe(combined_predictions, gene_density_dict_tis, gene_density_dict_TTS, genome_dict, \
         #                                           read_count_dict, total_mapped_list, wildcards, method, headers)
        msg.success("Potential ORFs detected: %s" % len(result_df))
        #msg.success("Combined ORFs detected: %s" % len(combined_result_df))
        #
    return result_df, combined_result_df

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="ORFBounder is a peak detection and annotation script for TIS and TTS data. It can be run with either TIS, TTS or both.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--fwd_file_tis", action="store", dest="fwd_wig_file_tis", default="", help="input forward wig file for TIS.")
    parser.add_argument("--rev_file_tis", action="store", dest="rev_wig_file_tis", default="", help="input reverse wig file for TIS.")
    parser.add_argument("--fwd_file_TTS", action="store", dest="fwd_wig_file_TTS", default="", help="input forward wig file for TTS.")
    parser.add_argument("--rev_file_TTS", action="store", dest="rev_wig_file_TTS", default="", help="input reverse wig file for TTS.")

    parser.add_argument("-a", "--annotation_file", action="store", dest="annotation_file", help="input annotation file.", required=True)
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", help="input sequence file.", required=True)

    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])

    parser.add_argument("--offset_tis", action="store", dest="offset_tis", default=15)
    parser.add_argument("--offset_TTS", action="store", dest="offset_TTS", default=15)

    parser.add_argument("--peak_height_calculation", action="store", dest="peak_height_calculation", default="max"
                                                   , help="{max,sum}:\n"\
                                                         +"'max': within the codon interval select the highest value (> min_peak_height)"\
                                                         +"'sum': within the codon interval sum all values (> min_peak_height)")
    parser.add_argument("--TTS_start_selection", action="store", dest="TTS_start_selection", default="furthest_inframe"\
                                               , help="{furthest_inframe, next_inframe}\n"\
                                                      "'furthest_inframe': select the furthest inframe start codon that, without overstepping the next inframe stop codon.\n"\
                                                      "'next_inframe': select the closest inframe start codon.")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=100\
                                          , help="The maximum ORF length to take into account when using the combination method for TIS+TTS.")
    parser.add_argument("--min_peak_height", action="store", dest="min_peak_height", default=5, type=int\
                                           , help="Minimum height value to be considered a peak. (max option)\n"\
                                                 +"Minimum height value to be added to the total peak value (sum option)")
    parser.add_argument("--bam_file_path", action="store", dest="bam_file_path", default=""\
                                         , help="(optional) bam file to calculate RPKM and TE values for the final results.")
    parser.add_argument("--output_basename", action="store", dest="output_basename", required=True\
                                           , help="the basename for all output files." )
    parser.add_argument("-o","--output_path", action="store", dest="output_path", required=True\
                                            , help="Output path to the result folder.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    args = parser.parse_args()

    result_df, combined_result_df = \
             run_ORFBounder(args.fwd_wig_file_tis, args.rev_wig_file_tis, args.fwd_wig_file_TTS, args.rev_wig_file_TTS, \
                            args.bam_file_path, args.annotation_file, args.genome_file, args.start_codons, args.stop_codons, \
                            args.output_path, args.output_basename, args.offset_tis, args.offset_TTS, args.TTS_start_selection, \
                            args.min_peak_height, args.split_gff, args.max_ORF_length, args.peak_height_calculation)

    io.write_results_to_gff(result_df, os.path.join(args.output_path, "result_tables"), args.output_basename, args.split_gff)
    io.write_results_to_table(result_df, os.path.join(args.output_path, "result_tables"), args.output_basename)

    if not combined_result_df.empty:
        io.write_results_to_gff(combined_result_df, os.path.join(args.output_path, "combined_results"), args.output_basename, args.split_gff)
        io.write_results_to_table(combined_result_df, os.path.join(args.output_path, "combined_results"), args.output_basename)
    msg.success("Success! Terminating...")

if __name__ == '__main__':
    main()
