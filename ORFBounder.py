#!/usr/bin/env python
import argparse
import re
import os
import pandas as pd

from lib.alignment_reader import PositionReader

import lib.io as io
import lib.misc as misc
import lib.predictions as pred
import lib.expression as expr
import lib.messaging as msg


def prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping_mode, \
                    start_codons, stop_codons, alignment_file, output_path, output_basename, \
                    offset_dict, method, longest_potential_orf, detected_orfs_dict, min_peak_height, \
                    peak_height_operator, min_read_count_dict):
    """
    execute the script for either TIS or TTS
    """

    if method == "TIS" or method == "RIBO":
        search_codons = start_codons
        match_codons = stop_codons
    else:
        search_codons = stop_codons
        match_codons = start_codons

    msg.message("Checking output folder...")
    if os.path.isfile(os.path.join(output_path, "result_tables", output_basename + ".csv")):
        msg.error("Error: Result table already found! Please ensure that prior output files with the same name are deleted.")

    msg.success("Done.")

    pr_object = PositionReader(alignment_file, read_length_dict, mapping_mode, offset_dict)
    pr_object.normalize_read_counts(normalization, min_read_count_dict)
    alignment_position_dict, _ = pr_object.output()

    #pr_object.to_wig(output_path)

    for chrom, genome_seq in genome_dict.items():
        msg.message("Current chromosome: %s" % chrom)
        if (chrom, "+") not in alignment_position_dict and (chrom, "-") not in alignment_position_dict:
            msg.warning("Warning: No valid entry found for chrom: %s" % chrom)
            msg.warning("Skipping...")
            continue

        codon_dict, gene_density_dict = {}, {}
        annotation_interlap_dict, gene_density_dict = misc.annotation_interlap(annotation_file)
        gene_density_dict = misc.calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        codon_interlap_dict, codon_dict = misc.create_codon_interlaps(chrom, genome_seq, search_codons)
        codon_dict = pred.screen_positions_for_tss(alignment_position_dict, codon_interlap_dict, codon_dict, min_peak_height, peak_height_operator)

        #io.write_codon_interval_gff(output_path, output_basename+f"_{method}_{chrom}_codon_intervals.gff", codon_dict)

        detected_orfs_dict = pred.detect_potential_orfs(codon_dict, genome_seq, search_codons, match_codons, method, detected_orfs_dict, longest_potential_orf)

    return detected_orfs_dict, gene_density_dict, codon_dict

def run_orfbounder(alignment_file_tis, alignment_file_tts, alignment_file_ribo, read_length_json, normalization, mapping, \
                annotation_file, genome_file, start_codons, stop_codons, output_path, output_basename, \
                offset_json, tts_start_selection, min_peak_height, max_ORF_length, \
                peak_height_operator, all_reads_rpkm, alignment_file_path, total_read_file_path):
    """
    run functions necessary to generate the final output of ORFBounder
    """

    method = io.parse_alignment_input(alignment_file_tis, alignment_file_tts)
    bam_files = io.check_alignment_path_input(alignment_file_path, alignment_file_tis, alignment_file_tts)
    read_length_dict = io.parse_read_lengths(read_length_json)
    min_read_count_dict = io.parse_total_reads(total_read_file_path, normalization)

    headers = ["TIS","TTS","RIBO"]
    if alignment_file_tis != "":
        headers[0] = re.split('_|\.', os.path.basename(alignment_file_tis))[0]
    if alignment_file_tts != "":
        headers[1] = re.split('_|\.', os.path.basename(alignment_file_tts))[0]
    if alignment_file_ribo != "":
        headers[2] = re.split('_|\.', os.path.basename(alignment_file_ribo))[0]

    offset_dict = io.parse_offset_json(offset_json)

    wildcards = []
    if bam_files == -1:
        msg.message("No valid bam files detected in the bam folder, skipping readcount calculation")
    else:
        for file in bam_files:
            wildcards.append(re.split('_|\.', os.path.basename(file))[0])

        wildcards, bam_files = (list(t) for t in zip(*sorted(zip(wildcards, bam_files))))

    msg.message("Fetching genome...")
    genome_dict = io.generate_genome_dict(genome_file)
    msg.success("Done.")

    read_count_dict, accepted_read_list = {}, []
    predictions, gene_density_tis_dict, gene_density_tts_dict, gene_density_ribo_dict = {}, {}, {}, {}
    combined_result_df = pd.DataFrame()
    if method == "TIS":
        predictions, gene_density_tis_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_tis, output_path, \
                                          output_basename, offset_dict, "TIS", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)

        if alignment_file_ribo != "":
            predictions, gene_density_ribo_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_ribo, output_path, \
                                          output_basename, offset_dict, "RIBO", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, accepted_read_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(predictions, gene_density_tis_dict, {}, gene_density_ribo_dict, genome_dict, read_count_dict, \
                                            accepted_read_list, wildcards, method, headers)
        msg.success("Potential ORFs detected: %s" % len(result_df))
    elif method == "TTS":
        predictions, gene_density_tts_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_tts, output_path, \
                                          output_basename, offset_dict, "TTS", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)

        if alignment_file_ribo != "":
            predictions, gene_density_ribo_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_ribo, output_path, \
                                          output_basename, offset_dict, "RIBO", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, accepted_read_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(predictions, {}, gene_density_tts_dict, {}, genome_dict, read_count_dict, \
                                            accepted_read_list, wildcards, method, headers)

        msg.success("Potential ORFs detected: %s" % len(result_df))
    else:
        predictions, gene_density_tis_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_tis, output_path, \
                                          output_basename, offset_dict, "TIS", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)

        if alignment_file_ribo != "":
            predictions, gene_density_ribo_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_ribo, output_path, \
                                          output_basename, offset_dict, "RIBO", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)

        predictions, gene_density_tts_dict, _ \
                        = prediction_call(annotation_file, genome_dict, read_length_dict, normalization, mapping, \
                                          start_codons, stop_codons, alignment_file_tts, output_path, \
                                          output_basename, offset_dict, "TTS", tts_start_selection, \
                                          predictions, min_peak_height, peak_height_operator, min_read_count_dict)

        #combined_predictions = pred.combined_data_detection(codon_dict_tis, codon_dict_tts, offset_tis, offset_tts, max_ORF_length)
        if bam_files != -1:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, total_mapped_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(predictions, gene_density_tis_dict, gene_density_tts_dict, gene_density_ribo_dict, genome_dict, \
                                            read_count_dict, total_mapped_list, wildcards, method, headers)

        #combined_result_df = misc.generate_result_dataframe(combined_predictions, gene_density_dict_tis, gene_density_dict_tts, genome_dict, \
         #                                           read_count_dict, total_mapped_list, wildcards, method, headers)
        msg.success("Potential ORFs detected: %s" % len(result_df))
        #msg.success("Combined ORFs detected: %s" % len(combined_result_df))
        #

    return result_df, combined_result_df

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="ORFBounder is a peak detection and annotation script for TIS and TTS data.\n"\
                                                +"It can be run with either TIS, TTS or both.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--alignment_file_tis", action="store", dest="alignment_file_tis", type=str, default="", help="input alignment file for TIS (sam/bam format).")
    parser.add_argument("--alignment_file_tts", action="store", dest="alignment_file_tts", type=str, default="", help="input alignment file for TTS (sam/bam format).")
    parser.add_argument("--alignment_file_ribo", action="store", dest="alignment_file_ribo", type=str, default="", help="input alignment file for RIBO (sam/bam format).")

    parser.add_argument("-l", "--read_length_json", action="store", dest="read_length_json", default="", help="JSON file containing read-length specifications per file.\n"\
                                                                                                             +"Ranges can be given using the - symbol (e.g. 15-20,31,33-35,39")
    parser.add_argument("-m", "--mapping_method", action="store", dest="mapping", type=str, default="threeprime", help="Read-Mapping to be used:\n"\
                                                                                                                      +"{threeprime, fiveprime, centered, global}\n"\
                                                                                                                      +"'threeprime': only threeprime end positions of each mapped read are used.\n"\
                                                                                                                      +"'fiveprime':  only fiveprime end positions of each mapped read are used.\n"\
                                                                                                                      +"'centered': the three middle nucleotide positions of each mapped read are used.\n"\
                                                                                                                      +"'global': all positions of each mapped read are used.")
    parser.add_argument("-n", "--normalization_method", action="store", dest="normalization", type=str
                                                      , help="Readcount-Normalization methods to be used:\n"\
                                                            +"{raw,min,mil}\n"\
                                                            +"'raw': unnormalized readcounts\n"\
                                                            +"'min': normalized by min #aligned reads / #aligned reads .\n"\
                                                            +"'mil': normalized by 1000000 / #aligned reads.")

    parser.add_argument("-a", "--annotation_file", action="store", dest="annotation_file", type=str, required=True\
                                                 , help="input annotation file.")
    parser.add_argument("-g", "--genome_file", action="store", dest="genome_file", type=str, required=True\
                                             , help="input sequence file.")

    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG","GTG","TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG","TAA","TGA"])

    parser.add_argument("--offset_json", action="store", dest="offset_json", type=str, default=""
                                       , help="A JSON file containing offsets for each file/read-length combination.\n"\
                                             +"Default value will be used for missing entries.")
    parser.add_argument("--peak_height_operator", action="store", dest="peak_height_operator", type=str, default="max"
                                                , help="{max,sum}:\n"\
                                                      +"'max': within the codon interval select the highest value (> min_peak_height)\n"\
                                                      +"'sum': within the codon interval sum up all values (> min_peak_height)")
    parser.add_argument("--tts_start_selection", action="store", dest="tts_start_selection", type=str, default="furthest_inframe"\
                                               , help="{furthest_inframe, next_inframe}\n"\
                                                      +"'furthest_inframe': select the furthest inframe start codon that, without overstepping the next inframe stop codon.\n"\
                                                      +"'next_inframe': select the closest inframe start codon.")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=150\
                                          , help="The maximum ORF length to take into account when using the combination method for TIS+TTS.")
    parser.add_argument("--min_peak_height", action="store", dest="min_peak_height", type=int, default=5\
                                           , help="Minimum height value to be considered a peak. (max option)\n"\
                                                 +"Minimum height value to be added to the total peak value (sum option)")
    parser.add_argument("--total_read_file_path", action="store", dest="total_read_file_path", default=-1
                                                , help="A tab seperated file containing the total read lengths for each sample.\n"\
                                                      +"This file is only necessary when using the min nomalization method.\n"\
                                                      +"We provide a script that creates these files for you.")
    parser.add_argument("--all_reads_rpkm", action="store_true", dest="all_reads_rpkm", type=bool
                                          , help="If set, all mapped reads will be used for the calculation of RPKM values.\n"\
                                                +"By default only mapped reads of the specified lengths will be used")
    parser.add_argument("--output_basename", action="store", dest="output_basename", type=str, required=True\
                                           , help="the basename for all output files.")
    parser.add_argument("--alignment_file_path", action="store", dest="alignment_file_path", type=str, default=""\
                                               , help="(optional) sam/bam files to calculate RPKM and TE values for the final results.")
    parser.add_argument("-o","--output_path", action="store", dest="output_path", type=str, required=True\
                                            , help="Output path to the result folder.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    args = parser.parse_args()

    result_df, _ = run_orfbounder(args.alignment_file_tis, args.alignment_file_tts, args.alignment_file_ribo, args.read_lengths, \
                        args.normalization, args.mapping, args.annotation_file, args.genome_file, args.start_codons, args.stop_codons, args.output_path, \
                        args.output_basename, args.offset_json, args.tts_start_selection, args.min_peak_height, \
                        args.max_ORF_length, args.peak_height_operator, args.all_reads_rpkm, args.alignment_file_path, args.total_read_file_path)

    io.write_results_to_gff(result_df, os.path.join(args.output_path, "result_tables"), args.output_basename, args.split_gff)
    io.write_results_to_table(result_df, os.path.join(args.output_path, "result_tables"), args.output_basename)

    # if not combined_result_df.empty:
    #     io.write_results_to_gff(combined_result_df, os.path.join(args.output_path, "combined_results"), args.output_basename, args.split_gff)
    #     io.write_results_to_table(combined_result_df, os.path.join(args.output_path, "combined_results"), args.output_basename)
    msg.success("Success! Terminating...")

if __name__ == '__main__':
    main()
