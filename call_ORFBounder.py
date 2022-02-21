#!/usr/bin/env python
import os
import re
import argparse
import math

import pandas as pd
from pathlib import Path

import ORFBounder as ob
import lib.merging as mg
import lib.io as io
import lib.messaging as msg

def is_empty(entry):
    """
    Check if the current table entry is empty
    """

    if type(entry) == str and entry == "":
        print("str")
        return True
    elif type(entry) == float and math.isnan(entry):
        print("float")
        return True
    else:
        return False

def check_config_sheet(config_sheet):
    """
    Check whether the config sheet is correctly formatted
    """

    config_df = pd.read_csv(config_sheet, sep="\t")
    expected_columns = ["Annotation", "Bam_folder", "Experiment", "Filepath_RIBO", "Filepath_TIS", "Filepath_TTS", "Genome", "Normalization", "Offset_JSON", "Read_lengths", "Start_codons", "Stop_codons"]
    input_columns = config_df.columns

    if list(set(expected_columns) - set(input_columns)) != []:
        msg.error_list("Config Sheet columns are incomplete:\n",\
                       "Ensure that the file is TAB seperated.",\
                       "Required columns:",\
                       expected_columns,\
                       input_columns)

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        file_path_tis = getattr(row, "Filepath_TIS")
        file_path_tts = getattr(row, "Filepath_TTS")
        file_path_ribo = getattr(row, "Filepath_RIBO")
        read_lengths = getattr(row, "Read_lengths")
        mapping_method = getattr(row, "Mapping_method")
        normalization = getattr(row, "Normalization")
        offset_json = getattr(row, "Offset_JSON")
        bam_folder = getattr(row, "Bam_folder")

        with_tis = True
        with_tts = True
        with_ribo = True
        if is_empty(experiment):
            msg.error("Empty entry found: Missing Experiment!")
        if is_empty(annotation):
            msg.error("Empty entry found: Missing Annotation!")
        if is_empty(genome):
            msg.error("Empty entry found: Missing Genome!")
        if is_empty(file_path_tis):
            with_tis = False
        if is_empty(file_path_tts):
            with_tts = False
        if is_empty(file_path_ribo):
            with_ribo = False
        if is_empty(mapping_method):
            msg.error("Empty entry found: Missing Mapping_method!")
        if is_empty(normalization):
            msg.error("Empty entry found: Missing Normalization!")
        if is_empty(offset_json):
            msg.error("Empty entry found: Missing Offset_JSON!")
        if is_empty(read_lengths):
            msg.error("Empty entry found: Missing Read_lengths!")

        if not with_tis and not with_tts:
            msg.error("No TIS or TTS path given! Specify atleast one.")

        if not Path(annotation).is_file():
            msg.error("Annotation file is not valid! Ensure to enter a correct file path!\n%s" % annotation)
        if not Path(genome).is_file():
            msg.error("Genome file is not valid! Ensure to enter a correct file path!\n%s" % genome)
        if with_tis and not Path(file_path_tis).is_dir():
            msg.error("Mapping TIS directory is not valid! Ensure to enter a correct path!\n%s" % file_path_tis)
        if with_tts and not Path(file_path_tts).is_dir():
            msg.error("Mapping TTS directory is not valid! Ensure to enter a correct path!\n%s" % file_path_tts)
        if with_ribo and not Path(file_path_ribo).is_dir():
            msg.error("Mapping RIBO directory is not valid! Ensure to enter a correct path!\n%s" % file_path_ribo)
        if not Path(offset_json).is_file():
            msg.error("Offsets file is not valid! Ensure to enter a correct file path!\n%s" % offset_json)
        if bam_folder != "" and isinstance(bam_folder, str):
            if not Path(bam_folder).is_dir():
                msg.error("Given bam_folder is non-existant, either provide no bamfolder or an existing one!\n%s" % bamfolder)

        for norm in normalization.split(","):
            if norm not in ["raw", "mil", "min"]:
                msg.error("Given normalization method is not allowed: %s.\n Choose from {raw, mil, min}." % norm)

    return config_df

def retrieve_bam_input_information(file_path_tis, file_path_tts, file_path_ribo):
    """
    Create a list of matching TIS/TTS condition+replicate files to run together.
    """
    if file_path_tis != "":
        _, _, files_tis = next(os.walk(file_path_tis))
    else:
        files_tis = []

    if file_path_tts != "":
        _, _, files_tts = next(os.walk(file_path_tts))
    else:
        files_tts = []

    if file_path_ribo != "":
        _, _, files_ribo = next(os.walk(file_path_ribo))
    else:
        files_ribo = []

    tt_files = []
    tt_files.extend([bam for bam in files_tis if ("TIS" in bam and not "RNA" in bam) and (bam.endswith(".bam"))])
    tt_files.extend([bam for bam in files_tts if ("TTS" in bam and not "RNA" in bam) and (bam.endswith(".bam"))])
    tt_files.extend([bam for bam in files_ribo if ("RIBO" in bam and not "RNA" in bam) and (bam.endswith(".bam"))])

    sample_dict = {}
    for file in tt_files:
        wildcard = re.split('_|\.', os.path.basename(file))[0]
        method, condition, replicate = wildcard.split("-")

        if (condition, replicate) not in sample_dict:
            if (method == "TIS"):
                sample_dict[(condition, replicate)] = [os.path.join(file_path_tis, file),"",""]
            elif (method == "TTS"):
                sample_dict[(condition, replicate)] = ["",os.path.join(file_path_tts, file),""]
            elif (method == "RIBO"):
                sample_dict[(condition, replicate)] = ["","",os.path.join(file_path_ribo, file)]

        else:
            if (method == "TIS"):
                sample_dict[(condition, replicate)][0] = os.path.join(file_path_tis, file)
            elif (method == "TTS"):
                sample_dict[(condition, replicate)][1] = os.path.join(file_path_tts, file)
            elif (method == "RIBO"):
                sample_dict[(condition, replicate)][2] = os.path.join(file_path_ribo, file)


    bam_input_list = []
    for key, val in sample_dict.items():
        if val[0] != "" or val[1] != "" or  val[2] != "":

            bam_input_list.append((val, "%s-%s" % key))

    return bam_input_list

def call_ORFBounder(config_df, tts_start_selection, min_peak_height, peak_height_operator, max_orf_length, split_gff, result_path, contrasts, all_reads_rpkm):
    """
    Run the ORFBounder experiments specified in the config sheet.
    """

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        file_path_tis = getattr(row, "Filepath_TIS")
        file_path_tts = getattr(row, "Filepath_TTS")
        file_path_ribo = getattr(row, "Filepath_RIBO")
        read_lengths = getattr(row, "Read_lengths")
        mapping_method = getattr(row, "Mapping_method").split(",")
        normalization = getattr(row, "Normalization").split(",")
        offset_json = getattr(row, "Offset_JSON")
        start_codons = getattr(row, "Start_codons")
        stop_codons = getattr(row, "Stop_codons")
        bam_folder = getattr(row, "Bam_folder")

        if max_orf_length == "" or math.isnan(max_orf_length):
            max_orf_length = 150
        else:
            max_orf_length = int(max_orf_length)

        if start_codons == "" or type(start_codons) != str:
            start_codons = ["ATG", "GTG", "TTG"]
        else:
            start_codons = [codon.strip(" ") for codon in start_codons.split(",")]

        if stop_codons == "" or type(stop_codons) != str:
            stop_codons = ["TAG", "TAA", "TGA"]
        else:
            stop_codons = [codon.strip(" ") for codon in stop_codons.split(",")]

        for mapping in mapping_method:
            for norm in normalization:
                meta_dict, dynamic_dict = {}, {}
                combined_meta_dict, combined_dynamic_dict = {}, {}

                bam_input_list = retrieve_bam_input_information(file_path_tis, file_path_tts, file_path_ribo)

                res_path = os.path.join(result_path, experiment, mapping, norm)
                for (file_tis, file_tts, file_ribo), wildcard in bam_input_list:
                    try:
                        res_df, combined_res_df = ob.run_orfbounder(file_tis, file_tts, file_ribo, read_lengths, \
                                                                norm, mapping, annotation, genome, \
                                                                start_codons, stop_codons, res_path, wildcard, \
                                                                offset_json, tts_start_selection, min_peak_height, \
                                                                max_orf_length, peak_height_operator, \
                                                                all_reads_rpkm, bam_folder)
                    except SystemExit:
                        msg.warning("Error encountered while calling ORFBounder! Moving to next run!")
                        continue

                    io.write_results_to_gff(res_df, os.path.join(res_path, "result_tables"), wildcard, split_gff)
                    io.write_results_to_table(res_df, os.path.join(res_path, "result_tables"), wildcard)

                    meta_dict, dynamic_dict = mg.extend_combined_dictionary(res_df, meta_dict, dynamic_dict)
                    if not combined_res_df.empty:
                        io.write_results_to_gff(combined_res_df, os.path.join(res_path, "combined_results"), wildcard, split_gff)
                        io.write_results_to_table(combined_res_df, os.path.join(res_path, "combined_results"), wildcard)
                        combined_meta_dict, combined_dynamic_dict = mg.extend_combined_dictionary(combined_res_df, combined_meta_dict, combined_dynamic_dict)

                if meta_dict:
                    mg.write_merged_table(meta_dict, dynamic_dict, os.path.join(res_path, "%s_final.xlsx" % experiment), contrasts)
                    mg.write_merged_gff(meta_dict, os.path.join(res_path, "%s_final.xlsx" % experiment))
                if combined_meta_dict:
                    mg.write_merged_table(combined_meta_dict, combined_dynamic_dict, os.path.join(res_path, "%s_combined_final.xlsx" % experiment), contrasts)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    parser.add_argument("--peak_height_operator", action="store", dest="peak_height_operator", default="max"
                                                , help="{max,sum}:\n"\
                                                      +"'max': within the codon interval select the highest value (> min_peak_height)\n"\
                                                      +"'sum': within the codon interval sum up all values (> min_peak_height)")
    parser.add_argument("--tts_start_selection", action="store", dest="tts_start_selection", default="furthest_inframe"\
                                               , help="{furthest_inframe, next_inframe}\n"\
                                                      "'furthest_inframe': select the furthest inframe start codon that, without overstepping the next inframe stop codon.\n"\
                                                      "'next_inframe': select the closest inframe start codon.")
    parser.add_argument("--min_peak_height", action="store", dest="min_peak_height", default=5, type=int\
                                           , help="Minimum height value to be considered a peak. (max option)\n"\
                                                 +"Minimum height value to be added to the total peak value (sum option)")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=100, help="The max length to take into account when using the combination method for TIS+TTS.")
    parser.add_argument("--all_reads_rpkm", action="store_true", dest="all_reads_rpkm", help="If set, all mapped reads will be used for the calculation of RPKM values.\n"\
                                                                                            +"By default only mapped reads of the specified lengths will be used")

    parser.add_argument("--log_fold_contrasts", nargs="+", default=[], help="List of contrasts for which log2fold change will be calculated. (e.g. TIS-A-1_RIBO-A-1")
    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")
    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_ORFBounder(config_df, args.tts_start_selection, args.min_peak_height, args.peak_height_operator, args.max_ORF_length, args.split_gff, args.result_path, args.log_fold_contrasts, args.all_reads_rpkm)


if __name__ == '__main__':
    main()
