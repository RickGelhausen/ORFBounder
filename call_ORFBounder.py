#!/usr/bin/env python
import os
import re
import json
import argparse
import math

import pandas as pd
from pathlib import Path
from collections import deque

import ORFBounder as ob
import lib.merging as mg
import lib.io as io
import lib.messaging as msg
import lib.misc as misc
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

    if sorted(config_df.columns) != ["Annotation", "Bam_folder", "Experiment", "Genome", "Mapping_TIS", "Mapping_TTS", "Normalization", "Offsets", "Start_codons", "Stop_codons"]:
        msg.error("Config Sheet columns are incomplete:\n\
                Required columns: Experiment,Annotation,Genome,Mapping_TIS,Mapping_TTS,Normalization,Offsets,Bam_folder,Start_codons,Stop_codons\n\
                Ensure that the file is TAB seperated.")

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        mapping_tis = getattr(row, "Mapping_TIS")
        mapping_tts = getattr(row, "Mapping_TTS")
        normalization = getattr(row, "Normalization")
        bamfolder = getattr(row, "Bam_folder")
        offsets = getattr(row, "Offsets")

        with_tis = True
        with_tts = True
        if is_empty(experiment):
            msg.error("Empty entry found: Missing Experiment!")
        if is_empty(annotation):
            msg.error("Empty entry found: Missing Annotation!")
        if is_empty(genome):
            msg.error("Empty entry found: Missing Genome!")
        if is_empty(mapping_tis):
            with_tis = False
        if is_empty(mapping_tts):
            with_tts = False
        if is_empty(normalization):
            msg.error("Empty entry found: Missing Normalization!")
        if is_empty(offsets):
            msg.error("Empty entry found: Missing Offsets!")

        if not with_tis and not with_tts:
            msg.error("No TIS or TTS path given! Specify atleast one.")

        if not Path(annotation).is_file():
            msg.error("Annotation file is not valid! Ensure to enter a correct file path!\n%s" % annotation)
        if not Path(genome).is_file():
            msg.error("Genome file is not valid! Ensure to enter a correct file path!\n%s" % genome)
        if with_tis and not Path(mapping_tis).is_dir():
            msg.error("Mapping TIS directory is not valid! Ensure to enter a correct path!\n%s" % mapping_tis)
        if with_tts and not Path(mapping_tts).is_dir():
            msg.error("Mapping TTS directory is not valid! Ensure to enter a correct path!\n%s" % mapping_tts)
        if not Path(offsets).is_file():
            msg.error("Offsets file is not valid! Ensure to enter a correct file path!\n%s" % offsets)
        if bamfolder != "" and isinstance(bamfolder, str):
            if not Path(bamfolder).is_dir():
                msg.error("Given bamfolder is non-existant, either provide no bamfolder or an existing one!\n%s" % bamfolder)

        for norm in normalization.split(","):
            if with_tis:
                norm_path = os.path.join(mapping_tis, norm)
                if not Path(norm_path).is_dir():
                    msg.error("Normalization path does not exist: %s" % norm_path)

            if with_tts:
                norm_path = os.path.join(mapping_tts, norm)
                if not Path(norm_path).is_dir():
                    msg.error("Normalization path does not exist: %s" % norm_path)

    return config_df

def retrieve_wig_information(wig_path_tis, wig_path_tts, offset_dict):
    """
    Create a list of matching TIS/TTS condition+replicate files to run together.
    """
    if wig_path_tis != "":
        _, _, wig_files_tis = next(os.walk(wig_path_tis))
    else:
        wig_files_tis = []

    if wig_path_tts != "":
        _, _, wig_files_tts = next(os.walk(wig_path_tts))
    else:
        wig_files_tts = []

    tt_files = []
    tt_files.extend([wig for wig in wig_files_tis if ("TIS" in wig and not "RNA" in wig) and (wig.endswith(".wig"))])
    tt_files.extend([wig for wig in wig_files_tts if ("TTS" in wig and not "RNA" in wig) and (wig.endswith(".wig"))])

    sample_dict = {}
    for file in tt_files:
        wildcard = re.split('_|\.', os.path.basename(file))[0]
        method, condition, replicate = wildcard.split("-")

        if (condition, replicate) not in sample_dict:
            if (method == "TIS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)] = [os.path.join(wig_path_tis, file),"","",""]
            elif (method == "TIS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)] = ["",os.path.join(wig_path_tis, file),"",""]
            elif (method == "TTS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)] = ["","",os.path.join(wig_path_tts, file),""]
            elif (method == "TTS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)] = ["","","",os.path.join(wig_path_tts, file)]

        else:
            if (method == "TIS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)][0] = os.path.join(wig_path_tis, file)
            elif (method == "TIS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)][1] = os.path.join(wig_path_tis, file)
            elif (method == "TTS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)][2] = os.path.join(wig_path_tts, file)
            elif (method == "TTS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)][3] = os.path.join(wig_path_tts, file)

    wig_list = []
    for key, val in sample_dict.items():
        if (val[0] != "" and val[1] != "") or (val[2] != "" and val[3] != ""):
            tis_offset, tts_offset = 15, 15

            if "TIS" in offset_dict:
                if "TIS-%s-%s" % key in offset_dict["TIS"]:
                    tis_offset = offset_dict["TIS"]["TIS-%s-%s" % key]
                elif "default" in offset_dict["TIS"]:
                    tis_offset = offset_dict["TIS"]["default"]
                else:
                    msg.warning("Wrongly formatted JSON file, missing default value! (using 15 instead)")

            if "TTS" in offset_dict:
                if "TTS-%s-%s" % key in offset_dict["TTS"]:
                    tts_offset = offset_dict["TTS"]["TTS-%s-%s" % key]
                elif "default" in offset_dict["TTS"]:
                    tts_offset = offset_dict["TTS"]["default"]
                else:
                    msg.warning("Wrongly formatted JSON file, missing default value! (using 15 instead)")

            wig_list.append((val, "%s-%s" % key, tis_offset, tts_offset))

    return wig_list

def call_ORFBounder(config_df, tts_start_selection, min_peak_height, peak_height_calculation, max_ORF_length, split_gff, result_path):
    """
    Run the ORFBounder experiments specified in the config sheet.
    """

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        mapping_tis = getattr(row, "Mapping_TIS")
        mapping_tts = getattr(row, "Mapping_TTS")
        normalization = getattr(row, "Normalization").split(",")
        offset_file = getattr(row, "Offsets")
        start_codons = getattr(row, "Start_codons")
        stop_codons = getattr(row, "Stop_codons")
        bamfolder = getattr(row, "Bam_folder")

        if max_ORF_length == "" or math.isnan(max_ORF_length):
            max_ORF_length = 100
        else:
            max_ORF_length = int(max_ORF_length)

        if start_codons == "" or math.isnan(start_codons):
            start_codons = ["ATG", "GTG", "TTG"]
        else:
            start_codons = [codon.strip(" ") for codon in start_codons.split(",")]

        if stop_codons == "" or math.isnan(stop_codons):
            stop_codons = ["TAG", "TAA", "TGA"]
        else:
            stop_codons = [codon.strip(" ") for codon in stop_codons.split(",")]

        offset_data = misc.build_offset_dictionary(offset_file, genome, mapping_tis, mapping_tts)

        for norm in normalization:
            meta_dict, dynamic_dict = {}, {}
            combined_meta_dict, combined_dynamic_dict = {}, {}

            try:
                wig_path_tis = os.path.join(mapping_tis, norm)
            except TypeError:
                wig_path_tis = ""

            try:
                wig_path_tts = os.path.join(mapping_tts, norm)
            except TypeError:
                wig_path_tts = ""

            wig_list = retrieve_wig_information(wig_path_tis, wig_path_tts, offset_data)

            res_path = os.path.join(result_path, experiment, norm)
            for (tis_fwd_wig, tis_rev_wig, tts_fwd_wig, tts_rev_wig), conrep, tis_offset, tts_offset in wig_list:
                try:
                    res_df, combined_res_df = ob.run_ORFBounder(tis_fwd_wig, tis_rev_wig, tts_fwd_wig, tts_rev_wig, bamfolder, \
                                                            annotation, genome, start_codons, stop_codons, res_path, conrep, \
                                                            tis_offset, tts_offset, tts_start_selection, min_peak_height, \
                                                            split_gff, max_ORF_length, peak_height_calculation)
                except SystemExit:
                    msg.warning("Error encountered while calling ORFBounder! Moving to next run!")
                    continue

                io.write_results_to_gff(res_df, os.path.join(res_path, "result_tables"), conrep, split_gff)
                io.write_results_to_table(res_df, os.path.join(res_path, "result_tables"), conrep)

                meta_dict, dynamic_dict = mg.extend_combined_dictionary(res_df, meta_dict, dynamic_dict)
                if not combined_res_df.empty:
                    io.write_results_to_gff(combined_res_df, os.path.join(res_path, "combined_results"), conrep, split_gff)
                    io.write_results_to_table(combined_res_df, os.path.join(res_path, "combined_results"), conrep)
                    combined_meta_dict, combined_dynamic_dict = mg.extend_combined_dictionary(combined_res_df, combined_meta_dict, combined_dynamic_dict)

            if meta_dict:
                mg.write_merged_table(meta_dict, dynamic_dict, os.path.join(res_path, "%s_final.xlsx" % experiment))
            if combined_meta_dict:
                mg.write_merged_table(combined_meta_dict, combined_dynamic_dict, os.path.join(res_path, "%s_combined_final.xlsx" % experiment))
def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    parser.add_argument("--peak_height_calculation", action="store", dest="peak_height_calculation", default="max"
                                                   , help="{max,sum}:\n"\
                                                         +"'max': within the codon interval select the highest value (> min_peak_height)\n"\
                                                         +"'sum': within the codon interval sum all values (> min_peak_height)")
    parser.add_argument("--tts_start_selection", action="store", dest="tts_start_selection", default="furthest_inframe"\
                                               , help="{furthest_inframe, next_inframe}\n"\
                                                      "'furthest_inframe': select the furthest inframe start codon that, without overstepping the next inframe stop codon.\n"\
                                                      "'next_inframe': select the closest inframe start codon.")
    parser.add_argument("--min_peak_height", action="store", dest="min_peak_height", default=5, type=int\
                                           , help="Minimum height value to be considered a peak. (max option)\n"\
                                                 +"Minimum height value to be added to the total peak value (sum option)")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=100, help="The max length to take into account when using the combination method for TIS+TTS.")
    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")
    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_ORFBounder(config_df, args.tts_start_selection, args.min_peak_height, args.peak_height_calculation, args.max_ORF_length, args.split_gff, args.result_path)


if __name__ == '__main__':
    main()
