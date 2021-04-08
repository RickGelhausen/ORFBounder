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

        if experiment == "":
            msg.error("Empty entry found: Missing Experiment!")
        if annotation == "":
            msg.error("Empty entry found: Missing Annotation!")
        if genome == "":
            msg.error("Empty entry found: Missing Genome!")
        if mapping_tis == "":
            msg.error("Empty entry found: Missing Mapping_TIS!")
        if mapping_tts == "":
            msg.error("Empty entry found: Missing Mapping_TTS!")
        if normalization == "":
            msg.error("Empty entry found: Missing Normalization!")
        if offsets == "":
            msg.error("Empty entry found: Missing Offsets!")

        if not Path(annotation).is_file():
            msg.error("Annotation file is not valid! Ensure to enter a correct file path!\n%s" % annotation)
        if not Path(genome).is_file():
            msg.error("Genome file is not valid! Ensure to enter a correct file path!\n%s" % genome)
        if not Path(mapping_tis).is_dir():
            msg.error("Mapping TIS directory is not valid! Ensure to enter a correct path!\n%s" % mapping_tis)
        if not Path(mapping_tts).is_dir():
            msg.error("Mapping TTS directory is not valid! Ensure to enter a correct path!\n%s" % mapping_tts)
        if not Path(offsets).is_file():
            msg.error("Offsets file is not valid! Ensure to enter a correct file path!\n%s" % offsets)
        if bamfolder != "" and isinstance(bamfolder, str):
            if not Path(bamfolder).is_dir():
                msg.error("Given bamfolder is non-existant, either provide no bamfolder or an existing one!\n%s" % bamfolder)

        for norm in normalization.split(","):
            norm_path = os.path.join(mapping_tis, norm)
            if not Path(norm_path).is_dir():
                msg.error("Normalization path does not exist: %s" % norm_path)

            norm_path = os.path.join(mapping_tts, norm)
            if not Path(norm_path).is_dir():
                msg.error("Normalization path does not exist: %s" % norm_path)

    return config_df

def retrieve_wig_information(wig_path_tis, wig_path_tts, offset_dict):
    """
    Create a list of matching TIS/TTS condition+replicate files to run together.
    """
    _, _, wig_files_tis = next(os.walk(wig_path_tis))
    _, _, wig_files_tts = next(os.walk(wig_path_tts))

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

def dictionary_depth(dic):
    """
    get the depth of a dictionary
    """
    queue = deque([(id(dic), dic, 0)])
    already_visited = set()
    while queue:
        id_, o, level = queue.popleft()
        if id_ in already_visited:
            continue
        already_visited.add(id_)
        if isinstance(o, dict):
            queue += ((id(v), v, level + 1) for v in o.values())
    return level

def base_mapping(mapping):
    """
    retrieve basename of mapping (only useful for HRIBO output)
    """

    if "fiveprime" in mapping:
        return "fiveprime"
    elif "threeprime" in mapping:
        return "threeprime"
    elif "global" in mapping:
        return "global"
    elif "centered" in mapping:
        return "centered"

def build_offset_dictionary(offset_file, genome, mapping_tis, mapping_tts):
    """
    read a json offset file and process
    """
    with open(offset_file, "r") as f:
        offset_data = json.load(f)
    print(offset_data)
    print(dictionary_depth(offset_data))
    if dictionary_depth(offset_data) == 2:
        return offset_data

    else:
        new_dict = {}
        chromosome_name = ""
        cur_length = 0
        for key, val in io.generate_genome_dict(genome).items():
            if len(val) > cur_length:
                chromosome_name = key
                cur_length = len(val)

        mapping_tis = base_mapping(os.path.basename(mapping_tis))
        mapping_tts = base_mapping(os.path.basename(mapping_tts))

        for method, norm_dict in offset_data.items():
            if method.lower() not in ["tis", "tts"]:
                continue

            total_offset = 0
            total_count = 0
            for norm, sample_dict in norm_dict.items():
                if norm.lower() != "raw":
                    continue

                for mapping, readlength_dict in sample_dict.items():
                    cur_chrom, cur_mapping = mapping.split("_")
                    if cur_mapping in [mapping_tis, mapping_tts] and chromsome == cur_chrom:
                        if "raw" in readlength_dict:
                            offset = int(read_dict["raw"].split(",")[0])

                        elif "mean" in readlength_dict:
                            offset = int(read_dict["mean"].split(",")[0])

                        elif len(readlength_dict) != 0:
                            counter = 0
                            offset = 0
                            for readlength, value in readlength_dict.items():

                                if readlength.lower() in ["raw", "mean"]:
                                    continue
                                counter += 1
                                offset += int(value.split(","))
                            offset = int(offset / counter)
                        else:
                            msg.error("Error: empty readlength data in JSON file")

                    else:
                        continue

                    if method in new_dict:
                        new_dict[method][sample] = offset
                    else:
                        new_dict[method] = {sample : offset}

                    total_count += 1
                    total_offset = offset
            new_dict[method]["default"] = int(total_offset / total_count)

    return new_dict

def call_ORFBounder(config_df, TTS_start_selection, min_peak_height, peak_height_calculation, max_ORF_length, split_gff, result_path):
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

        offset_data = build_offset_dictionary(offset_file, genome, mapping_tis, mapping_tts)

        for norm in normalization:
            meta_dict, dynamic_dict = {}, {}
            combined_meta_dict, combined_dynamic_dict = {}, {}

            wig_path_tis = os.path.join(mapping_tis, norm)
            wig_path_tts = os.path.join(mapping_tts, norm)

            wig_list = retrieve_wig_information(wig_path_tis, wig_path_tts, offset_data)

            res_path = os.path.join(result_path, experiment, norm)
            for (TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig), conrep, tis_offset, tts_offset in wig_list:
                try:
                    res_df, combined_res_df = ob.run_ORFBounder(TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig, bamfolder, \
                                                            annotation, genome, start_codons, stop_codons, res_path, conrep, \
                                                            tis_offset, tts_offset, TTS_start_selection, min_peak_height, \
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
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.")

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    parser.add_argument("--peak_height_calculation", action="store", dest="peak_height_calculation", default="max"
                                                   , help="{max,sum}:\n"\
                                                         +"'max': within the codon interval select the highest value (> min_peak_height)"\
                                                         +"'sum': within the codon interval sum all values (> min_peak_height)")
    parser.add_argument("--TTS_start_selection", action="store", dest="TTS_start_selection", default="furthest_inframe"\
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
    call_ORFBounder(config_df, args.TTS_start_selection, args.min_peak_height, args.peak_height_calculation, args.max_ORF_length, args.split_gff, args.result_path)


if __name__ == '__main__':
    main()
