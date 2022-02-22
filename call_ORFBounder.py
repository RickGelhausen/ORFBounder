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
    expected_columns = ["experiment_name",\
                        "annotation_file_path", "genome_file_path", "alignment_folder_path",\
                        "RIBO_file_path", "TIS_file_path", "TTS_file_path",\
                        "normalization_method", "mapping_method", "offset_file_path", "read_lengths",\
                        "min_peak_height", "peak_height_operator", "tts_start_selection",\
                        "log_fold_contrasts", "max_ORF_length", "rpkm_read_usage",\
                        "gff_output_mode", "start_codons", "stop_codons"]

    input_columns = config_df.columns

    if list(set(expected_columns) - set(input_columns)) != []:
        msg.error_list("Config Sheet columns are incomplete:\n",\
                       "Ensure that the file is TAB seperated.",\
                       "Required columns:",\
                       expected_columns,\
                       input_columns)

    for row in config_df.itertuples(index=False, name="Pandas"):
        # Required
        experiment = getattr(row, "experiment_name")
        annotation = getattr(row, "annotation_file_path")
        genome = getattr(row, "genome_file_path")
        file_path_tis = getattr(row, "TIS_file_path")
        file_path_tts = getattr(row, "TTS_file_path")
        file_path_ribo = getattr(row, "RIBO_file_path")
        mapping_method = getattr(row, "mapping_method")
        normalization = getattr(row, "normalization_method")
        offset_json = getattr(row, "offset_file_path")

        # Optional
        read_lengths = getattr(row, "read_lengths")
        bam_folder = getattr(row, "alignment_folder_path")
        min_peak_height = getattr(row, "min_peak_height")
        peak_height_operator = getattr(row, "peak_height_operator")
        tts_start_selection = getattr(row, "tts_start_selection")
        log_fold_contrasts = getattr(row, "log_fold_contrasts")
        max_ORF_length = getattr(row, "max_ORF_length")
        rpkm_read_usage = getattr(row, "rpkm_read_usage")
        gff_output_mode = getattr(row, "gff_output_mode")


        # TIS / TTS / RIBO check
        with_tis = True
        with_tts = True
        with_ribo = True
        if is_empty(file_path_tis):
            with_tis = False
        if is_empty(file_path_tts):
            with_tts = False
        if is_empty(file_path_ribo):
            with_ribo = False
        if not with_tis and not with_tts:
            msg.error("No TIS or TTS path given! Specify atleast one.")
        if with_tis and not Path(file_path_tis).is_dir():
            msg.error(f"Mapping TIS directory is not valid! Ensure to enter a correct path!\n{file_path_tis}")
        if with_tts and not Path(file_path_tts).is_dir():
            msg.error(f"Mapping TTS directory is not valid! Ensure to enter a correct path!\n{file_path_tts}")
        if with_ribo and not Path(file_path_ribo).is_dir():
            msg.error(f"Mapping RIBO directory is not valid! Ensure to enter a correct path!\n{file_path_ribo}")

        # Required parameter check
        if is_empty(experiment):
            msg.error("Empty entry found: Missing experiment_name!")

        if is_empty(annotation):
            msg.error("Empty entry found: Missing annotation_file_path!")
        if not Path(annotation).is_file():
            msg.error(f"Annotation file is not valid! Ensure to enter a correct file path!\n{annotation}")

        if is_empty(genome):
            msg.error("Empty entry found: Missing genome_file_path!")
        if not Path(genome).is_file():
            msg.error(f"Genome file is not valid! Ensure to enter a correct file path!\n{genome}")

        if is_empty(normalization):
            msg.error("Empty entry found: Missing normalization_method!")
        for norm in normalization.split(","):
            if norm not in ["raw", "mil", "min"]:
                msg.error(f"Given normalization method is not allowed: {norm}.\n Choose from [raw, mil, min].")

        if is_empty(mapping_method):
            msg.error("Empty entry found: Missing mapping_method!")
        for mapping in mapping_method.split(","):
            if mapping not in ["threeprime", "fiveprime", "centered", "global"]:
                msg.error(f"Given mapping method is not allowed: {mapping}.\n Choose from [fiveprime, threeprime, centered, global].")

        if is_empty(offset_json):
            msg.error("Empty entry found: Missing offset_file_path!")
        if not Path(offset_json).is_file():
            msg.error(f"Offsets file is not valid! Ensure to enter a correct file path!\n{offset_json}")


        # Optional parameters
        if is_empty(read_lengths):
            msg.warning("No read lengths specified, using default: -1 (all read lengths).")

        if bam_folder != "" and isinstance(bam_folder, str):
            if not Path(bam_folder).is_dir():
                msg.error(f"Error: Given alignment_folder_path does not exist, either provide no bamfolder or an existing one!\n{bam_folder}")

        if is_empty(min_peak_height):
            msg.warning("No minimum peak length specfied, using default: 5.")
        else:
            if not str(min_peak_height).isnumeric() or "." in str(min_peak_height):
                msg.error("Error: Non-numerical or float value given for min_peak_length!")
            else:
                if int(min_peak_height) < 0:
                    msg.error("Error: Negative min_peak_heigth given.")

        if is_empty(peak_height_operator):
            msg.warning("No peak_height_operator specified, using default: max.")
        else:
            if peak_height_operator not in ["sum", "max"]:
                msg.error(f"Error: Given peak_height_operator does not exist: {peak_height_operator}. Use [sum, max]")

        if is_empty(tts_start_selection):
            msg.warning("No tts_start_selection specified, using default: furthest_inframe.")
        else:
            if tts_start_selection not in ["furthest_inframe", "next_inframe"]:
                msg.error(f"Error: Given tts_start_selection does not exist {tts_start_selection}. Use [furthest_inframe, next_inframe]")

        if is_empty(rpkm_read_usage):
            msg.warning("No rpkm_read_usage specified, using default: all.")
        else:
            if rpkm_read_usage not in ["all", "specific"]:
                msg.error(f"Error: Given rpkm_read_usage does not exist {rpkm_read_usage}. Use [all, specific]")

        if is_empty(gff_output_mode):
            msg.warning("No gff_output_mode specified, using default: combined.")
        else:
            if gff_output_mode not in ["combined", "split"]:
                msg.error(f"Error: Given gff_output_mode does not exist {gff_output_mode}. Use [combined, split]")

        if is_empty(log_fold_contrasts):
            msg.warning("No log_fold_contrasts given, skipping!")

        if is_empty(max_ORF_length):
            msg.warning("No max_ORF_length specfied, using default: 150.")
        else:
            if not str(max_ORF_length).isnumeric() or "." in str(max_ORF_length):
                msg.error("Error: Non-numerical or float value given for max_ORF_length!")
            else:
                if int(max_ORF_length) < 0:
                    msg.error("Error: Negative max_ORF_length given.")

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

def call_ORFBounder(config_df, result_path):
    """
    Run the ORFBounder experiments specified in the config sheet.
    """

    for row in config_df.itertuples(index=False, name="Pandas"):
        # Required
        experiment = getattr(row, "experiment_name")
        annotation = getattr(row, "annotation_file_path")
        genome = getattr(row, "genome_file_path")
        file_path_tis = getattr(row, "TIS_file_path")
        file_path_tts = getattr(row, "TTS_file_path")
        file_path_ribo = getattr(row, "RIBO_file_path")
        read_lengths = getattr(row, "read_lengths")
        mapping_method = getattr(row, "mapping_method").split(",")
        normalization = getattr(row, "normalization_method").split(",")
        offset_json = getattr(row, "offset_file_path")
        # Optional
        read_lengths = getattr(row, "read_lengths")
        bam_folder = getattr(row, "alignment_folder_path")
        min_peak_height = getattr(row, "min_peak_height")
        peak_height_operator = getattr(row, "peak_height_operator")
        tts_start_selection = getattr(row, "tts_start_selection")
        log_fold_contrasts = getattr(row, "log_fold_contrasts")
        max_orf_length = getattr(row, "max_ORF_length")
        rpkm_read_usage = getattr(row, "rpkm_read_usage")
        gff_output_mode = getattr(row, "gff_output_mode")
        start_codons = getattr(row, "start_codons")
        stop_codons = getattr(row, "stop_codons")

        if is_empty(max_orf_length):
            max_orf_length = 150
        else:
            max_orf_length = int(max_orf_length)

        if is_empty(peak_height_operator):
            peak_height_operator = "max"

        if is_empty(tts_start_selection):
            tts_start_selection = "furthest_inframe"

        if is_empty(rpkm_read_usage):
            all_reads_rpkm = True
        else:
            if rpkm_read_usage == "specific":
                all_reads_rpkm = False
            else:
                all_reads_rpkm = True

        if is_empty(gff_output_mode):
            split_gff = False
        else:
            if gff_output_mode == "split":
                split_gff = True
            else:
                split_gff = False

        if is_empty(start_codons):
            start_codons = ["ATG", "GTG", "TTG"]
        else:
            start_codons = [codon.strip(" ") for codon in start_codons.split(",")]

        if is_empty(stop_codons):
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
                    mg.write_merged_table(meta_dict, dynamic_dict, os.path.join(res_path, "%s_final.xlsx" % experiment), log_fold_contrasts)
                    mg.write_merged_gff(meta_dict, os.path.join(res_path, "%s_final.xlsx" % experiment))
                if combined_meta_dict:
                    mg.write_merged_table(combined_meta_dict, combined_dynamic_dict, os.path.join(res_path, "%s_combined_final.xlsx" % experiment), log_fold_contrasts)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")

    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_ORFBounder(config_df, args.tts_start_selection, args.min_peak_height, args.peak_height_operator, args.max_ORF_length, args.split_gff, args.result_path, args.log_fold_contrasts, args.all_reads_rpkm)


if __name__ == '__main__':
    main()
