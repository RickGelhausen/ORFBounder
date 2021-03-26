#!/usr/bin/env python
import os,sys
import re
import json
import argparse
import math

import pandas as pd
from pathlib import Path

import ORFBounder as ob
import lib.merging as mg

def check_config_sheet(config_sheet):
    """
    Check whether the config sheet is correctly formatted
    """

    config_df = pd.read_csv(config_sheet, sep="\t")

    if sorted(config_df.columns) != ["Annotation", "Bamfolder", "Experiment", "Genome", "Mapping", "Normalization", "Offsets", "Readthreshold", "Startcodons", "Stopcodons"]:
        sys.exit("Config Sheet columns are incomplete:\n\
                Required columns: Experiment,Annotation,Genome,Mapping,Normalization,Offsets,Bamfolder,Startcodons,Stopcodons,Readthreshold\n\
                Ensure that the file is TAB seperated.")

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        mapping = getattr(row, "Mapping")
        normalization = getattr(row, "Normalization")
        bamfolder = getattr(row, "Bamfolder")
        offsets = getattr(row, "Offsets")

        if experiment == "":
            sys.exit("Empty entry found: Missing Experiment!")
        if annotation == "":
            sys.exit("Empty entry found: Missing Annotation!")
        if genome == "":
            sys.exit("Empty entry found: Missing Genome!")
        if mapping == "":
            sys.exit("Empty entry found: Missing Mapping!")
        if normalization == "":
            sys.exit("Empty entry found: Missing Normalization!")
        if offsets == "":
            sys.exit("Empty entry found: Missing Offsets!")

        if not Path(annotation).is_file():
            sys.exit("Annotation file is not valid! Ensure to enter a correct file path!\n%s" % annotation)
        if not Path(genome).is_file():
            sys.exit("Genome file is not valid! Ensure to enter a correct file path!\n%s" % genome)
        if not Path(mapping).is_dir():
            sys.exit("Mapping directory is not valid! Ensure to enter a correct path!\n%s" % mapping)
        if not Path(offsets).is_file():
            sys.exit("Offsets file is not valid! Ensure to enter a correct file path!\n%s" % offsets)
        if bamfolder != "" and not Path(bamfolder).is_dir():
            sys.exit("Given bamfolder is non-existant, either provide no bamfolder or an existing one!\n%s" % bamfolder)

        for norm in normalization.split(","):
            norm_path = os.path.join(mapping, norm)
            if not Path(norm_path).is_dir():
                sys.exit("Normalization path does not exist: %s" % norm_path)
    return config_df

def retrieve_wig_information(wig_path, offset_dict):
    """
    Create a list of matching TIS/TTS condition+replicate files to run together.
    """
    _, _, wig_files = next(os.walk(wig_path))

    tt_files = [wig for wig in wig_files if (("TIS" in wig or "TTS" in wig) and not "RNA" in wig) and (wig.endswith(".wig")) ]

    sample_dict = {}
    for file in tt_files:
        wildcard = re.split('_|\.', os.path.basename(file))[0]
        method, condition, replicate = wildcard.split("-")

        if (condition, replicate) not in sample_dict:
            if (method == "TIS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)] = [os.path.join(wig_path, file),"","",""]
            elif (method == "TIS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)] = ["",os.path.join(wig_path, file),"",""]
            elif (method == "TTS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)] = ["","",os.path.join(wig_path, file),""]
            elif (method == "TTS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)] = ["","","",os.path.join(wig_path, file)]

        else:
            if (method == "TIS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)][0] = os.path.join(wig_path, file)
            elif (method == "TIS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)][1] = os.path.join(wig_path, file)
            elif (method == "TTS") and ("fwd" in file or "forward" in file):
                sample_dict[(condition, replicate)][2] = os.path.join(wig_path, file)
            elif (method == "TTS") and ("rev" in file or "reverse" in file):
                sample_dict[(condition, replicate)][3] = os.path.join(wig_path, file)

    wig_list = []
    for key, val in sample_dict.items():
        if (val[0] != "" and val[1] != "") or (val[2] != "" and val[3] != ""):
            tis_offset, tts_offset = "", ""
            if "TIS" in offset_dict:
                if "TIS-%s-%s" % key in offset_dict["TIS"]:
                    tis_offset = offset_dict["TIS"]["TIS-%s-%s" % key]
                elif "default" in offset_dict["TIS"]:
                    tis_offset = offset_dict["TIS"]["default"]
                else:
                    sys.exit("Wrongly formatted JSON file, missing default value!")

            if "TTS" in offset_dict:
                if "TTS-%s-%s" % key in offset_dict["TTS"]:
                    tts_offset = offset_dict["TTS"]["TTS-%s-%s" % key]
                elif "default" in offset_dict["TTS"]:
                    tts_offset = offset_dict["TTS"]["default"]
                else:
                    sys.exit("Wrongly formatted JSON file, missing default value!")

            wig_list.append((val, "%s-%s" % key, tis_offset, tts_offset))

    return wig_list


def call_ORFBounder(config_df, use_longest_TTS_ORF, max_ORF_length, split_gff, result_path):
    """
    Run the ORFBounder experiments specified in the config sheet.
    """

    for row in config_df.itertuples(index=False, name="Pandas"):
        experiment = getattr(row, "Experiment")
        annotation = getattr(row, "Annotation")
        genome = getattr(row, "Genome")
        mapping = getattr(row, "Mapping")
        normalization = getattr(row, "Normalization").split(",")
        offsets = getattr(row, "Offsets")
        start_codons = getattr(row, "Startcodons")
        stop_codons = getattr(row, "Stopcodons")
        read_threshold = getattr(row, "Readthreshold")
        bamfolder = getattr(row, "Bamfolder")

        if read_threshold == "" or math.isnan(read_threshold):
            read_threshold = 5
        else:
            read_threshold = int(read_threshold)

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

        with open(offsets, "r") as f:
            offset_data = json.load(f)

        for norm in normalization:
            meta_dict, dynamic_dict = {}, {}

            wig_path = os.path.join(mapping, norm)

            wig_list = retrieve_wig_information(wig_path, offset_data)

            res_path = os.path.join(result_path, experiment, norm)
            for (TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig), conrep, tis_offset, tts_offset in wig_list:
                result_df = ob.run_ORFBounder(TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig, bamfolder, annotation, \
                                            genome, start_codons, stop_codons, res_path, conrep, tis_offset, tts_offset, \
                                            use_longest_TTS_ORF, read_threshold, split_gff, max_ORF_length)

                meta_dict, dynamic_dict = mg.extend_combined_dictionary(result_df, meta_dict, dynamic_dict)
            mg.write_merged_table(meta_dict, dynamic_dict, os.path.join(res_path, "%s_final.xlsx" % experiment))

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.")

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("--split_gff", action="store_true", dest="split_gff", help="Split gff into one for each gene_type.")
    parser.add_argument("--use_longest_TTS_ORF", action="store_true", dest="use_longest_TTS_ORF", help="Use the furthest possible inframe start codon for each detected stop codon to form the longest possible ORF that contains only one inframe stop codon. \
                                                                                                        Default uses the first detected start codon and may result in very short ORFs.")
    parser.add_argument("--max_ORF_length", action="store", dest="max_ORF_length", type=int, default=100, help="The max length to take into account when using the combination method for TIS+TTS.")

    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")
    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_ORFBounder(config_df, args.use_longest_TTS_ORF, args.max_ORF_length, args.split_gff, args.result_path)


if __name__ == '__main__':
    main()
