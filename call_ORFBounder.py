#!/usr/bin/env python
import os,sys
import re
import json
import argparse

import pandas as pd
from pathlib import Path

import ORFBounder as ob

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
            sys.exit("Annotation file is not valid! Ensure to enter a correct file path!")
        if not Path(genome).is_file():
            sys.exit("Genome file is not valid! Ensure to enter a correct file path!")
        if not Path(mapping).is_dir():
            sys.exit("Mapping directory is not valid! Ensure to enter a correct path!")
        if not Path(offsets).is_file():
            sys.exit("Offsets file is not valid! Ensure to enter a correct file path!")
        if bamfolder != "" and not Path(bamfolder).is_dir():
            sys.exit("Given bamfolder is not existing, either provide no bamfolder or an existing one!")

        for norm in normalization.split(",")
            norm_path = os.path.join(mapping, norm)
            if not Path(norm_path).is_dir():
                sys.exit("Normalization path does not exist: %s" % norm_path)
    return config_df

def build_wig_tuples(tt_files, wig_path):
    """
    Create a list of matching TIS/TTS condition+replicate files to run together.
    """

    sample_dict = {}
    for file in tt_files:
        wildcard = re.split('_|\.', os.path.basename(file))[0]
        method, condition, replicate = wildcard.split("-")

        if (condition, replicate) not in sample_dict:
            if (method == "TIS") and ("fwd" or "forward" in file):
                sample_dict[(condition, replicate)] = [os.path(wig_path, file),"","",""]
            elif (method == "TIS") and ("rev" or "reverse" in file):
                sample_dict[(condition, replicate)] = ["",os.path(wig_path, file),"",""]
            elif (method == "TTS") and ("fwd" or "forward" in file):
                sample_dict[(condition, replicate)] = ["","",os.path(wig_path, file),""]
            elif (method == "TTS") and ("rev" or "reverse" in file):
                sample_dict[(condition, replicate)] = ["","","",os.path(wig_path, file)]

        else:
            if (method == "TIS") and ("fwd" or "forward" in file):
                sample_dict[(condition, replicate)][0] = os.path(wig_path, file)
            elif (method == "TIS") and ("rev" or "reverse" in file):
                sample_dict[(condition, replicate)][1] = os.path(wig_path, file)
            elif (method == "TTS") and ("fwd" or "forward" in file):
                sample_dict[(condition, replicate)][2] = os.path(wig_path, file)
            elif (method == "TTS") and ("rev" or "reverse" in file):
                sample_dict[(condition, replicate)][3] = os.path(wig_path, file)

    wig_list = []
    for key, val in sample_dict.items():
        if (val[0] != "" and val[1] != "") or (val[2] != "" and val[3] != "")
            wig_list.append(val)

    return wig_list


def call_ORFBounder(config_df, result_path):
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

        read_threshold = 5 if read_threshold == "" else read_threshold = int(read_threshold)



        if start_codon =
        with open(offsets, "r") as f:
            offset_data = json.load(f)

        for norm in normalization:
            meta_dict, dynamic_dict = {}, {}

            wig_path = os.path.join(mapping, norm)
            _, _, wig_files = next(os.walk(wig_path))

            tt_files = [wig for wig in wig_files if ("TIS" or "TTS") and not "RNA" in wig]

            wig_list = build_wig_tuples(tt_files, wig_path)

            for TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig in wig_list:
                result_df = ob.run_ORFBounder(TIS_fwd_wig, TIS_rev_wig, TTS_fwd_wig, TTS_rev_wig, annotation, genome, \
                )


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.")

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")
    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_ORFBounder(config_df, args.result_path)



if __name__ == '__main__':
    main()
