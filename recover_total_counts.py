from numpy import diff
import pysam
import argparse
import math
import pandas as pd
import re
import os

from pathlib import Path
import lib.messaging as msg
import lib.io as io

def is_empty(entry):
    """
    Check if the current table entry is empty
    """

    return entry == "" or (isinstance(entry, float) and math.isnan(entry))

def check_config_sheet(config_sheet):
    """
    Check whether the config sheet is correctly formatted
    """

    config_df = pd.read_csv(config_sheet, sep="\t")
    expected_columns = ["experiment_name", "RIBO_folder_path", "TIS_folder_path", "TTS_folder_path"]

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

        file_path_tis = getattr(row, "TIS_folder_path")
        file_path_tts = getattr(row, "TTS_folder_path")
        file_path_ribo = getattr(row, "RIBO_folder_path")
        read_length_json = getattr(row, "read_length_json")

        msg.message(f"Checking config file for: {experiment}")
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

        if is_empty(read_length_json):
            msg.warning("No read lengths specified, using default: -1 (all read lengths).")

    return config_df

def count_mapped_reads(alignment_file_path, read_count_dict, read_length_dict):
    """
    Read alignment file and count the number of mapped reads.
    """
    sample = alignment_file_path.stem
    msg.message(f">Counting reads for {sample}.")

    wildcard = re.split('_|\.', os.path.basename(alignment_file_path))[0]
    if read_length_dict == -1:
        read_length_list = -1
    elif wildcard in read_length_dict:
        read_length_list = read_length_dict[wildcard]
    elif "default" in read_length_dict:
        read_length_list = read_length_dict["default"]
    else:
        msg.warning("Warning no default value given for read-lengths. Using all read lengths.")
        read_length_list = -1

    alignment_file = pysam.AlignmentFile(alignment_file_path)
    try:
        for read in alignment_file.fetch():
            chrom = read.reference_name

            if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                continue

            start, stop = read.reference_start, read.reference_end - 1
            read_length = stop - start + 1

            if read_length_list != -1:
                if str(read_length) not in read_length_list:
                    continue

            if (sample, chrom) in read_count_dict:
                read_count_dict[(sample, chrom)] += 1
            else:
                read_count_dict[(sample, chrom)] = 1

    except ValueError:
        msg.error("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

    return read_count_dict

def write_read_counts_to_file(read_count_dict, experiment_name, result_path):
    """
    Create a tsv table with total mapped read counts.
    """

    with open(Path(result_path).joinpath(f"{experiment_name}_mapped_reads.tsv"), "w", encoding="utf-8") as f:
        for sample, chrom in read_count_dict.keys():
            f.write(f"{sample}\t{chrom}\t{read_count_dict[(sample, chrom)]}\n")


def recover_read_information(config_df, result_path):
    """
    For every desired experiment create a total mapped read file.
    """

    for row in config_df.itertuples(index=False):
        experiment_name = getattr(row, "experiment_name")
        msg.message(f"Creating read count file for experiment: {experiment_name} ...")

        file_path_tis = getattr(row, "TIS_folder_path")
        file_path_tts = getattr(row, "TTS_folder_path")
        file_path_ribo = getattr(row, "RIBO_folder_path")
        read_length_json = getattr(row, "read_length_json")

        differing_paths = set([file_path_tis, file_path_tts, file_path_ribo])
        files = []
        for path in differing_paths:
            if is_empty(path):
                continue
            files.extend([entry for entry in Path(path).glob("*.bam") if entry.is_file()])
            files.extend([entry for entry in Path(path).glob("*.sam") if entry.is_file()])

        read_length_dict = io.parse_read_lengths(read_length_json)
        read_count_dict = {}
        for alignment_file in sorted(files):
            read_count_dict = count_mapped_reads(alignment_file, read_count_dict, read_length_dict)

        write_read_counts_to_file(read_count_dict, experiment_name, result_path)

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-c","--config_sheet", action="store", dest="config_sheet", required=True, help="Config sheet containing information on experiments to be run.")
    parser.add_argument("-r","--result_path", action="store", dest="result_path", required=True, help="Path of the result folder.")

    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    recover_read_information(config_df, args.result_path)


if __name__ == '__main__':
    main()
