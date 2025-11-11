#!/usr/bin/env python
"""
A wrapper script to call ORFBounder for multiple experiments based on a config sheet.
"""

from pathlib import Path

import os
import re
import argparse
import math

import pandas as pd

import orfbounder as ob
import lib.merging as mg
import lib.io as io
import lib.messaging as msg

# Constants
RESULT_TABLES_DIR = "result_tables"
COMBINED_RESULTS_DIR = "combined_results"
FINAL_OUTPUT_EXTENSION = ".xlsx"
FINAL_GFF_EXTENSION = ".gff"


def is_empty(entry: str | float) -> bool:
    """
    Check if the current table entry is empty.

    Args:
        entry: Table entry value to check

    Returns:
        True if entry is empty string or NaN, False otherwise
    """
    return entry == "" or (isinstance(entry, float) and math.isnan(entry))

def check_config_sheet(config_sheet: str) -> pd.DataFrame:
    """
    Check whether the config sheet is correctly formatted.

    Args:
        config_sheet: Path to the configuration file

    Returns:
        DataFrame containing the validated configuration

    Raises:
        SystemExit: If configuration is invalid
    """
    config_sheet_path = Path(config_sheet)
    config_df = pd.read_csv(config_sheet_path, sep="\t")

    expected_columns = [
        "experiment_name",
        "annotation_file_path", "genome_file_path", "alignment_folder_path",
        "RIBO_folder_path", "TIS_folder_path", "TTS_folder_path",
        "normalization_method", "mapped_counts_file_path", "mapping_method", "offset_file_path",
        "read_length_json", "min_peak_height", "peak_height_operator",
        "tts_start_selection", "rpkm_read_usage",
        "gff_output_mode", "start_codons", "stop_codons"
    ]

    input_columns = config_df.columns

    if list(set(expected_columns) - set(input_columns)):
        msg.error_list(
            "Config Sheet columns are incomplete:\n",
            "Ensure that the file is TAB seperated.",
            "Required columns:",
            expected_columns,
            input_columns
        )

    for row in config_df.itertuples(index=False, name="Pandas"):
        # Required
        experiment = getattr(row, "experiment_name")
        annotation = getattr(row, "annotation_file_path")
        genome = getattr(row, "genome_file_path")
        file_path_tis = getattr(row, "TIS_folder_path")
        file_path_tts = getattr(row, "TTS_folder_path")
        file_path_ribo = getattr(row, "RIBO_folder_path")
        mapping_method = getattr(row, "mapping_method")
        normalization = getattr(row, "normalization_method")
        offset_json = getattr(row, "offset_file_path")

        # Optional
        read_length_json = getattr(row, "read_length_json")
        bam_folder = getattr(row, "alignment_folder_path")
        min_peak_height = getattr(row, "min_peak_height")
        peak_height_operator = getattr(row, "peak_height_operator")
        tts_start_selection = getattr(row, "tts_start_selection")
        rpkm_read_usage = getattr(row, "rpkm_read_usage")
        gff_output_mode = getattr(row, "gff_output_mode")
        mapped_counts_file_path = getattr(row, "mapped_counts_file_path")

        msg.message(f"Checking config file for: {experiment}")

        # TIS / TTS / RIBO check
        with_tis = not is_empty(file_path_tis)
        with_tts = not is_empty(file_path_tts)
        with_ribo = not is_empty(file_path_ribo)

        if not with_tis and not with_tts:
            raise ValueError(msg.error("No TIS or TTS path given! Specify atleast one."))
        if with_tis and not Path(file_path_tis).is_dir():
            raise ValueError(msg.error(f"Mapping TIS directory is not valid! Ensure to enter a correct path!\n{file_path_tis}"))
        if with_tts and not Path(file_path_tts).is_dir():
            raise ValueError(msg.error(f"Mapping TTS directory is not valid! Ensure to enter a correct path!\n{file_path_tts}"))
        if with_ribo and not Path(file_path_ribo).is_dir():
            raise ValueError(msg.error(f"Mapping RIBO directory is not valid! Ensure to enter a correct path!\n{file_path_ribo}"))

        # Required parameter check
        if is_empty(experiment):
            raise ValueError(msg.error("Empty entry found: Missing experiment_name!"))

        if is_empty(annotation):
            raise ValueError(msg.error("Empty entry found: Missing annotation_file_path!"))
        if not Path(annotation).is_file():
            raise ValueError(msg.error(f"Annotation file is not valid! Ensure to enter a correct file path!\n{annotation}"))

        if is_empty(genome):
            raise ValueError(msg.error("Empty entry found: Missing genome_file_path!"))
        if not Path(genome).is_file():
            raise ValueError(msg.error(f"Genome file is not valid! Ensure to enter a correct file path!\n{genome}"))

        if is_empty(normalization):
            raise ValueError(msg.error("Empty entry found: Missing normalization_method!"))
        for norm in normalization.split(","):
            if norm not in ["raw", "mil", "min"]:
                raise ValueError(msg.error(f"Given normalization method is not allowed: {norm}.\n Choose from [raw, mil, min]."))

        if is_empty(mapping_method):
            raise ValueError(msg.error("Empty entry found: Missing mapping_method!"))
        for mapping in mapping_method.split(","):
            if mapping not in ["threeprime", "fiveprime", "centered", "global"]:
                raise ValueError(msg.error(f"Given mapping method is not allowed: {mapping}.\n Choose from [fiveprime, threeprime, centered, global]."))

        if is_empty(offset_json):
            raise ValueError(msg.error("Empty entry found: Missing offset_file_path!"))
        if not Path(offset_json).is_file():
            raise ValueError(msg.error(f"Offsets file is not valid! Ensure to enter a correct file path!\n{offset_json}"))

        # Optional parameters
        if is_empty(read_length_json):
            msg.warning("No read lengths specified, using default: None (all read lengths).")

        if bam_folder != "" and isinstance(bam_folder, str):
            if not Path(bam_folder).is_dir():
                raise ValueError(msg.error(f"Error: Given alignment_folder_path does not exist, either provide no bamfolder or an existing one!\n{bam_folder}"))

        if is_empty(min_peak_height):
            msg.warning("No minimum peak length specfied, using default: 5.")
        else:
            if not str(min_peak_height).isnumeric() or "." in str(min_peak_height):
                raise ValueError(msg.error("Error: Non-numerical or float value given for min_peak_length!"))
            else:
                if int(min_peak_height) < 0:
                    raise ValueError(msg.error("Error: Negative min_peak_height given."))

        if is_empty(peak_height_operator):
            msg.warning("No peak_height_operator specified, using default: max.")
        else:
            if peak_height_operator not in ["sum", "max"]:
                raise ValueError(msg.error(f"Error: Given peak_height_operator does not exist: {peak_height_operator}. Use [sum, max]"))

        if is_empty(tts_start_selection):
            msg.warning("No tts_start_selection specified, using default: furthest_inframe.")
        else:
            if tts_start_selection not in ["furthest_inframe", "next_inframe"]:
                raise ValueError(msg.error(f"Error: Given tts_start_selection does not exist {tts_start_selection}. Use [furthest_inframe, next_inframe]"))

        if is_empty(rpkm_read_usage):
            msg.warning("No rpkm_read_usage specified, using default: all.")
        else:
            if rpkm_read_usage not in ["all", "specific"]:
                raise ValueError(msg.error(f"Error: Given rpkm_read_usage does not exist {rpkm_read_usage}. Use [all, specific]"))

        if is_empty(gff_output_mode):
            msg.warning("No gff_output_mode specified, using default: combined.")
        else:
            if gff_output_mode not in ["combined", "split"]:
                raise ValueError(msg.error(f"Error: Given gff_output_mode does not exist {gff_output_mode}. Use [combined, split]"))

        if "min" in normalization.lower():
            if is_empty(mapped_counts_file_path):
                raise ValueError(msg.error(
                    "Error: min normalization given but no mapped_counts_file_path specified.\n"
                    "       Please specify a file or choose a different normalization.\n"
                    "       You can use our helper script to create the file."
                ))

    return config_df

def retrieve_bam_input_information(
    file_path_tis: Path | None,
    file_path_tts: Path | None,
    file_path_ribo: Path | None
) -> list[tuple[list[Path | None], str]]:
    """
    Create a list of matching TIS/TTS/RIBO condition+replicate files to run together.

    Args:
        file_path_tis: Path to TIS alignment files directory (None if not used)
        file_path_tts: Path to TTS alignment files directory (None if not used)
        file_path_ribo: Path to RIBO alignment files directory (None if not used)

    Returns:
        List of tuples containing ([tis_file, tts_file, ribo_file], wildcard)
        where each file is a Path object (None if not present for that sample)
    """
    def get_bam_files(dir_path: Path | None, method_prefix: str) -> list[tuple[Path, str]]:
        """Get BAM files from directory matching method prefix, excluding RNA files."""
        if not dir_path:
            return []
        return [
            (dir_path / f.name, f.stem.split('_')[0])
            for f in dir_path.iterdir()
            if f.is_file() and f.suffix == ".bam" and method_prefix in f.name and "RNA" not in f.name
        ]

    # Collect all BAM files with their method info
    tis_files = get_bam_files(file_path_tis, "TIS")
    tts_files = get_bam_files(file_path_tts, "TTS")
    ribo_files = get_bam_files(file_path_ribo, "RIBO")

    # Group files by (condition, replicate)
    sample_dict: dict[tuple[str, str], list[Path | None]] = {}

    for file_path, wildcard in tis_files + tts_files + ribo_files:
        method, condition, replicate = wildcard.split("-")
        key = (condition, replicate)

        if key not in sample_dict:
            sample_dict[key] = [None, None, None]

        method_index = {"TIS": 0, "TTS": 1, "RIBO": 2}[method]
        sample_dict[key][method_index] = Path(file_path)

    return [
        (files, f"{condition}-{replicate}")
        for (condition, replicate), files in sample_dict.items()
        if any(files)
    ]


def call_orfbounder(config_df: pd.DataFrame, result_path: Path) -> None:
    """
    Run the ORFBounder experiments specified in the config sheet.

    Args:
        config_df: DataFrame containing validated configuration
        result_path: Base path for all results
    """

    for row in config_df.itertuples(index=False, name="Pandas"):
        # Required
        experiment = getattr(row, "experiment_name")
        annotation = Path(getattr(row, "annotation_file_path"))
        genome = Path(getattr(row, "genome_file_path"))
        file_path_tis = getattr(row, "TIS_folder_path")
        file_path_tts = getattr(row, "TTS_folder_path")
        file_path_ribo = getattr(row, "RIBO_folder_path")
        mapping_method = getattr(row, "mapping_method").split(",")
        normalization = getattr(row, "normalization_method").split(",")
        offset_json = getattr(row, "offset_file_path")

        # Optional
        read_length_json = Path(getattr(row, "read_length_json")) if not is_empty(getattr(row, "read_length_json")) else None
        bam_folder = Path(getattr(row, "alignment_folder_path")) if not is_empty(getattr(row, "alignment_folder_path")) else None
        min_peak_height = getattr(row, "min_peak_height")
        peak_height_operator = getattr(row, "peak_height_operator")
        tts_start_selection = getattr(row, "tts_start_selection")
        rpkm_read_usage = getattr(row, "rpkm_read_usage")
        gff_output_mode = getattr(row, "gff_output_mode")
        start_codons = getattr(row, "start_codons")
        stop_codons = getattr(row, "stop_codons")
        mapped_counts_file_path = Path(getattr(row, "mapped_counts_file_path")) if not is_empty(getattr(row, "mapped_counts_file_path")) else None

        if is_empty(min_peak_height):
            min_peak_height = 5
        else:
            min_peak_height = int(min_peak_height)

        if is_empty(peak_height_operator):
            peak_height_operator = "max"

        if is_empty(tts_start_selection):
            tts_start_selection = "furthest_inframe"

        if is_empty(rpkm_read_usage):
            all_reads_rpkm = True
        else:
            all_reads_rpkm = rpkm_read_usage != "specific"

        if is_empty(gff_output_mode):
            split_gff = False
        else:
            split_gff = gff_output_mode == "split"

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
                meta_dict: dict = {}
                dynamic_dict: dict = {}
                combined_meta_dict: dict = {}
                combined_dynamic_dict: dict = {}

                # Convert to Path objects, None if empty
                path_tis = Path(file_path_tis) if file_path_tis and not is_empty(file_path_tis) else None
                path_tts = Path(file_path_tts) if file_path_tts and not is_empty(file_path_tts) else None
                path_ribo = Path(file_path_ribo) if file_path_ribo and not is_empty(file_path_ribo) else None

                bam_input_list = retrieve_bam_input_information(path_tis, path_tts, path_ribo)

                res_path = result_path / experiment / mapping / norm

                for (file_tis, file_tts, file_ribo), wildcard in bam_input_list:
                    try:
                        res_df, combined_res_df = ob.run_orfbounder(
                            file_tis,
                            file_tts,
                            file_ribo,
                            read_length_json,
                            norm, mapping, annotation, genome,
                            start_codons, stop_codons, res_path, wildcard,
                            offset_json, tts_start_selection, min_peak_height,
                            peak_height_operator,
                            all_reads_rpkm, bam_folder, mapped_counts_file_path
                        )
                    except SystemExit:
                        msg.warning("Error encountered while calling ORFBounder! Moving to next run!")
                        continue

                    result_tables_path = res_path / RESULT_TABLES_DIR
                    io.write_results_to_gff(res_df, result_tables_path, wildcard, split_gff)
                    io.write_results_to_table(res_df, result_tables_path, wildcard)

                    meta_dict, dynamic_dict = mg.extend_combined_dictionary(res_df, meta_dict, dynamic_dict)

                    if not combined_res_df.empty:
                        combined_results_path = res_path / COMBINED_RESULTS_DIR
                        io.write_results_to_gff(combined_res_df, combined_results_path, wildcard, split_gff)
                        io.write_results_to_table(combined_res_df, combined_results_path, wildcard)
                        combined_meta_dict, combined_dynamic_dict = mg.extend_combined_dictionary(
                            combined_res_df, combined_meta_dict, combined_dynamic_dict
                        )

                if meta_dict:
                    final_excel_path = res_path / f"{experiment}_final{FINAL_OUTPUT_EXTENSION}"
                    result_df = mg.write_merged_table(meta_dict, dynamic_dict, final_excel_path)

                    final_gff_path = res_path / f"{experiment}_final{FINAL_GFF_EXTENSION}"
                    mg.write_merged_gff(result_df, final_gff_path)

                if combined_meta_dict:
                    combined_final_excel_path = res_path / f"{experiment}_combined_final{FINAL_OUTPUT_EXTENSION}"
                    mg.write_merged_table(combined_meta_dict, combined_dynamic_dict, combined_final_excel_path)


def main() -> None:
    """Main entry point for ORFBounder wrapper."""
    parser = argparse.ArgumentParser(
        description="Wrapper for the ORFBounder.py, when running ORFBounder for multiple experiments.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-c", "--config_sheet", action="store", dest="config_sheet", type=Path, required=True,
        help="Config sheet containing information on experiments to be run."
    )
    parser.add_argument(
        "-r", "--result_path", action="store", dest="result_path", type=Path, required=True,
        help="Path of the result folder."
    )

    args = parser.parse_args()

    config_df = check_config_sheet(args.config_sheet)
    call_orfbounder(config_df, args.result_path)


if __name__ == '__main__':
    main()
