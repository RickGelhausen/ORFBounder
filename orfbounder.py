#!/usr/bin/env python

"""
ORFBounder: A peak detection and annotation script for TIS and TTS data.
It can be run with either TIS, TTS or both.
"""

from pathlib import Path

import argparse
import re
import pandas as pd

from lib.alignment_reader import PositionReader

from lib import io
from lib import misc
from lib import predictions as pred
from lib import expression as expr
from lib import messaging as msg


# Constants
RESULT_TABLES_DIR = "result_tables"
OUTPUT_FILE_EXTENSION = ".csv"


def prediction_call(
    annotation_file: Path,
    genome_dict: dict[str, str],
    read_length_dict: dict,
    normalization: str,
    mapping_mode: str,
    start_codons: list[str],
    stop_codons: list[str],
    alignment_file: Path,
    output_path: Path,
    output_basename: str,
    offset_dict: dict,
    method: str,
    longest_potential_orf: str,
    detected_orfs_dict: dict,
    min_peak_height: int,
    peak_height_operator: str,
    min_read_count_dict: dict,
) -> tuple[dict, dict, dict]:
    """
    Execute the script for either TIS or TTS.

    Args:
        annotation_file: Path to annotation file
        genome_dict: Dictionary mapping chromosome names to sequences
        read_length_dict: Dictionary of read length specifications
        normalization: Normalization method to use
        mapping_mode: Read mapping mode (threeprime, fiveprime, centered, global)
        start_codons: List of start codon sequences
        stop_codons: List of stop codon sequences
        alignment_file: Path to alignment file (SAM/BAM)
        output_path: Path to output directory
        output_basename: Base name for output files
        offset_dict: Dictionary of offsets per file/read-length
        method: Analysis method (TIS, TTS, or RIBO)
        longest_potential_orf: ORF selection strategy
        detected_orfs_dict: Dictionary to accumulate detected ORFs
        min_peak_height: Minimum height to consider a peak
        peak_height_operator: Peak height calculation method (max or sum)
        min_read_count_dict: Minimum read count thresholds

    Returns:
        Tuple of (detected_orfs_dict, gene_density_dict, codon_dict)
    """
    if method == "TIS" or method == "RIBO":
        search_codons = start_codons
        match_codons = stop_codons
    else:
        search_codons = stop_codons
        match_codons = start_codons

    msg.message("Checking output folder...")
    result_file = output_path / RESULT_TABLES_DIR / f"{output_basename}{OUTPUT_FILE_EXTENSION}"
    if result_file.is_file():
        raise FileExistsError(msg.error("Error: Result table already found! Please ensure that prior output files with the same name are deleted."))

    msg.success("Done.")

    pr_object = PositionReader(alignment_file, read_length_dict, mapping_mode, offset_dict)
    pr_object.normalize_read_counts(normalization, min_read_count_dict)
    alignment_position_dict, _ = pr_object.output()

    pr_object.to_wig(output_path / "coverage_files")

    for chrom, genome_seq in genome_dict.items():
        msg.message(f"Current chromosome: {chrom}")
        if (chrom, "+") not in alignment_position_dict and (chrom, "-") not in alignment_position_dict:
            msg.warning(f"Warning: No valid entry found for chrom: {chrom}")
            msg.warning("Skipping...")
            continue

        codon_dict, gene_density_dict = {}, {}
        annotation_interlap_dict, gene_density_dict = misc.annotation_interlap(annotation_file)
        gene_density_dict = misc.calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        codon_interlap_dict, codon_dict = misc.create_codon_interlaps(chrom, genome_seq, search_codons)
        codon_dict = pred.screen_positions_for_tss(alignment_position_dict, codon_interlap_dict, codon_dict, min_peak_height, peak_height_operator)

        detected_orfs_dict = pred.detect_potential_orfs(codon_dict, genome_seq, search_codons, match_codons, method, detected_orfs_dict, longest_potential_orf)

    return detected_orfs_dict, gene_density_dict, codon_dict

def run_orfbounder(
    alignment_file_tis: Path | None,
    alignment_file_tts: Path | None,
    alignment_file_ribo: Path | None,
    read_length_json: str,
    normalization: str,
    mapping: str,
    annotation_file_path: Path,
    genome_file_path: Path,
    start_codons: list[str],
    stop_codons: list[str],
    output_path: Path,
    output_basename: str,
    offset_json: str,
    tts_start_selection: str,
    min_peak_height: int,
    peak_height_operator: str,
    all_reads_rpkm: bool,
    alignment_file_path: Path,
    total_read_file_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run functions necessary to generate the final output of ORFBounder.

    Args:
        alignment_file_tis: Path to TIS alignment file
        alignment_file_tts: Path to TTS alignment file
        alignment_file_ribo: Path to RIBO alignment file
        read_length_json: Path to JSON with read length specifications
        normalization: Normalization method
        mapping: Read mapping method
        annotation_file: Path to annotation file
        genome_file: Path to genome sequence file
        start_codons: List of start codons
        stop_codons: List of stop codons
        output_path: Output directory path
        output_basename: Base name for output files
        offset_json: Path to JSON with offset specifications
        tts_start_selection: TTS start codon selection strategy
        min_peak_height: Minimum peak height threshold
        peak_height_operator: Peak height operator (max or sum)
        all_reads_rpkm: Whether to use all reads for RPKM calculation
        alignment_file_path: Path to directory with BAM files for read counting
        total_read_file_path: Path to file with total read counts

    Returns:
        Tuple of (result_df, combined_result_df)
    """


    method = io.parse_alignment_input(alignment_file_tis, alignment_file_tts)
    bam_files = io.check_alignment_path_input(alignment_file_path, alignment_file_tis, alignment_file_tts)
    read_length_dict = io.parse_read_lengths(read_length_json)
    min_read_count_dict = io.parse_total_reads(total_read_file_path, normalization)

    headers = ["TIS", "TTS", "RIBO"]
    if alignment_file_tis:
        headers[0] = re.split(r'_|\.', alignment_file_tis.name)[0]
    if alignment_file_tts:
        headers[1] = re.split(r'_|\.', alignment_file_tts.name)[0]
    if alignment_file_ribo:
        headers[2] = re.split(r'_|\.', alignment_file_ribo.name)[0]

    offset_dict = io.parse_offset_json(offset_json)

    wildcards = []
    if not bam_files:
        msg.message("No valid bam files detected in the bam folder, skipping readcount calculation")
    else:
        for file in bam_files:
            wildcards.append(re.split(r'_|\.', Path(file).name)[0])

        wildcards, bam_files = (list(t) for t in zip(*sorted(zip(wildcards, bam_files))))

    msg.message("Fetching genome...")
    genome_dict = io.generate_genome_dict(genome_file_path)
    msg.success("Done.")

    read_count_dict: dict = {}
    accepted_read_list: list = []
    predictions: dict = {}
    gene_density_tis_dict: dict = {}
    gene_density_tts_dict: dict = {}
    gene_density_ribo_dict: dict = {}
    combined_result_df = pd.DataFrame()

    if method == "TIS":
        predictions, gene_density_tis_dict, _ = prediction_call(
            annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
            start_codons, stop_codons, alignment_file_tis, output_path,
            output_basename, offset_dict, "TIS", tts_start_selection,
            predictions, min_peak_height, peak_height_operator, min_read_count_dict
        )

        if alignment_file_ribo:
            predictions, gene_density_ribo_dict, _ = prediction_call(
                annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
                start_codons, stop_codons, alignment_file_ribo, output_path,
                output_basename, offset_dict, "RIBO", tts_start_selection,
                predictions, min_peak_height, peak_height_operator, min_read_count_dict
            )

        if bam_files:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, accepted_read_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(
            predictions, gene_density_tis_dict, {}, gene_density_ribo_dict, genome_dict, read_count_dict,
            accepted_read_list, wildcards, headers
        )
        msg.success(f"Potential ORFs detected: {len(result_df)}")

    elif method == "TTS":
        predictions, gene_density_tts_dict, _ = prediction_call(
            annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
            start_codons, stop_codons, alignment_file_tts, output_path,
            output_basename, offset_dict, "TTS", tts_start_selection,
            predictions, min_peak_height, peak_height_operator, min_read_count_dict
        )

        if alignment_file_ribo:
            predictions, gene_density_ribo_dict, _ = prediction_call(
                annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
                start_codons, stop_codons, alignment_file_ribo, output_path,
                output_basename, offset_dict, "RIBO", tts_start_selection,
                predictions, min_peak_height, peak_height_operator, min_read_count_dict
            )

        if bam_files:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, accepted_read_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(
            predictions, {}, gene_density_tts_dict, {}, genome_dict, read_count_dict,
            accepted_read_list, wildcards, headers
        )
        msg.success(f"Potential ORFs detected: {len(result_df)}")

    else:  # Both TIS and TTS
        predictions, gene_density_tis_dict, _ = prediction_call(
            annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
            start_codons, stop_codons, alignment_file_tis, output_path,
            output_basename, offset_dict, "TIS", tts_start_selection,
            predictions, min_peak_height, peak_height_operator, min_read_count_dict
        )

        if alignment_file_ribo:
            predictions, gene_density_ribo_dict, _ = prediction_call(
                annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
                start_codons, stop_codons, alignment_file_ribo, output_path,
                output_basename, offset_dict, "RIBO", tts_start_selection,
                predictions, min_peak_height, peak_height_operator, min_read_count_dict
            )

        predictions, gene_density_tts_dict, _ = prediction_call(
            annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
            start_codons, stop_codons, alignment_file_tts, output_path,
            output_basename, offset_dict, "TTS", tts_start_selection,
            predictions, min_peak_height, peak_height_operator, min_read_count_dict
        )

        if bam_files:
            read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
            read_count_dict, accepted_read_list = expr.retrieve_read_counts(read_count_dict, bam_files, read_length_dict, all_reads_rpkm)

        result_df = misc.generate_result_dataframe(
            predictions, gene_density_tis_dict, gene_density_tts_dict, gene_density_ribo_dict, genome_dict,
            read_count_dict, accepted_read_list, wildcards, headers
        )
        msg.success(f"Potential ORFs detected: {len(result_df)}")

    return result_df, combined_result_df

def main() -> None:
    """Main entry point for ORFBounder."""
    parser = argparse.ArgumentParser(
        description="ORFBounder is a peak detection and annotation script for TIS and TTS data.\n"
                   +"It can be run with either TIS, TTS or both.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--alignment_file_tis", action="store", dest="alignment_file_tis", type=Path, default=None,
        help="input alignment file for TIS (sam/bam format)."
    )
    parser.add_argument(
        "--alignment_file_tts", action="store", dest="alignment_file_tts", type=Path, default=None,
        help="input alignment file for TTS (sam/bam format)."
    )
    parser.add_argument(
        "--alignment_file_ribo", action="store", dest="alignment_file_ribo", type=Path, default=None,
        help="input alignment file for RIBO (sam/bam format)."
    )

    parser.add_argument(
        "-l", "--read_length_json", action="store", dest="read_length_json", default=None,
        help="JSON file containing read-length specifications per file.\n"
            +"Ranges can be given using the - symbol (e.g. 15-20,31,33-35,39"
    )
    parser.add_argument(
        "-m", "--mapping_method", action="store", dest="mapping", type=str, default="threeprime",
        help="Read-Mapping to be used:\n"
            +"{threeprime, fiveprime, centered, global}\n"
            +"'threeprime': only threeprime end positions of each mapped read are used.\n"
            +"'fiveprime':  only fiveprime end positions of each mapped read are used.\n"
            +"'centered': the three middle nucleotide positions of each mapped read are used.\n"
            +"'global': all positions of each mapped read are used."
    )
    parser.add_argument(
        "-n", "--normalization_method", action="store", dest="normalization", type=str,
        help="Readcount-Normalization methods to be used:\n"
            +"{raw,min,mil}\n"
            +"'raw': unnormalized readcounts\n"
            +"'min': normalized by min #aligned reads / #aligned reads .\n"
            +"'mil': normalized by 1000000 / #aligned reads."
    )

    parser.add_argument(
        "-a", "--annotation_file", action="store", dest="annotation_file", type=Path, required=True,
        help="input annotation file."
    )
    parser.add_argument(
        "-g", "--genome_file", action="store", dest="genome_file", type=Path, required=True,
        help="input sequence file."
    )

    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=["ATG", "GTG", "TTG"])
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=["TAG", "TAA", "TGA"])

    parser.add_argument(
        "--offset_json", action="store", dest="offset_json", type=Path, default=None,
        help="A JSON file containing offsets for each file/read-length combination.\n"
            +"Default value will be used for missing entries."
    )
    parser.add_argument(
        "--peak_height_operator", action="store", dest="peak_height_operator", type=str, default="max",
        help="{max,sum}:\n"
            +"'max': within the codon interval select the highest value (> min_peak_height)\n"
            +"'sum': within the codon interval sum up all values (> min_peak_height)"
    )
    parser.add_argument(
        "--tts_start_selection", action="store", dest="tts_start_selection", type=str, default="furthest_inframe",
        help="{furthest_inframe, next_inframe}\n"
            +"'furthest_inframe': select the furthest inframe start codon that, without overstepping the next inframe stop codon.\n"
            +"'next_inframe': select the closest inframe start codon."
    )
    parser.add_argument(
        "--min_peak_height", action="store", dest="min_peak_height", type=int, default=5,
        help="Minimum height value to be considered a peak. (max option)\n"
            +"Minimum height value to be added to the total peak value (sum option)"
    )
    parser.add_argument(
        "--total_read_file_path", action="store", dest="total_read_file_path", default=None,
        help="A tab seperated file containing the total read lengths for each sample.\n"
            +"This file is only necessary when using the min nomalization method.\n"
            +"We provide a script that creates these files for you."
    )
    parser.add_argument(
        "--all_reads_rpkm", action="store_true", dest="all_reads_rpkm",
        help="If set, all mapped reads will be used for the calculation of RPKM values.\n"
            +"By default only mapped reads of the specified lengths will be used"
    )
    parser.add_argument(
        "--output_basename", action="store", dest="output_basename", type=str, required=True,
        help="the basename for all output files."
    )
    parser.add_argument(
        "--alignment_file_path", action="store", dest="alignment_file_path", type=Path, default=None,
        help="(optional) sam/bam files to calculate RPKM and TE values for the final results."
    )
    parser.add_argument(
        "-o", "--output_path", action="store", dest="output_path", type=Path, required=True,
        help="Output path to the result folder."
    )
    parser.add_argument(
        "--split_gff", action="store_true", dest="split_gff",
        help="Split gff into one for each gene_type."
    )
    args = parser.parse_args()

    result_df, _ = run_orfbounder(
        args.alignment_file_tis, args.alignment_file_tts, args.alignment_file_ribo, args.read_length_json,
        args.normalization, args.mapping, args.annotation_file, args.genome_file, args.start_codons,
        args.stop_codons, args.output_path, args.output_basename, args.offset_json, args.tts_start_selection,
        args.min_peak_height, args.peak_height_operator, args.all_reads_rpkm, args.alignment_file_path,
        args.total_read_file_path
    )

    output_result_dir = Path(args.output_path) / RESULT_TABLES_DIR
    io.write_results_to_gff(result_df, output_result_dir, args.output_basename, args.split_gff)
    io.write_results_to_table(result_df, output_result_dir, args.output_basename)

    msg.success("Success! Terminating...")


if __name__ == '__main__':
    main()
