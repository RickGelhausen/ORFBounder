#!/usr/bin/env python

"""
ORFBounder: A peak detection and annotation script for TIS and TTS data.
It can be run with either TIS, TTS or both.
"""

from pathlib import Path

import argparse
import re
import json
import math
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from numbers import Integral, Real
import pandas as pd

from lib.alignment_reader import AlignmentPolicy, PositionReader, validate_normalization_scope
from lib.genetic_code import resolve_genetic_code
from lib import run_output
from lib.config import parse_sample_name

from lib import io
from lib import misc
from lib import predictions as pred
from lib import expression as expr
from lib import messaging as msg


# Constants
RESULT_TABLES_DIR = "table_per_sample"
RESULT_GFF_DIR = "gff_per_sample"
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
    statistics_results: dict | None = None,
    background_width: int = 50,
    fdr: float = 0.05,
    fdr_method: str = "by",
    diagnostics: dict | None = None,
    normalization_scope: str = "contig",
    alignment_policy: AlignmentPolicy | None = None,
    circular_contigs: set[str] | None = None,
    annotation_features: tuple[misc.AnnotationFeature, ...] | None = None,
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
    circular_contigs = set(circular_contigs or ())

    msg.message("Checking output folder...")
    result_file = output_path / RESULT_TABLES_DIR / f"{output_basename}{OUTPUT_FILE_EXTENSION}"
    if result_file.is_file():
        raise FileExistsError("Result table already found. Choose a new output path or basename.")

    msg.success("Done.")

    pr_object = PositionReader(
        alignment_file, read_length_dict, mapping_mode, offset_dict, alignment_policy,
        circular_contigs,
    )
    if diagnostics is not None:
        diagnostics[alignment_file.stem] = {
            **pr_object.alignment_diagnostics,
            "accepted_reads_by_contig": pr_object.no_accepted_reads_dict,
        }
    if statistics_results is not None and method in {"TIS", "TTS"}:
        from lib.statistics import analyze_peak_enrichment
        statistics_results[method] = analyze_peak_enrichment(
            genome_dict, pr_object.reads_position_dict, search_codons,
            flank_width=background_width, fdr=fdr, correction=fdr_method,
            circular_contigs=circular_contigs,
        )
    pr_object.normalize_read_counts(normalization, min_read_count_dict, normalization_scope=normalization_scope)
    alignment_position_dict, _ = pr_object.output()

    pr_object.to_wig(output_path / "coverage_files" / output_basename)

    annotation_interlap_dict, gene_density_dict = misc.annotation_interlap(
        annotation_file, annotation_features=annotation_features,
    )
    gene_density_dict = misc.calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)
    all_codons = {}
    for chrom, genome_seq in genome_dict.items():
        msg.message(f"Current chromosome: {chrom}")
        if (chrom, "+") not in alignment_position_dict and (chrom, "-") not in alignment_position_dict:
            msg.warning(f"Warning: No valid entry found for chrom: {chrom}")
            msg.warning("Skipping...")
            continue

        is_circular = chrom in circular_contigs
        codon_interlap_dict, codon_dict = misc.create_codon_interlaps(
            chrom, genome_seq, search_codons, circular=is_circular,
        )
        codon_dict = pred.screen_positions_for_tss(alignment_position_dict, codon_interlap_dict, codon_dict, min_peak_height, peak_height_operator)

        detected_orfs_dict = pred.detect_potential_orfs(
            codon_dict, genome_seq, search_codons, match_codons, method,
            detected_orfs_dict, longest_potential_orf, circular=is_circular,
        )
        all_codons.update(codon_dict)

    return detected_orfs_dict, gene_density_dict, all_codons

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
    offset_json: Path,
    tts_start_selection: str,
    min_peak_height: int,
    peak_height_operator: str,
    all_reads_rpkm: bool,
    alignment_file_path: Path,
    total_read_file_path: Path,
    statistics: str = "none",
    background_width: int = 50,
    fdr: float = 0.05,
    fdr_method: str = "by",
    genetic_code: int = 11,
    normalization_scope: str = "contig",
    _fingerprint_cache: dict | None = None,
    _additional_inputs: tuple = (),
    _additional_input_records: tuple = (),
    alignment_policy: AlignmentPolicy | None = None,
    _reference_cache: dict | None = None,
    _normalization_samples: tuple[str, ...] = (),
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
        Tuple of the ORF table and a long table of sequence-defined statistical
        candidates (empty when local statistics are disabled).
    """


    validate_normalization_scope(normalization_scope, normalization)
    alignment_policy = alignment_policy or AlignmentPolicy()
    if not isinstance(alignment_policy, AlignmentPolicy):
        raise TypeError("alignment_policy must be an AlignmentPolicy instance.")
    io.parse_alignment_input(alignment_file_tis, alignment_file_tts)
    if alignment_file_ribo and not alignment_file_ribo.is_file():
        raise FileNotFoundError(f"RIBO alignment does not exist: {alignment_file_ribo}")
    sample_identities = []
    for expected_method, alignment in (
            ("TIS", alignment_file_tis), ("TTS", alignment_file_tts),
            ("RIBO", alignment_file_ribo)):
        if alignment is None:
            continue
        if Path(alignment).suffix.lower() not in {".sam", ".bam"}:
            raise ValueError(
                f"{expected_method} alignment must be a .sam or .bam file: {alignment}"
            )
        try:
            method, condition, replicate = parse_sample_name(alignment)
        except ValueError:
            # Direct runs historically allow custom sample labels. Validate
            # assay and sample identity wherever the documented pattern is
            # available, without changing that custom-label behavior.
            continue
        if method != expected_method:
            raise ValueError(
                f"{expected_method} input has {method} sample identity: {Path(alignment).name}"
            )
        sample_identities.append((condition, replicate))
    if len(set(sample_identities)) > 1:
        names = ", ".join(
            Path(path).stem.split("_", 1)[0]
            for path in (alignment_file_tis, alignment_file_tts, alignment_file_ribo)
            if path is not None
        )
        raise ValueError(
            "Direct-run TIS, TTS and RIBO inputs must share one condition and "
            f"replicate; got {names}."
        )
    if mapping not in {"fiveprime", "threeprime", "centered", "global"}:
        raise ValueError(f"Invalid mapping method: {mapping}")
    if normalization not in {"raw", "mil", "min"}:
        raise ValueError(f"Invalid normalization method: {normalization}")
    if (isinstance(min_peak_height, bool)
            or not isinstance(min_peak_height, Integral)
            or min_peak_height < 0):
        raise ValueError("min_peak_height must be a nonnegative integer.")
    if peak_height_operator not in {"max", "sum"}:
        raise ValueError("peak_height_operator must be max or sum.")
    if tts_start_selection not in {"furthest_inframe", "next_inframe"}:
        raise ValueError("tts_start_selection must be furthest_inframe or next_inframe.")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", output_basename) is None:
        raise ValueError("output_basename must be a filename without path separators, starting with a letter or digit.")
    if statistics not in {"none", "local"}:
        raise ValueError("statistics must be none or local.")
    if statistics == "local" and mapping not in {"fiveprime", "threeprime"}:
        raise ValueError("Local statistics require fiveprime or threeprime mapping (integer endpoint counts).")
    if (isinstance(fdr, bool) or not isinstance(fdr, Real)
            or not math.isfinite(fdr) or not 0 < fdr < 1
            or isinstance(background_width, bool)
            or not isinstance(background_width, Integral)
            or background_width < 1 or fdr_method not in {"bh", "by"}):
        raise ValueError("Statistical settings require 0 < fdr < 1, background_width >= 1, and fdr_method bh or by.")
    start_codons, stop_codons = resolve_genetic_code(genetic_code, start_codons, stop_codons)
    if read_length_json is not None and not Path(read_length_json).is_file():
        raise FileNotFoundError(f"Read-length JSON does not exist: {read_length_json}")
    if alignment_file_path is not None and not Path(alignment_file_path).is_dir():
        raise FileNotFoundError(f"Expression alignment folder does not exist: {alignment_file_path}")
    output_path = Path(output_path)
    run_output.assert_sample_output_available(output_path, output_basename)

    bam_files = io.check_alignment_path_input(alignment_file_path, alignment_file_tis, alignment_file_tts)
    alignments = {"TIS": alignment_file_tis, "TTS": alignment_file_tts, "RIBO": alignment_file_ribo}
    headers = [path.stem.split("_")[0] if path else method for method, path in alignments.items()]
    if len(set(headers)) != len(headers):
        raise ValueError("TIS, TTS and RIBO inputs must have distinct sample names. Use METHOD-condition-replicate filenames.")
    bam_files = sorted(bam_files or [], key=lambda path: path.stem)
    wildcards = [path.stem.split("_")[0] for path in bam_files]
    input_files = run_output.merge_input_records(
        run_output.fingerprint_inputs(
            [*alignments.values(), *bam_files, genome_file_path, annotation_file_path,
             read_length_json, offset_json, total_read_file_path, *_additional_inputs],
            cache=_fingerprint_cache,
        ),
        _additional_input_records,
    )
    read_length_dict = io.parse_read_lengths(Path(read_length_json) if read_length_json is not None else None)
    filtered_alignments = [path for path in alignments.values() if path]
    if not all_reads_rpkm:
        filtered_alignments.extend(bam_files)
    io.validate_read_length_samples(read_length_dict, filtered_alignments)
    # Batch callers add the other samples which share this normalization
    # target, but may never replace the alignments in the current call.  Keeping
    # the union here makes the private batching hint fail closed if it is stale
    # or incomplete.
    normalization_samples = {
        path.stem.split("_", 1)[0] for path in alignments.values() if path
    }
    normalization_samples.update(map(str, _normalization_samples))
    min_read_count_dict = io.parse_total_reads(
        total_read_file_path, normalization, normalization_scope,
        required_samples=normalization_samples,
    )
    offset_dict = io.parse_offset_json(offset_json)
    io.validate_offset_samples(offset_dict, alignments.values())
    reference_context = io.load_reference_context(
        genome_file_path, annotation_file_path, input_files, cache=_reference_cache,
    )
    genome_dict = reference_context.genome_dict()
    io.validate_reference_inputs(genome_dict, annotation_file_path,
                                 [path for path in alignments.values() if path] + bam_files,
                                 reference_context=reference_context)
    circular_contigs = set(reference_context.circular_contigs)
    reference_lengths = {chrom: len(sequence) for chrom, sequence in genome_dict.items()}

    predictions, densities, diagnostics = {}, {}, {}
    statistics_results = {} if statistics == "local" else None
    for method in ("TIS", "RIBO", "TTS"):
        alignment = alignments[method]
        if alignment is None:
            densities[method] = {}
            continue
        predictions, densities[method], _ = prediction_call(
            annotation_file_path, genome_dict, read_length_dict, normalization, mapping,
            start_codons, stop_codons, alignment, output_path, output_basename, offset_dict,
            method, tts_start_selection, predictions, min_peak_height, peak_height_operator,
            min_read_count_dict, statistics_results, background_width, fdr, fdr_method,
            diagnostics, normalization_scope, alignment_policy, circular_contigs,
            reference_context.annotation_features,
        )
    read_count_dict, accepted_read_list = {}, []
    expression_diagnostics = {}
    if bam_files:
        read_count_dict = expr.init_read_count_dict(read_count_dict, predictions)
        read_count_dict, accepted_read_list = expr.retrieve_read_counts(
            read_count_dict, bam_files, read_length_dict, all_reads_rpkm,
            alignment_policy=alignment_policy, diagnostics=expression_diagnostics,
            reference_lengths=reference_lengths, circular_contigs=circular_contigs,
        )
    result_df = misc.generate_result_dataframe(
        predictions, densities["TIS"], densities["TTS"], densities["RIBO"], genome_dict,
        read_count_dict, accepted_read_list, wildcards, headers,
        genetic_code=genetic_code, normalization_scope=normalization_scope,
        circular_contigs=circular_contigs,
    )
    candidate_frames = []
    if statistics_results is not None:
        from lib.statistics import attach_to_orfs
        statistics_path = output_path / "peak_statistics"
        statistics_path.mkdir(parents=True, exist_ok=True)
        for method, frame in statistics_results.items():
            sample = alignments[method].stem.split("_")[0]
            frame.to_csv(statistics_path / f"{output_basename}_{sample}.tsv", sep="\t", index=False)
            result_df = attach_to_orfs(
                result_df, frame, method, sample,
                reference_lengths=reference_lengths, circular_contigs=circular_contigs,
            )
            candidate_frames.append(frame.assign(assay=method, sample_name=sample))
    candidate_results = (pd.concat(candidate_frames, ignore_index=True)
                         if candidate_frames else pd.DataFrame())
    statistics_metadata = {}
    if statistics_results is not None:
        from lib.statistics import (
            LOCAL_ALTERNATIVE,
            LOCAL_MODEL_NAME,
            MAX_EXACT_COMBINATION_WORK,
            MAX_SUPPORTED_TOTAL_COUNT,
            MAX_TAIL_RECURRENCE_WORK,
            MIN_POSITIVE_P_VALUE,
        )
        statistics_metadata = {
            "statistics_model": LOCAL_MODEL_NAME,
            "statistics_alternative": LOCAL_ALTERNATIVE,
            "statistics_minimum_positive_p_value": MIN_POSITIVE_P_VALUE,
            "statistics_maximum_supported_total_count": MAX_SUPPORTED_TOTAL_COUNT,
            "statistics_max_tail_recurrence_work": MAX_TAIL_RECURRENCE_WORK,
            "statistics_max_exact_combination_work": MAX_EXACT_COMBINATION_WORK,
        }
    output_path.mkdir(parents=True, exist_ok=True)
    manifest = {
        "orfbounder_version": "2.0.0", "python_version": platform.python_version(),
        "dependencies": {package: version(package) for package in ("pysam", "pandas", "numpy", "biopython", "scipy")},
        "status": "analysis_complete", "analyzed_at": datetime.now(timezone.utc).isoformat(), "output_basename": output_basename,
        "genetic_code": genetic_code, "normalization_scope": normalization_scope, "input_files": input_files,
        "inputs": {method: str(path.resolve()) if path else None for method, path in alignments.items()},
        "sample_names": {method: path.stem.split("_")[0] for method, path in alignments.items() if path},
        "genome": str(Path(genome_file_path).resolve()), "annotation": str(Path(annotation_file_path).resolve()),
        "reference_lengths": reference_lengths, "circular_contigs": sorted(circular_contigs),
        "mapping": mapping, "normalization": normalization, "offsets": offset_dict,
        "read_lengths": read_length_dict, "minimum_counts": min_read_count_dict,
        "start_codons": start_codons, "stop_codons": stop_codons,
        "tts_start_selection": tts_start_selection, "min_peak_height": min_peak_height,
        "peak_height_operator": peak_height_operator, "all_reads_rpkm": all_reads_rpkm,
        "expression_alignments": [str(path.resolve()) for path in bam_files],
        "statistics": statistics, "background_width": background_width, "fdr": fdr,
        "fdr_method": fdr_method,
        **statistics_metadata,
        "alignment_policy": {
            "min_mapq": alignment_policy.min_mapq,
            "duplicates": "include" if alignment_policy.include_duplicates else "exclude",
            "multimappers": alignment_policy.multimappers,
            "count_unit": "alignment_record",
        },
        "alignment_diagnostics": diagnostics,
        "expression_alignment_diagnostics": expression_diagnostics,
        "orf_count": len(result_df),
    }
    run_output.verify_inputs_unchanged(input_files)
    run_output.write_json(output_path / f"{output_basename}.run.json", manifest)
    msg.success(f"Potential ORFs detected: {len(result_df)}")
    return result_df, candidate_results

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
            +"Use lengths/ranges such as 28-30,32. Omit this file to use all lengths."
    )
    parser.add_argument(
        "-m", "--mapping_method", action="store", dest="mapping", type=str, default="threeprime", choices=["threeprime", "fiveprime", "centered", "global"],
        help="Read-Mapping to be used:\n"
            +"{threeprime, fiveprime, centered, global}\n"
            +"'threeprime': only threeprime end positions of each mapped read are used.\n"
            +"'fiveprime':  only fiveprime end positions of each mapped read are used.\n"
            +"'centered': trim 11 query bases from each end and share one count over the remaining aligned bases.\n"
            +"'global': all aligned reference positions of each mapped read are used."
    )
    parser.add_argument(
        "-n", "--normalization_method", action="store", dest="normalization", type=str, default="raw", choices=["raw", "mil", "min"],
        help="Readcount-Normalization methods to be used:\n"
            +"{raw,min,mil}\n"
            +"'raw': unnormalized readcounts\n"
            +"'min': scale to the smallest supplied accepted-read count.\n"
            +"'mil': scale to one million accepted reads."
    )

    parser.add_argument(
        "-a", "--annotation_file", action="store", dest="annotation_file", type=Path, required=True,
        help="CDS annotation in GFF3 format, from the same reference assembly as the alignments."
    )
    parser.add_argument(
        "-g", "--genome_file", action="store", dest="genome_file", type=Path, required=True,
        help="Reference genome in FASTA format."
    )

    parser.add_argument("--start_codons", nargs="+", dest="start_codons", default=None, help="Optional start codons allowed by the selected genetic code.")
    parser.add_argument("--stop_codons", nargs="+", dest="stop_codons", default=None, help="Optional stop codons allowed by the selected genetic code.")

    parser.add_argument(
        "--offset_json", action="store", dest="offset_json", type=Path, required=True,
        help="A JSON file containing offsets for each file/read-length combination.\n"
            +"An explicitly supplied default is used when a sample or read-length entry is missing."
    )
    parser.add_argument(
        "--peak_height_operator", action="store", dest="peak_height_operator", type=str, default="max", choices=["max", "sum"],
        help="{max,sum}:\n"
            +"'max': within the codon interval select the highest value (> min_peak_height)\n"
            +"'sum': within the codon interval sum up all values (> min_peak_height)"
    )
    parser.add_argument(
        "--tts_start_selection", action="store", dest="tts_start_selection", type=str, default="furthest_inframe", choices=["furthest_inframe", "next_inframe"],
        help="{furthest_inframe, next_inframe}\n"
            +"'furthest_inframe': select the furthest in-frame start without crossing another in-frame stop.\n"
            +"'next_inframe': select the closest inframe start codon."
    )
    parser.add_argument(
        "--min_peak_height", action="store", dest="min_peak_height", type=int, default=5,
        help="Only positions strictly greater than this height contribute to a called peak (default: 5)."
    )
    parser.add_argument(
        "--total_read_file_path", action="store", dest="total_read_file_path", default=None,
        help="Headerless tab-separated accepted read counts: sample, chromosome, count.\n"
            +"Required only for min normalization."
    )
    parser.add_argument(
        "--all_reads_rpkm", action="store_true", dest="all_reads_rpkm",
        help="If set, all mapped reads will be used for the calculation of RPKM values.\n"
            +"By default only mapped reads of the specified lengths will be used"
    )
    parser.add_argument(
        "--output_basename", action="store", dest="output_basename", type=str, required=True,
        help="Sample/run name used for its report, completion record and result files."
    )
    parser.add_argument(
        "--alignment_file_path", action="store", dest="alignment_file_path", type=Path, default=None,
        help="Optional folder of matched SAM/BAM files for RPKM and translation efficiency."
    )
    parser.add_argument(
        "-o", "--output_path", action="store", dest="output_path", type=Path, required=True,
        help="Output path to the result folder."
    )
    parser.add_argument(
        "--split_gff", action="store_true", dest="split_gff",
        help="Also write a GFF file for each annotation category."
    )
    parser.add_argument("--genetic-code", type=int, default=11, help="NCBI genetic-code number (default: 11, bacterial).")
    parser.add_argument(
        "--normalization-scope", choices=["contig", "library"], default="contig",
        help=("Choose per-contig or whole-library normalization denominators. "
              "For min, library uses the smallest summed sample total (default: contig)."),
    )
    parser.add_argument("--statistics", choices=["none", "local"], default="none", help="Optional exploratory local peak enrichment (endpoint mapping only).")
    parser.add_argument("--background-width", type=int, default=50, help="Number of reference bases in each local background flank (default: 50).")
    parser.add_argument("--fdr", type=float, default=0.05, help="Adjusted p-value threshold for support flags (default: 0.05).")
    parser.add_argument("--fdr-method", choices=["by", "bh"], default="by", help="Multiple-testing correction; BY supports dependent tests under the model.")
    parser.add_argument("--min-mapq", type=int, default=0, help="Exclude alignments below this MAPQ (0-255; default: 0).")
    parser.add_argument("--duplicates", choices=["include", "exclude"], default="include", help="Include or exclude duplicate-flagged alignment records (default: include).")
    parser.add_argument("--multimappers", choices=["exclude", "include"], default="exclude", help="Exclude or include records with NH > 1 (default: exclude; records without NH cannot be identified as multimappers).")
    parser.add_argument("--version", action="version", version="ORFBounder 2.0.0")
    args = parser.parse_args()

    try:
        from lib.reporting import write_report
        alignment_policy = AlignmentPolicy(
            min_mapq=args.min_mapq,
            include_duplicates=args.duplicates == "include",
            multimappers=args.multimappers,
        )
        run_output.assert_sample_output_available(args.output_path, args.output_basename)
        with run_output.staged_output(args.output_path) as staging:
            result_df, _ = run_orfbounder(
                args.alignment_file_tis, args.alignment_file_tts, args.alignment_file_ribo, args.read_length_json,
                args.normalization, args.mapping, args.annotation_file, args.genome_file, args.start_codons,
                args.stop_codons, staging, args.output_basename, args.offset_json, args.tts_start_selection,
                args.min_peak_height, args.peak_height_operator, args.all_reads_rpkm, args.alignment_file_path,
                Path(args.total_read_file_path) if args.total_read_file_path else None,
                statistics=args.statistics, background_width=args.background_width,
                fdr=args.fdr, fdr_method=args.fdr_method, genetic_code=args.genetic_code,
                normalization_scope=args.normalization_scope, alignment_policy=alignment_policy,
            )
            basename = args.output_basename
            run_output.export_sample_results(result_df, staging, basename, args.split_gff)
            metadata = json.loads((staging / f"{basename}.run.json").read_text())
            # This report becomes visible only if all exports and publication succeed.
            metadata["status"] = "complete"
            links = run_output.report_links(staging, basename)
            index = run_output.write_file_index(
                staging, f"{basename}.files.html", f"Coverage and statistical tables: {basename}",
                [*staging.rglob("*.wig"), *staging.glob("peak_statistics/*.tsv")])
            links["Coverage and statistical tables"] = index.name
            write_report(result_df, staging / f"{basename}.report.html", f"ORFBounder results: {basename}",
                         run_metadata=[metadata], links=links)
            metadata = run_output.complete_run(staging, basename,
                [path for path in staging.rglob("*") if path.is_file()])
            run_output.write_completion(staging, f"{basename}.complete.json", [metadata])
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"ORFBounder error: {exc}\n")
    msg.success(f"Analysis complete. Open {Path(args.output_path) / (args.output_basename + '.report.html')}")


if __name__ == '__main__':
    main()
