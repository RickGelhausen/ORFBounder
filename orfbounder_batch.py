#!/usr/bin/env python
"""
A wrapper script to call ORFBounder for multiple experiments based on a config sheet.
"""

from pathlib import Path

import argparse
import math

import pandas as pd

import orfbounder as ob
import lib.merging as mg
import lib.messaging as msg
from lib import run_output
from lib.alignment_reader import AlignmentPolicy
from lib.config import (
    check_config_sheet, discover_samples as retrieve_bam_input_information,
    parse_sample_name,
)

# Constants
RESULT_TABLES_DIR = "table_per_sample"
RESULT_GFF_DIR = "gff_per_sample"
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

def call_orfbounder(config_df: pd.DataFrame, result_path: Path) -> None:
    """Analyze validated experiment rows and publish each complete result set."""
    from lib.reporting import write_report

    result_path = Path(result_path)
    config_snapshot = config_df.attrs.get("source_snapshot", [])
    run_output.verify_inputs_unchanged(config_snapshot)
    planned = [result_path / row.experiment_name / mapping / norm
               for row in config_df.itertuples(index=False)
               for mapping in row.mapping_method.split(",")
               for norm in row.normalization_method.split(",")]
    for path in planned:
        run_output.assert_output_directory_available(path)
    fingerprint_cache = {}
    additional_inputs = (config_df.attrs["source_path"],) if "source_path" in config_df.attrs else ()
    for row in config_df.itertuples(index=False):
        reference_cache = {}
        comparison_path = (
            Path(row.matched_comparison_file_path) if row.matched_comparison_file_path else None
        )
        if comparison_path:
            from lib.matched_statistics import read_matched_comparisons
            comparison_snapshot = run_output.fingerprint_inputs(
                [comparison_path], cache=fingerprint_cache,
            )
            comparisons = read_matched_comparisons(comparison_path)
            run_output.verify_inputs_unchanged(comparison_snapshot)
        else:
            comparisons = []
            comparison_snapshot = []
        row_additional_inputs = additional_inputs
        row_additional_input_records = tuple(comparison_snapshot)
        alignment_policy = AlignmentPolicy(
            min_mapq=int(row.min_mapq),
            include_duplicates=row.duplicates == "include",
            multimappers=row.multimappers,
        )
        def optional_path(column):
            value = getattr(row, column)
            return Path(value) if not is_empty(value) else None

        groups = retrieve_bam_input_information(*(optional_path(column) for column in
                    ("TIS_folder_path", "TTS_folder_path", "RIBO_folder_path")))
        normalization_samples = tuple(sorted({
            path.stem.split("_", 1)[0]
            for group, _ in groups for path in group if path
        }))
        for mapping in row.mapping_method.split(","):
            for norm in row.normalization_method.split(","):
                destination = result_path / row.experiment_name / mapping / norm
                with run_output.staged_output(destination) as staging:
                    resolved = row._asdict()
                    resolved["mapping_method"], resolved["normalization_method"] = mapping, norm
                    pd.DataFrame([resolved]).to_csv(staging / "config.resolved.tsv", sep="\t", index=False)
                    meta_dict, dynamic_dict, manifests, sample_statistics = {}, {}, [], {}
                    for (tis, tts, ribo), basename in groups:
                        before = set(staging.rglob("*"))
                        result, candidates = ob.run_orfbounder(
                            tis, tts, ribo, optional_path("read_length_json"), norm, mapping,
                            Path(row.annotation_file_path), Path(row.genome_file_path),
                            row.start_codons.split(","), row.stop_codons.split(","), staging, basename,
                            Path(row.offset_file_path), row.tts_start_selection, int(row.min_peak_height),
                            row.peak_height_operator, row.rpkm_read_usage == "all",
                            optional_path("alignment_folder_path"), optional_path("mapped_counts_file_path"),
                            statistics=row.statistics, background_width=int(row.background_width),
                            fdr=float(row.fdr), fdr_method=row.fdr_method,
                            genetic_code=int(row.genetic_code), normalization_scope=row.normalization_scope,
                            alignment_policy=alignment_policy,
                            _fingerprint_cache=fingerprint_cache, _additional_inputs=row_additional_inputs,
                            _additional_input_records=row_additional_input_records,
                            _reference_cache=reference_cache,
                            _normalization_samples=normalization_samples,
                        )
                        if comparisons:
                            requested_assays = {comparison.assay for comparison in comparisons}
                            for method, alignment in zip(("TIS", "TTS", "RIBO"), (tis, tts, ribo)):
                                if alignment is None or method not in requested_assays:
                                    continue
                                _, condition, replicate = parse_sample_name(alignment)
                                sample_statistics[(method, condition, replicate)] = (
                                    candidates.loc[candidates["assay"] == method]
                                    .drop(columns=["assay", "sample_name"])
                                    .reset_index(drop=True)
                                )
                        run_output.export_sample_results(result, staging, basename, row.gff_output_mode == "split")
                        files = [path for path in staging.rglob("*") if path.is_file() and path not in before]
                        manifests.append(run_output.complete_run(staging, basename, files))
                        meta_dict, dynamic_dict = mg.extend_combined_dictionary(result, meta_dict, dynamic_dict)
                    final_basename = f"{row.experiment_name}_final"
                    matched_results = []
                    if comparisons:
                        from lib.matched_statistics import (
                            MATCHED_ALTERNATIVE,
                            MATCHED_MODEL_NAME,
                            MATCHED_NUMERIC_BACKEND,
                            MATCHED_SUPPORT_RULE,
                            MAX_EXACT_CONVOLUTION_WORK,
                            MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL,
                            MIN_INFORMATIVE_PAIRS,
                            MIN_POSITIVE_P_VALUE,
                            analyze_matched_condition,
                        )
                        matched_path = staging / "matched_statistics"
                        matched_path.mkdir(parents=True, exist_ok=True)
                        comparison_metadata = []
                        for comparison in comparisons:
                            statistics, strata = analyze_matched_condition(
                                sample_statistics, comparison.assay,
                                comparison.numerator_condition, comparison.denominator_condition,
                                fdr=float(row.matched_fdr), correction=row.matched_fdr_method,
                                return_strata=True,
                            )
                            filename = f"{comparison.prefix}.tsv"
                            strata_filename = f"{comparison.prefix}_strata.tsv"
                            statistics.to_csv(matched_path / filename, sep="\t", index=False)
                            strata.to_csv(matched_path / strata_filename, sep="\t", index=False)
                            matched_results.append((comparison.prefix, comparison.assay, statistics))
                            replicate_ids = sorted({
                                replicate
                                for assay, condition, replicate in sample_statistics
                                if assay == comparison.assay
                                and condition == comparison.numerator_condition
                            })
                            comparison_metadata.append({
                                "comparison": comparison.name, "prefix": comparison.prefix,
                                "assay": comparison.assay,
                                "numerator_condition": comparison.numerator_condition,
                                "denominator_condition": comparison.denominator_condition,
                                "fdr": float(row.matched_fdr), "fdr_method": row.matched_fdr_method,
                                "model": MATCHED_MODEL_NAME,
                                "numeric_backend": MATCHED_NUMERIC_BACKEND,
                                "alternative": MATCHED_ALTERNATIVE,
                                "replicate_ids": replicate_ids,
                                "minimum_informative_pairs": MIN_INFORMATIVE_PAIRS,
                                "support_rule": MATCHED_SUPPORT_RULE,
                                "max_exact_convolution_work": MAX_EXACT_CONVOLUTION_WORK,
                                "maximum_supported_hypergeometric_total": (
                                    MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL
                                ),
                                "minimum_positive_p_value": MIN_POSITIVE_P_VALUE,
                                "candidate_table": f"matched_statistics/{filename}",
                                "strata_table": f"matched_statistics/{strata_filename}",
                            })
                        for manifest in manifests:
                            manifest["matched_comparisons"] = comparison_metadata
                            run_output.write_json(
                                staging / f"{manifest['output_basename']}.run.json", manifest,
                            )
                    result = mg.write_merged_table(
                        meta_dict, dynamic_dict, staging / f"{final_basename}.xlsx",
                        matched_results=matched_results,
                    )
                    mg.write_merged_gff(result, staging / f"{final_basename}.gff")
                    run_output.write_sequences(result, staging, final_basename)
                    links = run_output.report_links(staging, final_basename, merged=True)
                    file_index = run_output.write_file_index(staging, "files.html", "Coverage and statistical tables",
                        [*staging.rglob("*.wig"), *staging.glob("peak_statistics/*.tsv"),
                         *staging.glob("matched_statistics/*.tsv")])
                    links["Coverage and statistical tables"] = file_index.name
                    write_report(result, staging / "report.html", f"ORFBounder results: {row.experiment_name}",
                                 run_metadata=manifests, links=links)
                    run_output.verify_inputs_unchanged(config_snapshot)
                    run_output.write_completion(staging, "complete.json", manifests)
                msg.success(f"Analysis complete. Open {destination / 'report.html'}")


def main() -> None:
    """Main entry point for ORFBounder wrapper."""
    parser = argparse.ArgumentParser(
        description="Analyze experiments from a TSV configuration and create reports, tables, sequences and genome-browser files.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-c", "--config_sheet", action="store", dest="config_sheet", type=Path, required=True,
        help="Config sheet containing information on experiments to be run."
    )
    parser.add_argument(
        "-r", "--result_path", action="store", dest="result_path", type=Path,
        help="Path of the result folder."
    )

    parser.add_argument("--validate-only", action="store_true", help="Check all configuration rows and sample discovery without running analysis.")
    parser.add_argument("--version", action="version", version="ORFBounder 2.0.0")
    args = parser.parse_args()
    if not args.validate_only and args.result_path is None:
        parser.error("-r/--result_path is required unless --validate-only is used.")
    try:
        config_df = check_config_sheet(args.config_sheet)
        if args.validate_only:
            msg.success(f"Configuration valid: {len(config_df)} experiment(s).")
            return
        call_orfbounder(config_df, args.result_path)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"ORFBounder error: {exc}\n")


if __name__ == '__main__':
    main()
