#!/usr/bin/env python3
"""Recover primary mapped read totals for ORFBounder's min normalization.

Run from a checkout with ``python -m helpers.recover_total_counts -c CONFIG -r OUTPUT``.
The count-recovery preflight intentionally does not require an existing mapped
counts file, genome, annotation, or offsets; this allows bootstrapping a batch
configuration that requests min normalization.
"""

from pathlib import Path
import argparse
import re

import pandas as pd

from lib import io, messaging as msg, run_output
from lib.alignment_reader import (
    AlignmentPolicy, alignment_query_length, iter_primary_alignments,
)
from lib.config import discover_samples, sample_name


FOLDER_COLUMNS = ("TIS_folder_path", "TTS_folder_path", "RIBO_folder_path")


def check_config_sheet(config_sheet):
    """Accept a batch TSV or a minimal experiment/folder count-recovery TSV."""
    source_snapshot = run_output.fingerprint_inputs([config_sheet])
    config_df = pd.read_csv(config_sheet, sep="\t", dtype=str, keep_default_na=False)
    config_df = config_df.apply(lambda column: column.str.strip())
    if "experiment_name" not in config_df:
        raise ValueError("Config requires experiment_name and at least one TIS_folder_path or TTS_folder_path column.")
    if config_df.empty:
        raise ValueError("The config sheet contains no experiment rows.")
    defaults = {
        **{column: "" for column in FOLDER_COLUMNS},
        "read_length_json": "", "min_mapq": "0",
        "duplicates": "include", "multimappers": "exclude",
    }
    for column, default in defaults.items():
        if column not in config_df:
            config_df[column] = default
        else:
            config_df.loc[config_df[column] == "", column] = default
    if config_df.experiment_name.duplicated().any():
        raise ValueError("Each experiment_name must be unique.")
    for index, row in config_df.iterrows():
        try:
            if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row.experiment_name) is None:
                raise ValueError("experiment_name must be a safe filename starting with a letter or digit.")
            for column in FOLDER_COLUMNS:
                if row[column] and not Path(row[column]).is_dir():
                    raise ValueError(f"{column}: directory does not exist: {row[column]}")
            if row.read_length_json and not Path(row.read_length_json).is_file():
                raise ValueError(f"Read-length JSON does not exist: {row.read_length_json}")
            read_lengths = io.parse_read_lengths(
                Path(row.read_length_json) if row.read_length_json else None
            )
            if re.fullmatch(r"\d+", row.min_mapq) is None:
                raise ValueError("min_mapq must be an integer from 0 through 255.")
            if row.duplicates not in {"include", "exclude"}:
                raise ValueError("duplicates must be 'include' or 'exclude'.")
            AlignmentPolicy(
                min_mapq=int(row.min_mapq),
                include_duplicates=row.duplicates == "include",
                multimappers=row.multimappers,
            )
            groups = discover_samples(*(
                Path(row[column]) if row[column] else None
                for column in FOLDER_COLUMNS
            ))
            filtered_alignments = [
                path for group, _ in groups for path in group if path is not None
            ]
            io.validate_read_length_samples(read_lengths, filtered_alignments)
        except (OSError, ValueError) as exc:
            raise ValueError(f"Config row {index + 2} ({row.experiment_name}): {exc}") from exc
    config_df.attrs["source_snapshot"] = source_snapshot
    run_output.verify_inputs_unchanged(source_snapshot)
    return config_df


def count_mapped_reads(
    alignment_file_path,
    read_count_dict,
    read_length_dict,
    alignment_policy=None,
    diagnostics=None,
):
    """Count primary mapped records using the same query-length filter as peaks."""
    alignment_file_path = Path(alignment_file_path)
    sample = sample_name(alignment_file_path)
    msg.message(f"Counting reads for {sample}.")
    read_lengths = None
    if read_length_dict:
        read_lengths = read_length_dict.get(sample, read_length_dict.get("default"))
        if read_lengths is None:
            msg.warning(f"No read lengths or default specified for {sample}; using all read lengths.")
    for read in iter_primary_alignments(
            alignment_file_path, alignment_policy=alignment_policy, diagnostics=diagnostics):
        if read_lengths is not None:
            read_length = alignment_query_length(read)
            if str(read_length) not in read_lengths:
                continue
        if not read.get_reference_positions():
            continue
        key = (sample, read.reference_name)
        read_count_dict[key] = read_count_dict.get(key, 0) + 1
    return read_count_dict


def write_read_counts_to_file(
        read_count_dict, experiment_name, result_path, *, announce=True):
    """Write deterministic, headerless sample/contig/count rows."""
    output = Path(result_path) / f"{experiment_name}_mapped_reads.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        for (sample, chrom), count in sorted(read_count_dict.items()):
            handle.write(f"{sample}\t{chrom}\t{count}\n")
    if announce:
        msg.success(f"Wrote {output}")
    return output


def recover_read_information(config_df, result_path):
    """Count only samples that the batch runner will analyze for each experiment."""
    result_path = Path(result_path)
    config_snapshot = config_df.attrs.get("source_snapshot", [])
    run_output.verify_inputs_unchanged(config_snapshot)
    run_output.recover_stale_transactions(result_path)
    for name in config_df.experiment_name:
        output = result_path / f"{name}_mapped_reads.tsv"
        if output.exists() or output.is_symlink():
            raise FileExistsError(f"Output already exists: {output}. Choose a new result path.")
    recovered, input_snapshots = [], [config_snapshot]
    fingerprint_cache = {}
    for row in config_df.itertuples(index=False):
        groups = discover_samples(*(Path(getattr(row, column)) if getattr(row, column) else None
                                    for column in FOLDER_COLUMNS))
        files = sorted({path for group, _ in groups for path in group if path})
        read_length_path = Path(row.read_length_json) if row.read_length_json else None
        input_snapshots.append(run_output.fingerprint_inputs(
            [*files, read_length_path], cache=fingerprint_cache,
        ))
        read_lengths = io.parse_read_lengths(read_length_path)
        read_counts = {}
        policy = AlignmentPolicy(
            min_mapq=int(row.min_mapq),
            include_duplicates=row.duplicates == "include",
            multimappers=row.multimappers,
        )
        for alignment in files:
            count_mapped_reads(alignment, read_counts, read_lengths, policy)
        if not read_counts:
            raise ValueError(f"Experiment {row.experiment_name} has no accepted mapped reads. Check alignments and read lengths; min normalization requires positive counts.")
        recovered.append((row.experiment_name, read_counts))
    run_output.verify_inputs_unchanged(
        run_output.merge_input_records(*input_snapshots), verify_content=True,
    )
    with run_output.staged_output(result_path) as staging:
        for experiment_name, read_counts in recovered:
            write_read_counts_to_file(
                read_counts, experiment_name, staging, announce=False,
            )
    for experiment_name, _ in recovered:
        msg.success(f"Wrote {result_path / f'{experiment_name}_mapped_reads.tsv'}")


def main():
    parser = argparse.ArgumentParser(description="Recover primary mapped read totals for ORFBounder min normalization.")
    parser.add_argument("-c", "--config_sheet", type=Path, required=True, help="Batch TSV or minimal experiment/folder TSV.")
    parser.add_argument("-r", "--result_path", type=Path, required=True, help="Directory for EXPERIMENT_mapped_reads.tsv files.")
    args = parser.parse_args()
    try:
        config_df = check_config_sheet(args.config_sheet)
        recover_read_information(config_df, args.result_path)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"ORFBounder count recovery error: {exc}\n")


if __name__ == "__main__":
    main()
