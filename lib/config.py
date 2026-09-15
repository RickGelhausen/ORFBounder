"""Validation shared by the batch interface and input discovery."""

import math
import re
from pathlib import Path

import pandas as pd

from lib import io, messaging as msg
from lib.run_output import fingerprint_inputs, verify_inputs_unchanged
from lib.genetic_code import resolve_genetic_code
from lib.alignment_reader import AlignmentPolicy, validate_normalization_scope


REQUIRED_COLUMNS = (
    "experiment_name", "annotation_file_path", "genome_file_path", "offset_file_path",
)
DEFAULTS = {
    "TIS_folder_path": "", "TTS_folder_path": "", "RIBO_folder_path": "",
    "alignment_folder_path": "", "normalization_method": "raw",
    "mapping_method": "threeprime", "mapped_counts_file_path": "",
    "read_length_json": "", "min_peak_height": "5", "peak_height_operator": "max",
    "tts_start_selection": "furthest_inframe", "rpkm_read_usage": "all",
    "gff_output_mode": "combined", "start_codons": "",
    "stop_codons": "", "genetic_code": "11", "normalization_scope": "contig", "statistics": "none", "background_width": "50",
    "fdr": "0.05", "fdr_method": "by",
    "min_mapq": "0", "duplicates": "include", "multimappers": "exclude",
    "matched_comparison_file_path": "", "matched_fdr": "0.05",
    "matched_fdr_method": "by",
}
CHOICES = {
    "normalization_method": {"raw", "mil", "min"},
    "mapping_method": {"fiveprime", "threeprime", "centered", "global"},
    "peak_height_operator": {"max", "sum"},
    "tts_start_selection": {"furthest_inframe", "next_inframe"},
    "rpkm_read_usage": {"all", "specific"}, "gff_output_mode": {"combined", "split"},
    "statistics": {"none", "local"}, "fdr_method": {"bh", "by"},
    "normalization_scope": {"contig", "library"},
    "duplicates": {"include", "exclude"}, "multimappers": {"include", "exclude"},
    "matched_fdr_method": {"bh", "by"},
}
SAMPLE_PATTERN = re.compile(
    r"^(TIS|TTS|RIBO|RNA|RNATIS|RNATTS)-([A-Za-z0-9]+)-([1-9][0-9]*)(?:_[^.]+)?$"
)


def parse_sample_name(path: Path) -> tuple[str, str, str]:
    """Parse a sample identity, excluding an optional alignment suffix."""
    match = SAMPLE_PATTERN.fullmatch(Path(path).stem)
    if not match:
        raise ValueError(
            f"Invalid sample filename: {Path(path).name}. Expected "
            "METHOD-condition-replicate.bam (or .sam), e.g. TIS-WT-1.bam; "
            "condition must be alphanumeric and replicate a positive integer."
        )
    return match.group(1, 2, 3)


def sample_name(path: Path) -> str:
    return "-".join(parse_sample_name(path))


def discover_samples(tis: Path | None, tts: Path | None, ribo: Path | None) -> list:
    """Pair samples exactly and reject ambiguous duplicate inputs."""
    samples = {}
    for index, (folder, method) in enumerate(zip((tis, tts, ribo), ("TIS", "TTS", "RIBO"))):
        if folder is None:
            continue
        for path in sorted(Path(folder).iterdir()):
            if not path.is_file() or path.suffix.lower() not in {".sam", ".bam"}:
                continue
            if not path.name.startswith(f"{method}-"):
                continue
            _, condition, replicate = parse_sample_name(path)
            key = (condition, replicate)
            group = samples.setdefault(key, [None, None, None])
            if group[index] is not None:
                raise ValueError(f"Duplicate {method}-{condition}-{replicate} inputs: {group[index]} and {path}")
            group[index] = path
    result = []
    for (condition, replicate), group in sorted(samples.items()):
        if group[0] is None and group[1] is None:
            msg.warning(f"Skipping RIBO-{condition}-{replicate}: no matching TIS or TTS input.")
            continue
        result.append((group, f"{condition}-{replicate}"))
    if not result:
        raise ValueError("No TIS/TTS .bam or .sam samples found. Check folders and METHOD-condition-replicate filenames.")
    return result


def codons(value: str) -> list[str]:
    result = [entry.strip().upper() for entry in value.split(",")]
    if not result or any(re.fullmatch("[ACGT]{3}", entry) is None for entry in result):
        raise ValueError(f"Codons must be comma-separated DNA triplets, e.g. ATG,GTG,TTG; got {value!r}.")
    return list(dict.fromkeys(result))


def nonnegative_integer(value: str) -> int:
    if re.fullmatch(r"\d+", str(value)) is None:
        raise ValueError(f"Expected a nonnegative integer, got {value!r}.")
    return int(value)


def check_config_sheet(config_sheet: Path) -> pd.DataFrame:
    """Load a TSV, insert optional defaults, and validate all runs before output."""
    fingerprint_cache = {}
    reference_cache = {}
    source_snapshot = fingerprint_inputs([config_sheet], cache=fingerprint_cache)
    config_df = pd.read_csv(config_sheet, sep="\t", dtype=str, keep_default_na=False)
    config_df = config_df.apply(lambda col: col.str.strip())
    missing = sorted(set(REQUIRED_COLUMNS) - set(config_df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}. Save the configuration as a tab-separated TSV.")
    if config_df.empty:
        raise ValueError("The config sheet contains no experiment rows.")
    obsolete = {"read_lengths", "max_ORF_length", "log_fold_contrasts"} & set(config_df.columns)
    for column in sorted(obsolete):
        if config_df[column].ne("").any():
            hint = " Use read_length_json pointing to a JSON file." if column == "read_lengths" else " This option is not implemented."
            raise ValueError(f"Unsupported config column {column!r}.{hint}")
    unknown = set(config_df.columns) - set(REQUIRED_COLUMNS) - set(DEFAULTS) - obsolete
    if unknown:
        raise ValueError(f"Unknown config columns: {', '.join(sorted(unknown))}. Check spelling against the template.")
    for column, default in DEFAULTS.items():
        if column not in config_df:
            config_df[column] = default
        else:
            config_df.loc[config_df[column] == "", column] = default
    if config_df.experiment_name.duplicated().any():
        raise ValueError("Each experiment_name must be unique to prevent output collisions.")
    for index, row in config_df.iterrows():
        try:
            if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row.experiment_name) is None:
                raise ValueError("experiment_name must start with a letter or digit and contain only letters, digits, _, - or .")
            for column in (
                    "annotation_file_path", "genome_file_path", "offset_file_path",
                    "read_length_json", "mapped_counts_file_path", "matched_comparison_file_path"):
                value = row[column]
                if column in REQUIRED_COLUMNS or value:
                    if not value or not Path(value).is_file():
                        raise ValueError(f"{column}: file does not exist: {value!r}. Relative paths are resolved from the working directory.")
            for column in ("TIS_folder_path", "TTS_folder_path", "RIBO_folder_path", "alignment_folder_path"):
                if row[column] and not Path(row[column]).is_dir():
                    raise ValueError(f"{column}: directory does not exist: {row[column]}")
            for column, allowed in CHOICES.items():
                values = [part.strip() for part in row[column].split(",")]
                if column not in {"mapping_method", "normalization_method"} and len(values) != 1:
                    raise ValueError(f"{column} accepts a single value.")
                if not set(values) <= allowed:
                    raise ValueError(f"{column}: choose from {', '.join(sorted(allowed))}; got {row[column]!r}.")
                config_df.at[index, column] = ",".join(dict.fromkeys(values))
            nonnegative_integer(row.min_peak_height)
            min_mapq = nonnegative_integer(row.min_mapq)
            AlignmentPolicy(
                min_mapq=min_mapq,
                include_duplicates=row.duplicates == "include",
                multimappers=row.multimappers,
            )
            if nonnegative_integer(row.background_width) < 1:
                raise ValueError("background_width must be at least 1.")
            fdr = float(row.fdr)
            if not math.isfinite(fdr) or not 0 < fdr < 1:
                raise ValueError("fdr must be between 0 and 1, exclusive.")
            matched_fdr = float(row.matched_fdr)
            if not math.isfinite(matched_fdr) or not 0 < matched_fdr < 1:
                raise ValueError("matched_fdr must be between 0 and 1, exclusive.")
            if row.statistics == "local" and set(config_df.at[index, "mapping_method"].split(",")) - {"fiveprime", "threeprime"}:
                raise ValueError("Local statistics require fiveprime or threeprime mapping (integer endpoint counts).")
            if row.matched_comparison_file_path and row.statistics != "local":
                raise ValueError("matched_comparison_file_path requires statistics=local.")
            code = nonnegative_integer(row.genetic_code)
            starts, stops = resolve_genetic_code(
                code, codons(row.start_codons) if row.start_codons else None,
                codons(row.stop_codons) if row.stop_codons else None)
            config_df.at[index, "start_codons"] = ",".join(starts)
            config_df.at[index, "stop_codons"] = ",".join(stops)
            for normalization in config_df.at[index, "normalization_method"].split(","):
                validate_normalization_scope(row.normalization_scope, normalization)
            read_lengths = io.parse_read_lengths(
                Path(row.read_length_json) if row.read_length_json else None
            )
            offsets = io.parse_offset_json(Path(row.offset_file_path))
            groups = discover_samples(*(Path(row[c]) if row[c] else None for c in ("TIS_folder_path", "TTS_folder_path", "RIBO_folder_path")))
            calling_samples = {
                Path(path).stem.split("_", 1)[0]
                for group, _ in groups for path in group if path
            }
            io.parse_total_reads(
                Path(row.mapped_counts_file_path) if row.mapped_counts_file_path else None,
                "min" if "min" in config_df.at[index, "normalization_method"].split(",") else "raw",
                row.normalization_scope,
                required_samples=calling_samples,
            )
            if row.matched_comparison_file_path:
                from lib.matched_statistics import read_matched_comparisons
                availability = {}
                for group, _ in groups:
                    for path in filter(None, group):
                        method, condition, replicate = parse_sample_name(path)
                        availability.setdefault((method, condition), set()).add(replicate)
                for comparison in read_matched_comparisons(Path(row.matched_comparison_file_path)):
                    numerator = availability.get((comparison.assay, comparison.numerator_condition), set())
                    denominator = availability.get((comparison.assay, comparison.denominator_condition), set())
                    if numerator != denominator:
                        raise ValueError(
                            f"Matched comparison {comparison.name!r} requires identical replicate IDs; "
                            f"{comparison.numerator_condition} has {sorted(numerator)}, "
                            f"{comparison.denominator_condition} has {sorted(denominator)}."
                        )
                    if len(numerator) < 2:
                        raise ValueError(
                            f"Matched comparison {comparison.name!r} requires at least two matched replicates."
                        )
            genome_path = Path(row.genome_file_path)
            annotation_path = Path(row.annotation_file_path)
            reference_inputs = fingerprint_inputs(
                [genome_path, annotation_path], cache=fingerprint_cache,
            )
            reference_context = io.load_reference_context(
                genome_path, annotation_path, reference_inputs, cache=reference_cache,
            )
            genome = reference_context.genome_dict()
            expression_files = set()
            if row.alignment_folder_path:
                for group, _ in groups:
                    expression_files.update(io.check_alignment_path_input(
                        Path(row.alignment_folder_path), group[0], group[1]) or [])
            filtered_alignments = [path for group, _ in groups for path in group if path]
            if row.rpkm_read_usage == "specific":
                filtered_alignments.extend(expression_files)
            io.validate_read_length_samples(read_lengths, filtered_alignments)
            io.validate_offset_samples(
                offsets,
                [path for group, _ in groups for path in group if path],
            )
            io.validate_reference_inputs(
                genome,
                annotation_path,
                [path for group, _ in groups for path in group if path] + list(expression_files),
                reference_context=reference_context,
            )
        except (ValueError, OSError, TypeError) as exc:
            raise ValueError(f"Config row {index + 2} ({row.experiment_name}): {exc}") from exc
    config_df.attrs["source_path"] = str(Path(config_sheet).resolve())
    config_df.attrs["source_snapshot"] = source_snapshot
    verify_inputs_unchanged(source_snapshot)
    return config_df
