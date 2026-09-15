#!/usr/bin/env python

"""
Module to merge multiple output tables into one final table
"""

from collections import Counter, OrderedDict
from numbers import Integral, Real
from pathlib import Path
from urllib.parse import quote
import csv

import argparse
import pandas as pd
import numpy as np

from lib import io, run_output
import lib.expression as expr
from lib.matched_statistics import ATTACH_FIELDS
from lib.statistics import STATISTIC_FIELDS
from lib.coordinates import parse_coordinate_id

class OrderedCounter(Counter, OrderedDict):
    pass


_REFERENCE_TOPOLOGY_KEY = object()

_METRIC_INDICES = {"peak_height": 0, "relative_density": 1, "rpkm": 2, "TE": 3}
_EVIDENCE_FIELDS = (*STATISTIC_FIELDS, *ATTACH_FIELDS)
# Prefer the longest suffix. For example, a matched
# ``numerator_peak_count`` must not be mistaken for a local ``peak_count``.
_DYNAMIC_FIELDS = tuple(sorted(
    {*_METRIC_INDICES, *_EVIDENCE_FIELDS}, key=len, reverse=True,
))
_EXACT_NONNEGATIVE_INTEGER_FIELDS = {
    "peak_count", "background_count", "peak_tests_in_family",
    "matched_pair_count", "informative_pair_count", "concordant_pair_count",
    "discordant_pair_count", "tied_pair_count", "numerator_peak_count",
    "numerator_background_count", "denominator_peak_count",
    "denominator_background_count", "tests_in_family",
}
_BOOLEAN_FIELDS = {"peak_significant", "significant"}
_TEXT_COLUMNS = {
    "Type", "Identifier", "Genome", "Strand", "Locus_tag", "Evidence",
    "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq",
    "Amino_Acid_Seq",
}
_MAX_EXACT_BINARY64_INTEGER = 1 << 53


def _split_dynamic_column(column):
    """Return a result column's sample/prefix and known dynamic field."""
    for field in _DYNAMIC_FIELDS:
        suffix = f"_{field}"
        if column.endswith(suffix):
            return column[:-len(suffix)], field
    return None


def _parse_optional_nonnegative_integer(value, column):
    """Parse an exported count without allowing a float round trip."""
    if value is None or value is pd.NA:
        return pd.NA
    if isinstance(value, bool):
        raise ValueError(f"{column} must contain nonnegative decimal integers or blanks.")
    if isinstance(value, Integral):
        if value < 0:
            raise ValueError(f"{column} must contain nonnegative decimal integers or blanks.")
        return int(value)
    if isinstance(value, str):
        if value == "":
            return pd.NA
        if not value.isascii() or not value.isdecimal():
            raise ValueError(f"{column} must contain nonnegative decimal integers or blanks.")
        return int(value)
    if isinstance(value, Real):
        # Spreadsheet numeric cells are binary64 values. Accept only integral
        # values inside its exact-integer range; larger values must have been
        # written as text to be auditable.
        if (not np.isfinite(value) or not float(value).is_integer()
                or value < 0 or value > _MAX_EXACT_BINARY64_INTEGER):
            raise ValueError(f"{column} must contain exact nonnegative integers or blanks.")
        return int(value)
    raise ValueError(f"{column} must contain nonnegative decimal integers or blanks.")


def _parse_optional_text(value):
    """Keep exported identifiers and labels as text, including leading zeros."""
    if value is None or value is pd.NA or pd.isna(value) or value == "":
        return pd.NA
    return str(value)


def _parse_optional_boolean(value, column):
    """Normalize exported support flags without accepting arbitrary scalars."""
    if value is None or value is pd.NA:
        return pd.NA
    missing = pd.isna(value)
    if isinstance(missing, (bool, np.bool_)) and missing:
        return pd.NA
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, Integral) and not isinstance(value, bool):
        if value in (0, 1):
            return bool(value)
    elif isinstance(value, Real):
        if np.isfinite(value) and value in (0.0, 1.0):
            return bool(value)
    elif isinstance(value, str):
        if value == "":
            return pd.NA
        normalized = value.lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    raise ValueError(
        f"{column} must contain booleans (true/false or 1/0) or blanks; "
        f"got {value!r}."
    )


def _result_converters(columns):
    """Preserve text fields and parse raw-count/audit columns exactly."""
    converters = {
        column: _parse_optional_text for column in columns if column in _TEXT_COLUMNS
    }
    for column in columns:
        dynamic = _split_dynamic_column(column)
        if dynamic is None:
            continue
        if dynamic[1] in _EXACT_NONNEGATIVE_INTEGER_FIELDS:
            def parse_integer(value, *, current_column=column):
                return _parse_optional_nonnegative_integer(value, current_column)

            converters[column] = parse_integer
        elif dynamic[1] in _BOOLEAN_FIELDS:
            def parse_boolean(value, *, current_column=column):
                return _parse_optional_boolean(value, current_column)

            converters[column] = parse_boolean
    return converters


def _read_result_table(path):
    """Read one result table without losing literal labels or integer evidence.

    Pandas' default NA vocabulary includes valid annotation labels such as
    ``NA`` and ``NULL``. ORFBounder exports missing cells as blanks, so only
    blanks should become missing on re-import.
    """
    if path.suffix.lower() in {".csv", ".tsv"}:
        columns = pd.read_csv(path, sep="\t", nrows=0).columns
        return pd.read_csv(
            path, sep="\t", converters=_result_converters(columns),
            keep_default_na=False, na_values=[""],
        )
    columns = pd.read_excel(path, sheet_name="CDS", nrows=0).columns
    return pd.read_excel(
        path, sheet_name="CDS", converters=_result_converters(columns),
        keep_default_na=False, na_values=[""],
    )


def _coalesce_measurement(existing, incoming, unique_id, sample, field):
    """Merge one repeated sample field without making input order meaningful."""
    existing_missing = pd.isna(existing)
    incoming_missing = pd.isna(incoming)
    if not isinstance(existing_missing, (bool, np.bool_)) or not isinstance(
            incoming_missing, (bool, np.bool_)):
        raise ValueError(
            f"Merged measurement values must be scalar for ORF {unique_id!r}, "
            f"sample {sample!r}, field {field!r}."
        )
    if existing_missing:
        return existing if incoming_missing else incoming
    if incoming_missing:
        return existing
    try:
        equal = existing == incoming
    except (TypeError, ValueError):
        equal = False
    if not isinstance(equal, (bool, np.bool_)) or not bool(equal):
        raise ValueError(
            f"Conflicting measurements for ORF {unique_id!r}, sample {sample!r}, "
            f"field {field!r}: {existing!r} and {incoming!r}."
        )
    return existing


def extend_combined_dictionary(
    xlsx_df: pd.DataFrame,
    meta_dict: dict,
    dynamic_dict: dict
) -> tuple[dict, dict]:
    """
    collect data from current file and add it to the existing dictionary
    { chrom:start-stop:strand : { wildcard : peak_height, relative_density, RPKM, TE} }
    { chrom:start-stop:strand : metadata }
    """
    reference_lengths, circular_contigs = io.reference_topology_from_dataframe(xlsx_df)
    topology = meta_dict.setdefault(_REFERENCE_TOPOLOGY_KEY, {})
    for chrom, length in reference_lengths.items():
        value = (int(length), chrom in circular_contigs)
        if chrom in topology and topology[chrom] != value:
            raise ValueError(f"Conflicting reference topology for {chrom!r} across result tables.")
        topology[chrom] = value

    column_map = []
    for i, val in enumerate(xlsx_df.columns):
        dynamic = _split_dynamic_column(val)
        if dynamic is not None:
            sample, field = dynamic
            column_map.append((sample, field, i))

    for row in xlsx_df.itertuples(index=False, name=None):
        values = dict(zip(xlsx_df.columns, row))
        gene_type, unique_id = values["Type"], values["Identifier"]
        parsed_chrom, parsed_start, parsed_stop, parsed_strand = parse_coordinate_id(unique_id)
        topology_columns = {"Reference_length", "Is_circular"}.intersection(values)
        if topology_columns:
            chrom = values.get("Genome", parsed_chrom)
            start, stop = int(values.get("Start", parsed_start)), int(values.get("Stop", parsed_stop))
            strand = values.get("Strand", parsed_strand)
            if (str(chrom), start, stop, strand) != (parsed_chrom, parsed_start, parsed_stop, parsed_strand):
                raise ValueError(f"Identifier and coordinate columns disagree for {unique_id!r}.")
        else:
            # Historical merge inputs treated Identifier as authoritative.
            chrom, start, stop, strand = parsed_chrom, parsed_start, parsed_stop, parsed_strand
        gene_name, codon_count = values["Locus_tag"], values["Codon_count"]
        start_codon, stop_codon = values["Start_codon"], values["Stop_codon"]
        nt_upstream, nt_seq, aa_seq = values["15nt_window"], values["Nucleotide_Seq"], values["Amino_Acid_Seq"]
        fiveprime, threeprime = values["5'-distance"], values["3'-distance"]

        if topology_columns and topology_columns != {"Reference_length", "Is_circular"}:
            raise ValueError("Result tables must provide both Reference_length and Is_circular.")
        if topology_columns:
            reference_length = values["Reference_length"]
            is_circular = values["Is_circular"]
            if (isinstance(reference_length, bool) or pd.isna(reference_length)
                    or int(reference_length) != float(reference_length) or int(reference_length) < 1):
                raise ValueError(f"Invalid Reference_length for {unique_id!r}.")
            reference_length = int(reference_length)
            if not isinstance(is_circular, (bool, np.bool_)):
                raise ValueError(f"Is_circular must be boolean for {unique_id!r}.")
            is_circular = bool(is_circular)
            if not 1 <= start <= reference_length or stop < start:
                raise ValueError(f"Invalid coordinates for {unique_id!r}.")
            if (is_circular and stop - start + 1 > reference_length) or (
                    not is_circular and stop > reference_length):
                raise ValueError(f"Coordinates exceed the reference for {unique_id!r}.")
        else:
            reference_length, is_circular = pd.NA, False
        metadata = (
            gene_type, gene_name, codon_count, start_codon, stop_codon,
            nt_upstream, nt_seq, aa_seq, fiveprime, threeprime,
            chrom, start, stop, strand, reference_length, is_circular,
        )
        previous = meta_dict.get(unique_id)
        if previous is not None:
            # The same coordinate may carry different measurements, never a
            # different reference topology or biological identity.
            def same(left, right):
                if pd.isna(left) and pd.isna(right):
                    return True
                return bool(left == right) if not (pd.isna(left) or pd.isna(right)) else False
            if (len(previous) != len(metadata)
                    or any(not same(left, right) for left, right in zip(previous, metadata))):
                raise ValueError(f"Conflicting metadata for ORF {unique_id!r} across result tables.")
        meta_dict[unique_id] = metadata
        if unique_id not in dynamic_dict:
            dynamic_dict[unique_id] = {}

        for sample, field, column in column_map:
            record = dynamic_dict[unique_id].setdefault(sample, [np.nan, np.nan, np.nan, np.nan, {}])
            if len(record) < 5:
                record.append({})
            elif not isinstance(record[4], dict):
                record[4] = {}
            if field in _METRIC_INDICES:
                index = _METRIC_INDICES[field]
                record[index] = _coalesce_measurement(
                    record[index], row[column], unique_id, sample, field,
                )
                if field == "peak_height":
                    record[4]["__peak_height__"] = True
            else:
                if field in record[4]:
                    record[4][field] = _coalesce_measurement(
                        record[4][field], row[column], unique_id, sample, field,
                    )
                else:
                    record[4][field] = row[column]

    return meta_dict, dynamic_dict


def get_log_fc_contrast(wildcards: list[str]) -> list[tuple[str, str]]:
    """
    Create log fold change contrast list
    """
    contrasts = []

    wildcard_dict = {}
    for card in wildcards:
        if "-" not in card:
            continue

        parts = card.split("-", 2)
        if len(parts) != 3:
            continue
        method, condition, replicate = parts
        if "ribo" in method.lower():
            wildcard_dict[f"{condition}-{replicate}"] = []

    for card in wildcards:
        if "-" not in card:
            continue

        parts = card.split("-", 2)
        if len(parts) != 3:
            continue
        method, condition, replicate = parts
        if method.lower() in ["tts", "tis"] and f"{condition}-{replicate}" in wildcard_dict:
            wildcard_dict[f"{condition}-{replicate}"].append(card)

    for ribo, tt_list in wildcard_dict.items():
        for val in tt_list:
            contrasts.append((f"RIBO-{ribo}", val))

    return contrasts

def calculate_fold_changes(row: pd.Series, tts: str, ribo: str, min_val: float) -> float:
    """
    Calculate the contrast between two columns
    """
    tts_heights = np.float64(row[f"{tts}_peak_height"])
    ribo_heights = np.float64(row[f"{ribo}_peak_height"])

    if np.isnan(tts_heights) or tts_heights <= 0:
        return np.nan
    if (np.isnan(ribo_heights) or ribo_heights == 0) and tts_heights > 0:
        return np.log2(tts_heights / min_val)

    return np.log2(tts_heights / ribo_heights)


def _calculate_fold_changes_vectorized(
    tts_heights: np.ndarray,
    ribo_heights: np.ndarray,
    min_val: float,
) -> np.ndarray:
    """Calculate one contrast without constructing a Series for every ORF."""
    tts_values = np.asarray(tts_heights, dtype=np.float64)
    ribo_values = np.asarray(ribo_heights, dtype=np.float64)
    if tts_values.shape != ribo_values.shape:
        raise ValueError("Contrasted peak-height arrays must have the same shape.")

    result = np.full(tts_values.shape, np.nan, dtype=np.float64)
    valid_tts = ~np.isnan(tts_values) & (tts_values > 0)
    use_minimum = valid_tts & (np.isnan(ribo_values) | (ribo_values == 0))
    use_observed = valid_tts & ~use_minimum
    with np.errstate(divide="ignore", invalid="ignore"):
        result[use_minimum] = np.log2(tts_values[use_minimum] / min_val)
        result[use_observed] = np.log2(
            tts_values[use_observed] / ribo_values[use_observed]
        )
    return result

def build_merged_dataframe(
    meta_dict: dict,
    dynamic_dict: dict
) -> tuple[pd.DataFrame, list[str]]:
    """Merge sample measurements without changing their statistical families."""
    wildcards = sorted({card for records in dynamic_dict.values() for card in records})

    def extras(record):
        return record[4] if len(record) > 4 and isinstance(record[4], dict) else {}

    peak_cards = [card for card in wildcards if any(
        extras(records.get(card, [])).get("__peak_height__", False)
        for records in dynamic_dict.values()
    ) or card.split("-")[0] in {"TIS", "TTS", "RIBO"}]
    te_cards = expr.get_te_header(wildcards)
    statistic_columns = [(card, field) for card in wildcards for field in _EVIDENCE_FIELDS if any(
        field in extras(records.get(card, [])) for records in dynamic_dict.values()
    )]
    identity_header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count"]
    peak_header = [f"{card}_peak_height" for card in peak_cards]
    metadata_header = ["Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", "5'-distance", "3'-distance"]
    header = identity_header + peak_header + ["Evidence"] + metadata_header \
        + [f"{card}_relative_density" for card in peak_cards] \
        + [f"{card}_rpkm" for card in wildcards] \
        + [f"{card}_TE" for card in te_cards] \
        + [f"{card}_{field}" for card, field in statistic_columns] \
        + ["Reference_length", "Is_circular"]

    rows = []
    for unique_id, metadata in meta_dict.items():
        if unique_id is _REFERENCE_TOPOLOGY_KEY:
            continue
        chrom, start, stop, strand = metadata[10:14]
        sample_records = dynamic_dict[unique_id]

        def measurement(card, index):
            record = sample_records.get(card)
            return record[index] if record is not None else np.nan

        row = dict(zip(identity_header, [metadata[0], unique_id, chrom, start, stop, strand, metadata[1], metadata[2]]))
        row.update(zip(metadata_header, metadata[3:]))
        row["Reference_length"], row["Is_circular"] = metadata[14:16]
        row["Evidence"] = ",".join(card for card in peak_cards
            if card.split("-")[0] != "RIBO" and pd.notna(measurement(card, 0)) and measurement(card, 0) > 0)
        for card in peak_cards:
            row[f"{card}_peak_height"] = measurement(card, 0)
            row[f"{card}_relative_density"] = measurement(card, 1)
        for card in wildcards:
            row[f"{card}_rpkm"] = measurement(card, 2)
        for card in te_cards:
            row[f"{card}_TE"] = measurement(card, 3)
        for card, field in statistic_columns:
            row[f"{card}_{field}"] = extras(sample_records.get(card, [])).get(field, pd.NA)
        rows.append(row)

    result_df = pd.DataFrame.from_records(rows, columns=header)
    topology = meta_dict.get(_REFERENCE_TOPOLOGY_KEY, {})
    result_df.attrs["reference_lengths"] = {
        chrom: values[0] for chrom, values in topology.items()
    }
    result_df.attrs["circular_contigs"] = sorted(
        chrom for chrom, values in topology.items() if values[1]
    )
    contrasts = [(ribo, tis) for ribo, tis in get_log_fc_contrast(wildcards)
                 if ribo in peak_cards and tis in peak_cards]
    if contrasts and not result_df.empty:
        ts_columns = [f"{card}_peak_height" for card in peak_cards if card.split("-")[0] in {"TIS", "TTS"}]
        result_df = result_df.loc[(result_df[ts_columns] > 0).any(axis="columns")].copy()
        heights = result_df[peak_header]
        min_val = heights.where(heights > 0).min().min() * 0.9
        contrast_header = []
        for ribo, tis in contrasts:
            name = f"{tis}_{ribo}_log2FC"
            result_df[name] = _calculate_fold_changes_vectorized(
                result_df[f"{tis}_peak_height"].to_numpy(
                    dtype=np.float64, na_value=np.nan,
                ),
                result_df[f"{ribo}_peak_height"].to_numpy(
                    dtype=np.float64, na_value=np.nan,
                ),
                min_val,
            )
            contrast_header.append(name)
        insertion = len(identity_header) + len(peak_header)
        result_df = result_df[header[:insertion] + contrast_header + header[insertion:]]

    return result_df, wildcards


def screen_input_tables(table_list: list[Path]) -> tuple[dict, dict]:
    """
    screen over the input tables
    """
    meta_dict, dynamic_dict = {}, {}
    for table in table_list:
        table = Path(table)
        xlsx_df = _read_result_table(table)
        meta_dict, dynamic_dict = extend_combined_dictionary(xlsx_df, meta_dict, dynamic_dict)

    return meta_dict, dynamic_dict


def write_merged_table(
    meta_dict: dict,
    dynamic_dict: dict,
    output_path: Path,
    matched_results=(),
) -> pd.DataFrame:
    """
    create final merged table and write it to xlsx/csv file
    """

    df_res, _ = build_merged_dataframe(meta_dict, dynamic_dict)
    if matched_results:
        from lib.matched_statistics import attach_matched_to_orfs
        for prefix, assay, statistics in matched_results:
            df_res = attach_matched_to_orfs(df_res, statistics, assay, prefix)

    required_columns = {"Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count",
                        "Evidence", "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq",
                        "5'-distance", "3'-distance", "Reference_length", "Is_circular"}
    empty_optional = [column for column in df_res if column not in required_columns and df_res[column].isna().all()]
    df_res.drop(columns=empty_optional, inplace=True)
    df_res = df_res.sort_values(by=["Genome", "Start", "Stop", "Strand"])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    csv_path = output_path.with_suffix(".csv")
    # GFF3 percent-decoding can legitimately place tabs, newlines, or quotes in
    # display labels.  Quote only fields that need it so the merged TSV remains
    # both ordinary-looking and round-trippable.
    df_res.to_csv(csv_path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)

    io.excel_writer(output_path, {"CDS": df_res})
    return df_res


def write_merged_gff(res_df: pd.DataFrame, output_path: Path) -> None:
    """
    Write complete ORF CDS features with phase 0 and escaped GFF3 attributes.

    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """
    def escape(value):
        return "".join(f"%{ord(char):02X}" if char in "%;&=," or ord(char) < 32 or ord(char) == 127 else char
                       for char in str(value))

    reference_lengths, circular_contigs = io.reference_topology_from_dataframe(res_df)
    lines = ["##gff-version 3\n", *io.gff3_topology_lines(reference_lengths, circular_contigs)]
    for row in res_df.to_dict(orient="records"):
        attributes = [f"ID={escape(row['Identifier'])}"]
        for name, column in (("Name", "Locus_tag"), ("type", "Type")):
            if pd.notna(row[column]) and str(row[column]):
                attributes.append(f"{name}={escape(row[column])}")
        evidence = row.get("Evidence")
        if pd.notna(evidence) and str(evidence):
            attributes.append("evidence=" + ",".join(escape(sample) for sample in str(evidence).split(",")))
        for column, value in row.items():
            if column.endswith("_log2FC") and pd.notna(value):
                attributes.append(f"{escape(column.lower())}={escape(value)}")
        chrom = str(row["Genome"])
        segments = io.gff3_feature_segments(
            int(row["Start"]), int(row["Stop"]), row["Strand"],
            reference_lengths.get(chrom), chrom in circular_contigs,
        )
        for segment_start, segment_stop, phase in segments:
            fields = [quote(chrom, safe=".:^*$@!+_?-|"), "ORFBounder", "CDS",
                      str(segment_start), str(segment_stop), ".", row["Strand"], phase,
                      ";".join(attributes)]
            lines.append("\t".join(fields) + "\n")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gff_path = output_path.with_suffix(".gff")
    gff_path.write_text("".join(lines), encoding="utf-8")


def merge_tables(table_list: list[Path], output_path: Path) -> None:
    """
    collect information from all input tables and merge them into one final output table
    """
    output_path = Path(output_path)
    if output_path.suffix.lower() != ".xlsx":
        raise ValueError("Merged output path must use the .xlsx extension.")
    meta_dict, dynamic_dict = screen_input_tables(table_list)
    with run_output.staged_output(output_path.parent) as staging:
        staged_path = staging / output_path.name
        res_df = write_merged_table(meta_dict, dynamic_dict, staged_path)
        write_merged_gff(res_df, staged_path)


def main() -> None:
    """
    Command line interface for merging multiple result tables into one final table
    """
    parser = argparse.ArgumentParser(description='Merge all results for different files into one xlsx and csv file.')
    parser.add_argument("-t", "--tables", nargs="+", dest="table_list", type=Path,
                        required=True, help="list of input tables.")
    parser.add_argument("-o", "--output_path", action="store", dest="output_path",
                        type=Path, required=True, help="Output file .xlsx format")
    args = parser.parse_args()

    try:
        merge_tables(args.table_list, args.output_path)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"ORFBounder merge error: {exc}\n")


if __name__ == '__main__':
    main()
