#!/usr/bin/env python3
"""Export named ORF result columns to GFF3, colored by TIS/TTS evidence.

Unknown or missing assay evidence is gray. Log2 fold changes, including zero,
are retained under lowercase attributes derived from their full column names.
Escaping and CDS phase follow the Sequence Ontology GFF3 specification:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

import argparse
import csv
from pathlib import Path
import re
import sys
from urllib.parse import quote

import pandas as pd

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lib import io, run_output
from lib.coordinates import parse_coordinate_id


COLOR_DICT = {"TIS_only": "0,0,255", "TTS_only": "255,0,255", "both": "180,180,180", "unknown": "128,128,128"}
REQUIRED_COLUMNS = {"Identifier", "Genome", "Start", "Stop", "Strand"}
TEXT_COLUMNS = {"Identifier", "Genome", "Strand", "Evidence"}


def _escape_attribute(value):
    return "".join(f"%{ord(char):02X}" if char in "%;=,&" or ord(char) < 32 or ord(char) == 127 else char
                   for char in str(value))


def _assays(row):
    evidence = row.get("Evidence")
    samples = [] if pd.isna(evidence) or not str(evidence).strip() else str(evidence).split(",")
    if not samples:
        samples = [column[:-12] for column, value in row.items()
                   if column.endswith("_peak_height") and pd.notna(value) and float(value) > 0]
    return {re.split("[-_]", sample.strip(), maxsplit=1)[0].upper() for sample in samples}


def _topology(df):
    return io.reference_topology_from_dataframe(df)


def _feature_segments(start, stop, strand, reference_length=None, circular=False):
    return io.gff3_feature_segments(start, stop, strand, reference_length, circular)


def _optional_text(value):
    """Preserve result labels verbatim while representing blank cells as missing."""
    if value is None or value is pd.NA or pd.isna(value) or value == "":
        return pd.NA
    return str(value)


def _input_columns(input_path):
    """Read unmangled column names so ambiguous tables can be rejected."""
    if input_path.suffix.lower() in {".csv", ".tsv"}:
        with input_path.open("r", encoding="utf-8-sig", newline="") as stream:
            columns = next(csv.reader(stream, delimiter="\t"), None)
    else:
        header = pd.read_excel(
            input_path, sheet_name=0, header=None, nrows=1,
            keep_default_na=False,
        )
        columns = None if header.empty else header.iloc[0].tolist()
    if not columns:
        raise ValueError("Input table has no header row.")
    return [str(column) for column in columns]


def _read_result_table(input_path):
    """Read a result table without losing labels or accepting duplicate columns."""
    columns = _input_columns(input_path)
    duplicates = sorted({column for column in columns if columns.count(column) > 1})
    if duplicates:
        raise ValueError(f"Input table has duplicate columns: {', '.join(duplicates)}")
    converters = {
        column: _optional_text for column in columns if column in TEXT_COLUMNS
    }
    if input_path.suffix.lower() in {".csv", ".tsv"}:
        return pd.read_csv(
            input_path, sep="\t", converters=converters,
            keep_default_na=False, na_values=[""],
        )
    return pd.read_excel(
        input_path, sheet_name=0, converters=converters,
        keep_default_na=False, na_values=[""],
    )


def create_gff(df, output):
    """Write 1-based inclusive ORF coordinates using labels, including empty tables."""
    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(f"Input table is missing required columns: {', '.join(sorted(missing))}")
    reference_lengths, circular_contigs = _topology(df)
    lines = ["##gff-version 3\n"]
    for chrom in sorted(circular_contigs):
        genome = quote(chrom, safe=".:^*$@!+_?-|")
        length = reference_lengths[chrom]
        identifier = _escape_attribute(f"{chrom}:region")
        lines.extend([
            f"##sequence-region {genome} 1 {length}\n",
            f"{genome}\tORFBounder\tregion\t1\t{length}\t.\t.\t.\t"
            f"ID={identifier};Is_circular=true\n",
        ])
    for row in df.to_dict(orient="records"):
        if any(pd.isna(row[column]) for column in REQUIRED_COLUMNS):
            raise ValueError("ORF identifiers, coordinates and strands must not be missing.")
        start, stop = int(row["Start"]), int(row["Stop"])
        if start < 1 or start > stop or start != float(row["Start"]) or stop != float(row["Stop"]):
            raise ValueError(f"Invalid 1-based inclusive coordinates for {row['Identifier']}.")
        if reference_lengths:
            length = reference_lengths[str(row["Genome"])]
            circular = str(row["Genome"]) in circular_contigs
            if start > length or ((not circular and stop > length)
                                  or (circular and stop - start + 1 > length)):
                raise ValueError(f"Invalid coordinates for {row['Identifier']}.")
        if row["Strand"] not in {"+", "-"}:
            raise ValueError(f"Invalid ORF strand for {row['Identifier']}: {row['Strand']}")
        try:
            identifier_coordinates = parse_coordinate_id(row["Identifier"])
        except ValueError:
            # Historical result tables can use opaque feature identifiers.
            # When an identifier does encode coordinates, however, it must not
            # silently contradict the columns used to create the GFF feature.
            pass
        else:
            table_coordinates = (str(row["Genome"]), start, stop, row["Strand"])
            if identifier_coordinates != table_coordinates:
                raise ValueError(
                    f"Identifier and coordinate columns disagree for "
                    f"{row['Identifier']!r}."
                )
        assays = _assays(row)
        category = "both" if {"TIS", "TTS"}.issubset(assays) else (
            "TIS_only" if "TIS" in assays else "TTS_only" if "TTS" in assays else "unknown"
        )
        identifier = _escape_attribute(row["Identifier"])
        attributes = [f"ID={identifier}", f"Name={identifier}", f"color={COLOR_DICT[category]}"]
        for column, value in row.items():
            if column.endswith("_log2FC") and pd.notna(value):
                name = "log2fc_" + re.sub(r"[^a-z0-9_]", "_", column[:-7].lower())
                attributes.append(f"{name}={float(value):.4f}")
        genome = quote(str(row["Genome"]), safe=".:^*$@!+_?-|")
        length = reference_lengths.get(str(row["Genome"]))
        circular = str(row["Genome"]) in circular_contigs
        for segment_start, segment_stop, phase in _feature_segments(
                start, stop, row["Strand"], length, circular):
            fields = [genome, "ORFBounder", "CDS", str(segment_start), str(segment_stop),
                      ".", row["Strand"], phase, ";".join(attributes)]
            lines.append("\t".join(fields) + "\n")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(lines), encoding="utf-8")


def convert_table_to_gff(input_path, output_path):
    """Read one result table and exclusively publish its complete GFF3 file."""
    input_path, output_path = Path(input_path), Path(output_path)
    frame = _read_result_table(input_path)
    with run_output.staged_output(output_path.parent) as staging:
        create_gff(frame, staging / output_path.name)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Create GFF3 colored by TIS/TTS evidence; unknown evidence is gray.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Excel or tab-separated CSV/TSV result table.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output GFF3 file.")
    args = parser.parse_args(argv)
    try:
        convert_table_to_gff(args.input, args.output)
    except (OSError, ValueError) as exc:
        parser.error(f"Cannot create GFF: {exc}")


if __name__ == "__main__":
    main()
