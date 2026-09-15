"""Exercise legacy helper CLIs and GFF conversion against current table schemas."""

from pathlib import Path
import subprocess
import sys
from urllib.parse import quote

import pandas as pd
import pytest

import helpers.final_to_gff as final_to_gff
from helpers.final_to_gff import create_gff
from lib.merging import write_merged_gff


REPO = Path(__file__).resolve().parents[1]


def result_table():
    return pd.DataFrame([{
        "Identifier": "chr:1-12:+", "Genome": "chr", "Start": 1, "Stop": 12,
        "Strand": "+", "Type": "annotated", "Locus_tag": "gene", "Codon_count": 4,
        "Start_codon": "ATG", "Stop_codon": "TAA", "15nt_window": "AAAAAA",
        "Nucleotide_Seq": "ATGAAAAAATAA", "Amino_Acid_Seq": "MKK*",
        "5'-distance": 0, "3'-distance": 0, "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_q_value": 0.01, "TIS-a-1_peak_correction": "by",
        "TIS-a-1_peak_tests_in_family": 100,
    }])


def test_legacy_merge_cli_accepts_one_file_and_preserves_statistics(tmp_path):
    source = tmp_path / "sample.tsv"
    result_table().to_csv(source, sep="\t", index=False)
    target = tmp_path / "nested" / "merged.xlsx"
    completed = subprocess.run([
        sys.executable, str(REPO / "helpers/merge_tables.py"), "-i", str(source), "-o", str(target),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    merged = pd.read_csv(target.with_suffix(".csv"), sep="\t")
    assert merged.loc[0, "TIS-a-1_peak_q_value"] == 0.01
    assert merged.loc[0, "TIS-a-1_peak_correction"] == "by"
    assert merged.loc[0, "TIS-a-1_peak_tests_in_family"] == 100
    assert target.exists()
    assert target.with_suffix(".gff").exists()


def test_legacy_merge_help_explains_replacement_and_filename_change(tmp_path):
    completed = subprocess.run([sys.executable, str(REPO / "helpers/merge_tables.py"), "--help"],
                               cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0
    assert "Deprecated entry point" in completed.stdout
    assert "filenames no longer create experiment namespaces" in " ".join(completed.stdout.split())


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
def test_final_to_gff_cli_uses_named_columns_and_keeps_zero_log2fc(tmp_path, suffix):
    table = result_table()
    table["Evidence"] = "TIS-a-1,TTS-a-1"
    table["TIS-a-1_RIBO-a-1_log2FC"] = 0.0
    source = tmp_path / f"sample{suffix}"
    table = table[list(reversed(table.columns))]
    if suffix == ".xlsx":
        table.to_excel(source, index=False)
    else:
        table.to_csv(source, sep="\t", index=False)
    target = tmp_path / "nested" / "results.gff"
    completed = subprocess.run([
        sys.executable, str(REPO / "helpers/final_to_gff.py"), "-i", str(source), "-o", str(target),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    fields = target.read_text().splitlines()[1].split("\t")
    assert len(fields) == 9
    assert fields[3:5] == ["1", "12"]
    assert fields[7] == "0"
    assert "color=180,180,180" in fields[8]
    assert "log2fc_tis_a_1_ribo_a_1=0.0000" in fields[8]


def test_final_gff_missing_evidence_uses_peak_assay_then_unknown_gray(tmp_path):
    table = result_table()
    output = tmp_path / "result.gff"
    create_gff(table, output)
    assert "color=0,0,255" in output.read_text()
    table = table.drop(columns=["TIS-a-1_peak_height"])
    table["Evidence"] = pd.NA
    create_gff(table, output)
    assert "color=128,128,128" in output.read_text()


def test_final_gff_empty_table_writes_header(tmp_path):
    target = tmp_path / "empty.gff"
    create_gff(result_table().iloc[:0], target)
    assert target.read_text() == "##gff-version 3\n"


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
def test_final_gff_accepts_merged_legacy_linear_topology_without_lengths(tmp_path, suffix):
    table = result_table()
    table["Reference_length"] = pd.NA
    table["Is_circular"] = False
    source = tmp_path / f"merged-legacy{suffix}"
    if suffix == ".xlsx":
        table.to_excel(source, index=False)
    else:
        table.to_csv(source, sep="\t", index=False)
    target = tmp_path / "merged-legacy.gff"

    final_to_gff.convert_table_to_gff(source, target)

    lines = target.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert len(lines) == 2
    assert lines[1].split("\t")[3:5] == ["1", "12"]


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
@pytest.mark.parametrize(
    "genome,identifier",
    [("001", "0007"), ("NA", "NULL"), ("NULL", "N/A")],
)
def test_final_gff_import_preserves_literal_and_numeric_looking_labels(
        tmp_path, suffix, genome, identifier):
    table = result_table()
    table["Genome"] = genome
    table["Identifier"] = identifier
    # A literal label that pandas normally treats as missing must not fall back
    # to the positive TIS peak in this row.
    table["Evidence"] = "NA"
    source = tmp_path / f"labels{suffix}"
    if suffix == ".xlsx":
        table.to_excel(source, index=False)
    else:
        table.to_csv(source, sep="\t", index=False)
    target = tmp_path / "labels.gff"

    final_to_gff.convert_table_to_gff(source, target)

    fields = target.read_text().splitlines()[1].split("\t")
    assert fields[0] == quote(genome, safe=".:^*$@!+_?-|")
    assert f"ID={identifier}" in fields[8]
    assert "color=128,128,128" in fields[8]


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
def test_final_gff_import_keeps_blank_evidence_missing(tmp_path, suffix):
    table = result_table()
    table["Evidence"] = pd.NA
    source = tmp_path / f"blank{suffix}"
    if suffix == ".xlsx":
        table.to_excel(source, index=False)
    else:
        table.to_csv(source, sep="\t", index=False)
    target = tmp_path / "blank.gff"

    final_to_gff.convert_table_to_gff(source, target)

    assert "color=0,0,255" in target.read_text()


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
def test_final_gff_import_rejects_duplicate_columns(tmp_path, suffix):
    table = result_table()
    table.insert(3, "Genome", "conflicting", allow_duplicates=True)
    source = tmp_path / f"duplicate{suffix}"
    if suffix == ".xlsx":
        table.to_excel(source, index=False)
    else:
        table.to_csv(source, sep="\t", index=False)
    target = tmp_path / "duplicate.gff"

    with pytest.raises(ValueError, match="duplicate columns: Genome"):
        final_to_gff.convert_table_to_gff(source, target)

    assert not target.exists()


def test_final_gff_conversion_failure_does_not_publish_partial_output(
        tmp_path, monkeypatch):
    source = tmp_path / "sample.tsv"
    result_table().to_csv(source, sep="\t", index=False)
    target = tmp_path / "nested" / "result.gff"

    def fail_after_write(frame, output):
        Path(output).write_text("partial")
        raise OSError("GFF write failed")

    monkeypatch.setattr(final_to_gff, "create_gff", fail_after_write)
    with pytest.raises(OSError, match="GFF write failed"):
        final_to_gff.convert_table_to_gff(source, target)

    assert not target.exists()
    assert not list(target.parent.glob(".orfbounder-stage-*"))


def test_final_to_gff_cli_never_overwrites_its_input(tmp_path):
    source = tmp_path / "sample.tsv"
    result_table().to_csv(source, sep="\t", index=False)
    original = source.read_bytes()

    completed = subprocess.run([
        sys.executable, str(REPO / "helpers/final_to_gff.py"),
        "-i", str(source), "-o", str(source),
    ], cwd=tmp_path, capture_output=True, text=True)

    assert completed.returncode == 2
    assert "Output already exists" in completed.stderr
    assert "Traceback" not in completed.stderr
    assert source.read_bytes() == original
    assert not list(tmp_path.glob(".orfbounder-stage-*"))


@pytest.mark.parametrize("field,value", [("Start", 0), ("Start", 1.5), ("Stop", 0), ("Strand", "wrong")])
def test_final_gff_rejects_invalid_coordinates_or_strand(tmp_path, field, value):
    table = result_table().astype({field: "object"})
    table.loc[0, field] = value
    with pytest.raises(ValueError, match="Invalid"):
        create_gff(table, tmp_path / "invalid.gff")


@pytest.mark.parametrize(
    "field,value",
    [("Genome", "other"), ("Start", 2), ("Stop", 11), ("Strand", "-")],
)
def test_final_gff_rejects_coordinate_identifier_disagreement(
        tmp_path, field, value):
    table = result_table().astype({field: "object"})
    table.loc[0, field] = value

    with pytest.raises(ValueError, match="Identifier and coordinate columns disagree"):
        create_gff(table, tmp_path / "contradictory.gff")


def test_final_gff_coordinate_identifier_is_parsed_from_the_right(tmp_path):
    table = result_table()
    table["Genome"] = "plasmid:alpha-beta"
    table["Identifier"] = "plasmid:alpha-beta:1-12:+"
    target = tmp_path / "right-safe.gff"

    create_gff(table, target)

    assert target.read_text().splitlines()[1].startswith(
        "plasmid:alpha-beta\tORFBounder\tCDS\t1\t12\t"
    )


def test_merged_gff_escapes_attributes_and_keeps_evidence_list(tmp_path):
    table = result_table()
    table["Identifier"] = "id;unsafe=1%"
    table["Locus_tag"] = "name&with,delimiters"
    table["Genome"] = "chr space"
    table["Evidence"] = "TIS-a&b-1,TTS-a&b-1"
    # This is a valid Python identifier, so positional namedtuple attributes
    # cannot be relied on to retrieve the value.
    table["comparison_log2FC"] = 0.0
    target = tmp_path / "merged.gff"
    write_merged_gff(table, target)
    fields = target.read_text().splitlines()[1].split("\t")
    assert fields[0] == "chr%20space"
    assert fields[7] == "0"
    assert "ID=id%3Bunsafe%3D1%25" in fields[8]
    assert "Name=name%26with%2Cdelimiters" in fields[8]
    assert "evidence=TIS-a%26b-1,TTS-a%26b-1" in fields[8]
    assert "comparison_log2fc=0.0" in fields[8]
