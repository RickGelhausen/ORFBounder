"""Regression tests for merging sparse measurements and statistical evidence."""

from pathlib import Path
import csv
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from lib import merging


REPOSITORY = Path(__file__).resolve().parents[1]


def sample_table(**measurements):
    return pd.DataFrame([{
        "Type": "annotated", "Identifier": "chr:1-12:+", "Genome": "chr",
        "Start": 1, "Stop": 12, "Strand": "+", "Locus_tag": "gene",
        "Codon_count": 4, "Start_codon": "ATG", "Stop_codon": "TAA",
        "15nt_window": "AAAAA", "Nucleotide_Seq": "ATGAAAAAATAA",
        "Amino_Acid_Seq": "MKK*", "5'-distance": 0, "3'-distance": 0,
        **measurements,
    }])


def test_rna_only_measurement_is_stored_in_rpkm_slot():
    table = sample_table(**{"TIS-a-1_peak_height": 20, "RNA-a-1_rpkm": 125})
    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})
    assert dynamic["chr:1-12:+"]["RNA-a-1"][2] == 125
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.loc[0, "RNA-a-1_rpkm"] == 125
    assert "RNA-a-1_peak_height" not in result.columns


def test_te_only_measurement_is_stored_in_te_slot():
    table = sample_table(**{"TIS-a-1_TE": 3.5, "RNATIS-a-1_rpkm": 10})
    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})
    assert dynamic["chr:1-12:+"]["TIS-a-1"][3] == 3.5
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.loc[0, "TIS-a-1_TE"] == 3.5


def test_reordered_columns_and_plain_sample_names_keep_peak_and_statistics():
    table = sample_table(**{
        "sample_peak_height": 20, "sample_relative_density": 4,
        "sample_peak_count": 26, "sample_background_count": 10,
        "sample_peak_fold_enrichment": 5.2,
        "sample_peak_p_value": 0.001, "sample_peak_q_value": 0.05,
        "sample_peak_significant": True, "sample_peak_correction": "by",
        "sample_peak_tests_in_family": 100,
    })
    metadata, dynamic = merging.extend_combined_dictionary(table[sorted(table.columns)], {}, {})
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.loc[0, "sample_peak_height"] == 20
    assert result.loc[0, "sample_relative_density"] == 4
    assert result.loc[0, "sample_peak_fold_enrichment"] == 5.2
    assert result.loc[0, "sample_peak_q_value"] == 0.05
    assert result.loc[0, "sample_peak_tests_in_family"] == 100
    assert result.loc[0, "sample_peak_correction"] == "by"
    assert result.loc[0, "Evidence"] == "sample"


def test_missing_statistics_in_another_sample_stay_missing():
    first = sample_table(**{"TIS-a-1_peak_height": 20, "TIS-a-1_peak_q_value": 0.02})
    second = sample_table(**{"TIS-b-1_peak_height": 30})
    second["Identifier"] = "chr:20-31:+"
    metadata, dynamic = merging.extend_combined_dictionary(first, {}, {})
    metadata, dynamic = merging.extend_combined_dictionary(second, metadata, dynamic)
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.loc[0, "TIS-a-1_peak_q_value"] == 0.02
    assert pd.isna(result.loc[1, "TIS-a-1_peak_q_value"])


@pytest.mark.parametrize(
    "field,left,right",
    [
        ("TIS-a-1_peak_height", 20, 30),
        ("TIS-a-1_peak_q_value", 0.02, 0.03),
        ("TIS-a-1_peak_correction", "by", "bh"),
    ],
)
def test_repeated_sample_fields_reject_conflicting_nonmissing_values(
        field, left, right):
    first_values = {"TIS-a-1_peak_height": 20, field: left}
    second_values = {"TIS-a-1_peak_height": 20, field: right}
    metadata, dynamic = merging.extend_combined_dictionary(
        sample_table(**first_values), {}, {},
    )

    with pytest.raises(
            ValueError,
            match=r"ORF 'chr:1-12:\+', sample 'TIS-a-1', field",
    ):
        merging.extend_combined_dictionary(
            sample_table(**second_values), metadata, dynamic,
        )


@pytest.mark.parametrize("missing_first", [False, True])
def test_repeated_sample_fields_coalesce_missing_values_order_independently(
        missing_first):
    complete = {
        "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_q_value": 0.02,
    }
    missing = {
        "TIS-a-1_peak_height": np.nan,
        "TIS-a-1_peak_q_value": pd.NA,
    }
    first, second = (missing, complete) if missing_first else (complete, missing)
    metadata, dynamic = merging.extend_combined_dictionary(
        sample_table(**first), {}, {},
    )
    metadata, dynamic = merging.extend_combined_dictionary(
        sample_table(**second), metadata, dynamic,
    )

    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert result.loc[0, "TIS-a-1_peak_height"] == 20
    assert result.loc[0, "TIS-a-1_peak_q_value"] == 0.02
    assert result.loc[0, "Evidence"] == "TIS-a-1"


def test_repeated_equal_sample_values_are_idempotent_and_evidence_is_unique():
    table = sample_table(**{
        "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_significant": True,
    })
    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})
    metadata, dynamic = merging.extend_combined_dictionary(table, metadata, dynamic)

    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert result.loc[0, "TIS-a-1_peak_height"] == 20
    assert bool(result.loc[0, "TIS-a-1_peak_significant"])
    assert result.loc[0, "Evidence"] == "TIS-a-1"


@pytest.mark.parametrize("reverse_order", [False, True])
@pytest.mark.parametrize(
    "field,conflicting_value",
    [
        ("15nt_window", "CCCCC"),
        ("Amino_Acid_Seq", "MQQ*"),
        ("5'-distance", 7),
        ("3'-distance", 9),
    ],
)
def test_shared_orf_metadata_conflicts_are_rejected_in_either_input_order(
        field, conflicting_value, reverse_order):
    original = sample_table(**{"TIS-a-1_peak_height": 20})
    conflicting = original.copy()
    conflicting.loc[0, field] = conflicting_value
    first, second = (
        (conflicting, original) if reverse_order else (original, conflicting)
    )
    metadata, dynamic = merging.extend_combined_dictionary(first, {}, {})

    with pytest.raises(ValueError, match="Conflicting metadata for ORF"):
        merging.extend_combined_dictionary(second, metadata, dynamic)


def test_statistical_columns_do_not_count_as_height_evidence():
    table = sample_table(**{
        "TIS-a-1_peak_height": np.nan, "RIBO-a-1_peak_height": 12,
        "TIS-a-1_peak_q_value": 0.9, "TIS-a-1_peak_count": 3,
    })
    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.empty


def test_zero_height_does_not_create_invalid_fold_change():
    row = pd.Series({"TIS-a-1_peak_height": 0, "RIBO-a-1_peak_height": 0})
    assert np.isnan(merging.calculate_fold_changes(row, "TIS-a-1", "RIBO-a-1", 1))


def test_vectorized_fold_changes_match_scalar_randomized_edge_cases():
    rng = np.random.default_rng(20260911)
    tis_heights = rng.normal(loc=10, scale=15, size=10_000)
    ribo_heights = rng.normal(loc=8, scale=12, size=10_000)
    edge_values = np.array([np.nan, 0.0, -0.0, -1.0, np.inf, -np.inf])
    tis_heights[:len(edge_values)] = edge_values
    ribo_heights[:len(edge_values)] = edge_values[::-1]
    tis_heights[100:106] = 20.0
    ribo_heights[100:106] = edge_values
    minimum = 0.45

    with np.errstate(divide="ignore", invalid="ignore"):
        expected = np.array([
            merging.calculate_fold_changes(
                pd.Series({
                    "TIS-a-1_peak_height": tis,
                    "RIBO-a-1_peak_height": ribo,
                }),
                "TIS-a-1", "RIBO-a-1", minimum,
            )
            for tis, ribo in zip(tis_heights, ribo_heights)
        ])
    observed = merging._calculate_fold_changes_vectorized(
        tis_heights, ribo_heights, minimum,
    )

    np.testing.assert_array_equal(observed, expected)


def test_vectorized_fold_changes_reject_misaligned_arrays():
    with pytest.raises(ValueError, match="same shape"):
        merging._calculate_fold_changes_vectorized(
            np.array([1.0]), np.array([1.0, 2.0]), 0.5,
        )


def test_merged_contrast_values_and_column_position_match_scalar_contract():
    tables = []
    heights = [(20.0, 10.0), (10.0, 0.0), (5.0, np.nan), (0.0, 10.0)]
    for index, (tis_height, ribo_height) in enumerate(heights):
        table = sample_table(**{
            "TIS-a-1_peak_height": tis_height,
            "RIBO-a-1_peak_height": ribo_height,
        })
        table.loc[0, "Identifier"] = f"chr:{index * 20 + 1}-{index * 20 + 12}:+"
        tables.append(table)
    metadata, dynamic = merging.extend_combined_dictionary(
        pd.concat(tables, ignore_index=True), {}, {},
    )

    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    contrast = "TIS-a-1_RIBO-a-1_log2FC"
    assert result.columns[8:12].tolist() == [
        "RIBO-a-1_peak_height", "TIS-a-1_peak_height", contrast, "Evidence",
    ]
    assert len(result) == 3
    minimum = result[[
        "RIBO-a-1_peak_height", "TIS-a-1_peak_height",
    ]].where(lambda values: values > 0).min().min() * 0.9
    expected = [
        merging.calculate_fold_changes(row, "TIS-a-1", "RIBO-a-1", minimum)
        for _, row in result.iterrows()
    ]
    np.testing.assert_array_equal(result[contrast].to_numpy(), expected)


def test_empty_merge_writes_valid_tables_and_gff(tmp_path):
    path = tmp_path / "empty.xlsx"
    result = merging.write_merged_table({}, {}, path)
    merging.write_merged_gff(result, path)
    assert result.empty
    assert "Evidence" in result.columns
    assert path.exists()
    assert "Identifier" in pd.read_csv(path.with_suffix(".csv"), sep="\t").columns
    assert path.with_suffix(".gff").read_text() == "##gff-version 3\n"


def test_tab_separated_result_tables_are_mergeable(tmp_path):
    path = tmp_path / "sample.csv"
    sample_table(**{"sample_peak_height": 20}).to_csv(path, sep="\t", index=False)
    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert result.loc[0, "sample_peak_height"] == 20


@pytest.mark.parametrize("suffix", [".csv", ".xlsx"])
@pytest.mark.parametrize("label", ["NA", "NULL", "N/A"])
def test_merge_round_trips_literal_na_annotation_labels(tmp_path, suffix, label):
    table = sample_table(**{"TIS-a-1_peak_height": 20})
    table.loc[0, "Locus_tag"] = label
    path = tmp_path / f"sample{suffix}"
    if suffix == ".xlsx":
        merging.io.excel_writer(path, {"CDS": table})
    else:
        table.to_csv(path, sep="\t", index=False)

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert result.at[0, "Locus_tag"] == label


@pytest.mark.parametrize("suffix", [".csv", ".xlsx"])
def test_merge_round_trips_numeric_looking_reference_and_locus_labels(
        tmp_path, suffix):
    table = sample_table(**{
        "TIS-a-1_peak_height": 20,
        "Reference_length": 30,
        "Is_circular": False,
    })
    table.loc[0, "Identifier"] = "001:1-12:+"
    table.loc[0, "Genome"] = "001"
    table.loc[0, "Locus_tag"] = "0007"
    path = tmp_path / f"sample{suffix}"
    if suffix == ".xlsx":
        merging.io.excel_writer(path, {"CDS": table})
    else:
        table.to_csv(path, sep="\t", index=False)

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert result.at[0, "Identifier"] == "001:1-12:+"
    assert result.at[0, "Genome"] == "001"
    assert result.at[0, "Locus_tag"] == "0007"
    gff_path = tmp_path / "merged.xlsx"
    merging.write_merged_gff(result, gff_path)
    gff = gff_path.with_suffix(".gff").read_text()
    assert "001\tORFBounder\tCDS\t1\t12\t.\t+\t0\tID=001:1-12:+;Name=0007" in gff


@pytest.mark.parametrize("suffix", [".csv", ".xlsx"])
def test_merge_keeps_blank_labels_and_all_missing_measurements_missing(
        tmp_path, suffix):
    table = sample_table(**{"TIS-a-1_peak_height": np.nan})
    table.loc[0, "Locus_tag"] = ""
    path = tmp_path / f"sample{suffix}"
    if suffix == ".xlsx":
        merging.io.excel_writer(path, {"CDS": table})
    else:
        # Force a quoted empty field to establish that only literal NA tokens
        # survive re-import, not the exporter representation of missing data.
        table.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_ALL)

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert pd.isna(result.at[0, "Locus_tag"])
    assert pd.isna(result.at[0, "TIS-a-1_peak_height"])


@pytest.mark.parametrize("suffix", [".csv", ".tsv"])
def test_tabular_merge_preserves_nullable_counts_above_float_precision(
        tmp_path, suffix):
    exact_count = (1 << 53) + 1
    first = sample_table(**{
        "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_count": exact_count,
    })
    first.loc[0, "Locus_tag"] = 'gene\tline\n"quoted"'
    second = sample_table(**{
        "TIS-a-1_peak_height": 30,
        "TIS-a-1_peak_count": pd.NA,
    })
    second.loc[0, ["Identifier", "Start", "Stop"]] = ["chr:21-32:+", 21, 32]
    table = pd.concat([first, second], ignore_index=True)
    table["TIS-a-1_peak_count"] = pd.array(
        [exact_count, pd.NA], dtype=object,
    )
    path = tmp_path / f"sample{suffix}"
    table.to_csv(path, sep="\t", index=False)

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    result = result.set_index("Identifier")

    assert result.at["chr:1-12:+", "TIS-a-1_peak_count"] == exact_count
    assert isinstance(result.at["chr:1-12:+", "TIS-a-1_peak_count"], int)
    assert pd.isna(result.at["chr:21-32:+", "TIS-a-1_peak_count"])
    assert result.at["chr:1-12:+", "Locus_tag"] == 'gene\tline\n"quoted"'


def test_excel_merge_preserves_large_counts_written_as_text(tmp_path):
    exact_count = (1 << 63) + 1
    table = sample_table(**{
        "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_count": exact_count,
    })
    table["TIS-a-1_peak_count"] = pd.array([exact_count], dtype=object)
    path = tmp_path / "sample.xlsx"
    merging.io.excel_writer(path, {"CDS": table})

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    assert result.at[0, "TIS-a-1_peak_count"] == exact_count


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
@pytest.mark.parametrize(
    "column",
    ["TIS-a-1_peak_significant", "treated_vs_control_TIS_significant"],
)
def test_merge_normalizes_nullable_statistical_booleans(tmp_path, suffix, column):
    table = pd.concat([
        sample_table(**{"TIS-a-1_peak_height": 20}),
        sample_table(**{"TIS-a-1_peak_height": 20}),
        sample_table(**{"TIS-a-1_peak_height": 20}),
    ], ignore_index=True)
    table["Identifier"] = ["chr:1-12:+", "chr:20-31:+", "chr:40-51:+"]
    table[column] = pd.array([True, None, False], dtype="boolean")
    path = tmp_path / f"sample{suffix}"
    if suffix == ".xlsx":
        merging.io.excel_writer(path, {"CDS": table})
    else:
        table.to_csv(path, sep="\t", index=False)

    metadata, dynamic = merging.screen_input_tables([path])
    result, _ = merging.build_merged_dataframe(metadata, dynamic)
    result = result.set_index("Identifier")

    assert result.at["chr:1-12:+", column] is True
    assert pd.isna(result.at["chr:20-31:+", column])
    assert result.at["chr:40-51:+", column] is False


@pytest.mark.parametrize("suffix", [".tsv", ".xlsx"])
@pytest.mark.parametrize(
    "column",
    ["TIS-a-1_peak_significant", "treated_vs_control_TIS_significant"],
)
def test_merge_rejects_invalid_statistical_booleans(tmp_path, suffix, column):
    table = sample_table(**{"TIS-a-1_peak_height": 20, column: 2})
    path = tmp_path / f"invalid{suffix}"
    if suffix == ".xlsx":
        merging.io.excel_writer(path, {"CDS": table})
    else:
        table.to_csv(path, sep="\t", index=False)

    with pytest.raises(ValueError) as error:
        merging.screen_input_tables([path])

    assert column in str(error.value)
    expected_value = "'2'" if suffix == ".tsv" else "2"
    assert f"got {expected_value}." in str(error.value)


@pytest.mark.parametrize("invalid", ["1.5", "1e3", "NaN", "inf", "True", "-1"])
def test_tabular_merge_rejects_non_decimal_statistical_counts(
        tmp_path, invalid):
    path = tmp_path / "invalid.tsv"
    sample_table(**{
        "TIS-a-1_peak_height": 20,
        "TIS-a-1_peak_count": invalid,
    }).to_csv(path, sep="\t", index=False)

    with pytest.raises(ValueError, match="nonnegative decimal integers"):
        merging.screen_input_tables([path])


def test_standalone_merge_preserves_attached_matched_evidence():
    prefix = "treated_vs_control_TIS"
    table = sample_table(**{
        "TIS-a-1_peak_height": 20,
        f"{prefix}_matched_pair_count": 2,
        f"{prefix}_numerator_peak_count": 9,
        f"{prefix}_common_odds_ratio": 4.5,
        f"{prefix}_p_value": 0.01,
        f"{prefix}_q_value": 0.02,
        f"{prefix}_significant": True,
        f"{prefix}_correction": "by",
        f"{prefix}_tests_in_family": 100,
        f"{prefix}_fit_status": "ok",
    })

    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})
    result, _ = merging.build_merged_dataframe(metadata, dynamic)

    for column in table.columns:
        if column.startswith(prefix):
            assert result.at[0, column] == table.at[0, column]


def test_standalone_merge_failure_does_not_publish_partial_output(
        tmp_path, monkeypatch):
    source = tmp_path / "sample.csv"
    sample_table(**{"sample_peak_height": 20}).to_csv(
        source, sep="\t", index=False,
    )
    output = tmp_path / "merged.xlsx"

    def fail_gff(*args, **kwargs):
        raise OSError("GFF write failed")

    monkeypatch.setattr(merging, "write_merged_gff", fail_gff)

    with pytest.raises(OSError, match="GFF write failed"):
        merging.merge_tables([source], output)

    assert source.is_file()
    assert not output.exists()
    assert not output.with_suffix(".csv").exists()
    assert not output.with_suffix(".gff").exists()
    assert not list(tmp_path.glob(".orfbounder-stage-*"))


def test_standalone_merge_rejects_non_xlsx_output_before_reading_inputs(
        tmp_path):
    output = tmp_path / "nested" / "merged.csv"

    with pytest.raises(ValueError, match=r"must use the \.xlsx extension"):
        merging.merge_tables([tmp_path / "missing.xlsx"], output)

    assert not output.parent.exists()


def test_supported_merge_cli_reports_input_failure_without_traceback(tmp_path):
    output = tmp_path / "nested" / "merged.xlsx"
    completed = subprocess.run(
        [
            sys.executable, "-m", "lib.merging", "-t",
            str(tmp_path / "missing.xlsx"), "-o", str(output),
        ],
        cwd=REPOSITORY,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert "ORFBounder merge error:" in completed.stderr
    assert "Traceback" not in completed.stderr
    assert not output.parent.exists()


def test_standalone_merge_never_overwrites_an_input_companion(tmp_path):
    source = tmp_path / "sample.csv"
    sample_table(**{"sample_peak_height": 20}).to_csv(
        source, sep="\t", index=False,
    )
    original = source.read_bytes()

    with pytest.raises(FileExistsError, match="Output already exists"):
        merging.merge_tables([source], source.with_suffix(".xlsx"))

    assert source.read_bytes() == original
    assert not source.with_suffix(".xlsx").exists()
    assert not source.with_suffix(".gff").exists()


def test_merged_tsv_round_trips_quoted_annotation_labels(tmp_path):
    label = 'gene\tline\n"quoted"'
    table = sample_table(**{"TIS-a-1_peak_height": 20})
    table.loc[0, "Locus_tag"] = label
    metadata, dynamic = merging.extend_combined_dictionary(table, {}, {})

    output = tmp_path / "merged.xlsx"
    merging.write_merged_table(metadata, dynamic, output)

    restored = pd.read_csv(output.with_suffix(".csv"), sep="\t")
    assert restored.loc[0, "Locus_tag"] == label
