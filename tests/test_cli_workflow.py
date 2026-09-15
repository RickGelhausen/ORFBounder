"""Configuration and complete-output regression tests using synthetic BAMs."""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import orfbounder as ob
from helpers.toy_example_generation import ToyExampleGenerator
from lib.config import check_config_sheet, discover_samples
from lib.io import check_alignment_path_input, parse_offset_json, parse_read_lengths
from orfbounder_batch import call_orfbounder


@pytest.fixture
def example(tmp_path):
    generator = ToyExampleGenerator(tmp_path / "example")
    generator.generate_all()
    return generator


def minimal_config(example, **overrides):
    row = {
        "experiment_name": "test", "annotation_file_path": str(example.data_dir / "annotation.gff"),
        "genome_file_path": str(example.data_dir / "genome.fa"), "TIS_folder_path": str(example.bam_dir),
        "offset_file_path": str(example.config_dir / "offsets.json"),
    }
    row.update(overrides)
    path = example.config_dir / "minimal.tsv"
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    return path


def test_optional_columns_and_mixed_blank_integer_values(example):
    path = minimal_config(example, min_peak_height="3")
    rows = pd.read_csv(path, sep="\t", dtype=str).to_dict("records")
    rows.append(dict(rows[0], experiment_name="second", min_peak_height=""))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    config = check_config_sheet(path)
    assert config.min_peak_height.tolist() == ["3", "5"]
    assert config.normalization_method.tolist() == ["raw", "raw"]
    assert config.read_length_json.tolist() == ["", ""]


@pytest.mark.parametrize("overrides, message", [
    ({"experiment_name": "../escape"}, "experiment_name"),
    ({"min_peak_height": "2.5"}, "integer"),
    ({"normalization_method": "min"}, "mapped_counts_file_path"),
    ({"mapping_method": "centered", "statistics": "local"}, "endpoint"),
    ({"fdr": "nan"}, "between"),
    ({"start_codons": "AXG"}, "triplets"),
    ({"start_codons": "TAA"}, "overlap"),
    ({"read_lengths": "30"}, "read_length_json"),
    ({"read_length_json": "/does/not/exist"}, "does not exist"),
    ({"normalisation_method": "raw"}, "Unknown config columns"),
    ({"min_mapq": "256"}, "min_mapq"),
    ({"duplicates": "sometimes"}, "duplicates"),
    ({"multimappers": "fractional"}, "multimappers"),
])
def test_invalid_configs_have_actionable_errors(example, overrides, message):
    with pytest.raises(ValueError, match=message):
        check_config_sheet(minimal_config(example, **overrides))


def test_batch_min_count_file_requires_every_discovered_calling_sample(example):
    counts = example.config_dir / "incomplete_counts.tsv"
    counts.write_text(f"TIS-WT-1\t{example.chrom}\t10\n")

    with pytest.raises(ValueError, match="missing samples: TIS-WT-2"):
        check_config_sheet(minimal_config(
            example,
            normalization_method="min",
            mapped_counts_file_path=str(counts),
        ))


def test_batch_rows_can_share_counts_with_suffix_aliases_and_references(
        example, tmp_path):
    first_folder = tmp_path / "first"
    second_folder = tmp_path / "second"
    first_folder.mkdir()
    second_folder.mkdir()
    (first_folder / "TIS-WT-1_aligned.bam").write_bytes(
        (example.bam_dir / "TIS-WT-1.bam").read_bytes()
    )
    (second_folder / "TIS-WT-2_filtered.bam").write_bytes(
        (example.bam_dir / "TIS-WT-2.bam").read_bytes()
    )
    counts = example.config_dir / "shared_counts.tsv"
    counts.write_text(
        f"TIS-WT-1\t{example.chrom}\t24\n"
        f"TIS-WT-2\t{example.chrom}\t20\n"
        f"TIS-reference-1\t{example.chrom}\t18\n"
    )
    common = {
        "annotation_file_path": str(example.data_dir / "annotation.gff"),
        "genome_file_path": str(example.data_dir / "genome.fa"),
        "offset_file_path": str(example.config_dir / "offsets.json"),
        "normalization_method": "min",
        "mapped_counts_file_path": str(counts),
    }
    config = example.config_dir / "shared-counts-config.tsv"
    pd.DataFrame([
        {**common, "experiment_name": "first", "TIS_folder_path": str(first_folder)},
        {**common, "experiment_name": "second", "TIS_folder_path": str(second_folder)},
    ]).to_csv(config, sep="\t", index=False)

    resolved = check_config_sheet(config)

    assert resolved.experiment_name.tolist() == ["first", "second"]


def test_validation_detects_reference_mismatch(example):
    path = minimal_config(example)
    (example.data_dir / "genome.fa").write_text(
        (example.data_dir / "genome.fa").read_text().replace("NC_000913.3", "wrong_assembly"))
    with pytest.raises(ValueError, match="same assembly"):
        check_config_sheet(path)


def test_discovery_is_exact_and_duplicate_safe(tmp_path):
    for name in ("TIS-WT-1.sam", "TIS-WT-10.bam", "RNATIS-WT-1.bam", "RIBO-only-1.bam"):
        (tmp_path / name).touch()
    groups = discover_samples(tmp_path, None, tmp_path)
    assert [wildcard for _, wildcard in groups] == ["WT-1", "WT-10"]
    assert groups[0][0][0].suffix == ".sam"
    (tmp_path / "TIS-WT-1_aligned.bam").touch()
    with pytest.raises(ValueError, match="Duplicate TIS-WT-1"):
        discover_samples(tmp_path, None, None)


def test_expression_discovery_does_not_mix_replicates(tmp_path):
    for name in ("TIS-WT-1.bam", "RNA-WT-1.sam", "TIS-WT-10.bam", "RNA-WT-10.bam"):
        (tmp_path / name).touch()
    result = check_alignment_path_input(tmp_path, tmp_path / "TIS-WT-1.bam", None)
    assert {path.name for path in result} == {"TIS-WT-1.bam", "RNA-WT-1.sam"}


def test_optional_lengths_and_invalid_json_schema(tmp_path):
    assert parse_read_lengths(None) is None
    path = tmp_path / "offsets.json"
    path.write_text('{"default": {"30": 1.5}}')
    with pytest.raises(ValueError, match="integer"):
        parse_offset_json(path)


def run_example(
    example, output_path, normalization="raw", statistics="none", threshold=1,
    tts=None, read_length_json=None,
):
    return ob.run_orfbounder(
        example.bam_dir / "TIS-WT-1.bam", tts, example.bam_dir / "RIBO-WT-1.bam",
        read_length_json, normalization, "fiveprime", example.data_dir / "annotation.gff",
        example.data_dir / "genome.fa", ["ATG"], ["TAG", "TAA", "TGA"], output_path,
        "WT-1", example.config_dir / "offsets.json", "furthest_inframe", threshold, "sum",
        True, example.bam_dir, None, statistics=statistics,
    )[0]


def test_direct_run_rejects_mixed_sample_identities_and_wrong_assay(
        example, tmp_path):
    mismatched_tts = tmp_path / "TTS-treated-2.bam"
    mismatched_tts.write_bytes((example.bam_dir / "TIS-WT-1.bam").read_bytes())

    mixed_output = tmp_path / "mixed-sample"
    with pytest.raises(ValueError, match="must share one condition and replicate"):
        run_example(example, mixed_output, tts=mismatched_tts)
    assert not mixed_output.exists()

    wrong_method_output = tmp_path / "wrong-method"
    with pytest.raises(ValueError, match="TTS input has RIBO sample identity"):
        run_example(
            example, wrong_method_output,
            tts=example.bam_dir / "RIBO-WT-1.bam",
        )
    assert not wrong_method_output.exists()


def test_direct_run_retains_custom_calling_filename_without_expression_folder(
        example, tmp_path):
    custom_alignment = tmp_path / "custom_start.bam"
    custom_alignment.write_bytes((example.bam_dir / "TIS-WT-1.bam").read_bytes())

    frame, _ = ob.run_orfbounder(
        custom_alignment, None, None, None, "raw", "fiveprime",
        example.data_dir / "annotation.gff", example.data_dir / "genome.fa",
        ["ATG"], ["TAG", "TAA", "TGA"], tmp_path / "custom-output",
        "custom", example.config_dir / "offsets.json", "furthest_inframe",
        1, "sum", True, None, None,
    )

    assert len(frame) == 3
    assert "custom_peak_height" in frame.columns


def test_read_length_filter_requires_each_calling_sample_or_default(example, tmp_path):
    lengths = example.config_dir / "typo_read_lengths.json"
    lengths.write_text(json.dumps({"TIS-WT-9": "28-30"}))

    with pytest.raises(ValueError, match="no entry or default.*TIS-WT-1, TIS-WT-2"):
        check_config_sheet(minimal_config(example, read_length_json=str(lengths)))
    with pytest.raises(ValueError, match="no entry or default.*RIBO-WT-1.*TIS-WT-1"):
        run_example(example, tmp_path / "rejected", read_length_json=lengths)
    assert not (tmp_path / "rejected").exists()

    lengths.write_text(json.dumps({"default": "28-30"}))
    assert len(check_config_sheet(
        minimal_config(example, read_length_json=str(lengths))
    )) == 1


def test_direct_run_validates_all_offset_samples_before_writing(example, tmp_path):
    offsets = example.config_dir / "incomplete_offsets.json"
    offsets.write_text(json.dumps({
        "TIS-WT-1": {"default": 0},
        "RIBO-WT-1": {"default": 0},
    }))
    tts = example.data_dir / "TTS-WT-1.sam"
    tts.write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:{example.chrom}\tLN:1000\n"
        f"stop\t0\t{example.chrom}\t160\t60\t30M\t*\t0\t0\t{'A' * 30}\t{'I' * 30}\n"
    )
    output = tmp_path / "rejected"

    with pytest.raises(ValueError, match="no entry or default.*TTS-WT-1"):
        ob.run_orfbounder(
            example.bam_dir / "TIS-WT-1.bam", tts,
            example.bam_dir / "RIBO-WT-1.bam", None, "raw", "fiveprime",
            example.data_dir / "annotation.gff", example.data_dir / "genome.fa",
            ["ATG"], ["TAG", "TAA", "TGA"], output, "WT-1", offsets,
            "furthest_inframe", 1, "sum", True, example.bam_dir, None,
        )

    assert not output.exists()


def test_specific_expression_lengths_require_expression_sample_entries(example):
    lengths = example.config_dir / "calling_only_read_lengths.json"
    lengths.write_text(json.dumps({"TIS-WT-1": "28-30", "TIS-WT-2": "28-30"}))
    settings = {
        "read_length_json": str(lengths),
        "alignment_folder_path": str(example.bam_dir),
    }

    assert len(check_config_sheet(minimal_config(example, **settings))) == 1
    with pytest.raises(ValueError, match="no entry or default.*RIBO-WT-1, RIBO-WT-2"):
        check_config_sheet(minimal_config(
            example, **settings, rpkm_read_usage="specific",
        ))


def test_fresh_example_produces_three_annotated_orfs_and_manifest(example, tmp_path):
    frame = run_example(example, tmp_path / "run")
    assert frame.Start.tolist() == [100, 300, 500]
    assert frame.Stop.tolist() == [162, 362, 562]
    assert frame.Type.tolist() == ["Annotated"] * 3
    assert frame.Codon_count.tolist() == [21] * 3
    manifest = json.loads((tmp_path / "run" / "WT-1.run.json").read_text())
    assert manifest["orf_count"] == 3
    assert manifest["read_lengths"] is None
    assert manifest["alignment_diagnostics"]["TIS-WT-1"]["accepted_reads_by_contig"] == {example.chrom: 24}
    assert manifest["alignment_policy"] == {
        "min_mapq": 0, "duplicates": "include", "multimappers": "exclude",
        "count_unit": "alignment_record",
    }
    with pytest.raises(FileExistsError):
        run_example(example, tmp_path / "run")


def test_statistics_are_independent_of_normalization(example, tmp_path):
    raw = run_example(example, tmp_path / "raw", statistics="local")
    mil = run_example(example, tmp_path / "mil", normalization="mil", statistics="local")
    assert raw.Identifier.tolist() == mil.Identifier.tolist()
    assert raw["TIS-WT-1_peak_q_value"].tolist() == mil["TIS-WT-1_peak_q_value"].tolist()
    assert raw["TIS-WT-1_peak_fold_enrichment"].tolist() == mil[
        "TIS-WT-1_peak_fold_enrichment"
    ].tolist()
    assert raw["TIS-WT-1_peak_fold_enrichment"].gt(1).all()
    assert raw["TIS-WT-1_peak_fit_status"].eq("ok").all()
    assert raw["TIS-WT-1_peak_count"].tolist() == [10, 8, 6]
    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "raw/peak_statistics/WT-1_TIS-WT-1.tsv", sep="\t"),
        pd.read_csv(tmp_path / "mil/peak_statistics/WT-1_TIS-WT-1.tsv", sep="\t"),
    )
    manifest = json.loads((tmp_path / "raw/WT-1.run.json").read_text())
    assert manifest["statistics_model"] == "conditional_binomial_local_endpoint_count"
    assert manifest["statistics_alternative"] == "greater"
    assert manifest["statistics_minimum_positive_p_value"] > 0
    assert manifest["statistics_maximum_supported_total_count"] == (1 << 31) - 1
    assert manifest["statistics_max_tail_recurrence_work"] == 1_000_000
    assert manifest["statistics_max_exact_combination_work"] == 100_000


def test_combined_tis_tts_keeps_both_assays_and_statistical_families(example, tmp_path):
    tts = example.data_dir / "TTS-WT-1.sam"
    tts.write_text(f"@HD\tVN:1.6\n@SQ\tSN:{example.chrom}\tLN:1000\n" + "".join(
        f"stop{stop}_{i}\t0\t{example.chrom}\t{stop - 2}\t60\t30M\t*\t0\t0\t{'A' * 30}\t{'I' * 30}\n"
        for stop in (162, 362, 562) for i in range(8)
    ))
    frame = run_example(example, tmp_path / "combined", statistics="local", tts=tts)
    assert len(frame) == 3
    assert frame["TTS-WT-1_peak_height"].tolist() == [8] * 3
    assert frame["TTS-WT-1_peak_count"].tolist() == [8] * 3
    assert frame["TIS-WT-1_peak_count"].tolist() == [10, 8, 6]
    assert frame["TTS-WT-1_peak_q_value"].notna().all()
    assert len(list((tmp_path / "combined/peak_statistics").glob("*.tsv"))) == 2


def test_expression_duplicates_fail_during_validation(example):
    duplicate = example.bam_dir / "RIBO-WT-1_copy.bam"
    duplicate.write_bytes((example.bam_dir / "RIBO-WT-1.bam").read_bytes())
    config = minimal_config(example, alignment_folder_path=str(example.bam_dir))
    with pytest.raises(ValueError, match="Duplicate expression alignments"):
        check_config_sheet(config)


def test_fully_empty_alignments_are_valid_zero_call_results(example, tmp_path):
    import pysam
    for path in example.bam_dir.glob("*.bam"):
        with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": example.chrom, "LN": 1000}]}):
            pass
        Path(str(path) + ".bai").unlink()
    frame = run_example(example, tmp_path / "empty", statistics="local")
    assert frame.empty
    candidates = pd.read_csv(tmp_path / "empty/peak_statistics/WT-1_TIS-WT-1.tsv", sep="\t")
    assert len(candidates) > 0
    assert candidates.q_value.eq(1).all()


def test_batch_emits_empty_merged_outputs(example, tmp_path):
    config = check_config_sheet(minimal_config(example, min_peak_height="999999"))
    out = tmp_path / "out"
    call_orfbounder(config, out)
    final = out / "test/threeprime/raw/test_final.csv"
    assert pd.read_csv(final, sep="\t").empty
    assert final.with_suffix(".gff").read_text().startswith("##gff-version 3")
    assert final.with_suffix(".xlsx").is_file()


def test_cli_validation_and_failure_exit_codes(example):
    config = minimal_config(example)
    command = [sys.executable, str(Path(ob.__file__).with_name("orfbounder_batch.py")), "-c", str(config)]
    result = subprocess.run(command + ["--validate-only"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Configuration valid" in result.stdout
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 2
    assert "result_path is required" in result.stderr
    config.write_text("wrong\theader\na\tb\n")
    result = subprocess.run(command + ["--validate-only"], capture_output=True, text=True)
    assert result.returncode == 2
    assert "Missing required columns" in result.stderr
    assert "Traceback" not in result.stderr
