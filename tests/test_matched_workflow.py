"""Batch integration for explicit matched-condition statistical comparisons."""

import json
import os
from pathlib import Path
import shutil

import pandas as pd
import pytest

from helpers.toy_example_generation import ToyExampleGenerator
from lib.config import check_config_sheet
from orfbounder_batch import call_orfbounder


@pytest.fixture
def matched_example(tmp_path):
    example = ToyExampleGenerator(tmp_path / "example")
    example.generate_all()
    alignments = tmp_path / "matched-alignments"
    alignments.mkdir()
    for condition, sources in {
        "treated": ("TIS-WT-1.bam", "TIS-WT-2.bam"),
        "control": ("RIBO-WT-1.bam", "RIBO-WT-2.bam"),
    }.items():
        for replicate, source in enumerate(sources, 1):
            shutil.copyfile(example.bam_dir / source, alignments / f"TIS-{condition}-{replicate}.bam")

    comparisons = tmp_path / "comparisons.tsv"
    comparisons.write_text(
        "comparison\tassay\tnumerator_condition\tdenominator_condition\n"
        "treated_vs_control\tTIS\ttreated\tcontrol\n"
    )
    config = tmp_path / "config.tsv"
    pd.DataFrame([{
        "experiment_name": "matched", "annotation_file_path": example.data_dir / "annotation.gff",
        "genome_file_path": example.data_dir / "genome.fa", "TIS_folder_path": alignments,
        "offset_file_path": example.config_dir / "offsets.json", "mapping_method": "fiveprime",
        "statistics": "local", "min_peak_height": "1",
        "matched_comparison_file_path": comparisons,
    }]).to_csv(config, sep="\t", index=False)
    return example, alignments, comparisons, config


def test_matched_comparison_validates_replicates_and_requires_local_statistics(matched_example):
    _, alignments, _, config = matched_example
    check_config_sheet(config)

    frame = pd.read_csv(config, sep="\t", dtype=str)
    frame.loc[0, "statistics"] = "none"
    frame.to_csv(config, sep="\t", index=False)
    with pytest.raises(ValueError, match="requires statistics=local"):
        check_config_sheet(config)

    frame.loc[0, "statistics"] = "local"
    (alignments / "TIS-control-2.bam").unlink()
    frame.to_csv(config, sep="\t", index=False)
    with pytest.raises(ValueError, match="identical replicate IDs"):
        check_config_sheet(config)


def test_batch_exports_and_attaches_matched_condition_evidence(matched_example, tmp_path):
    _, _, comparisons, config = matched_example
    output = tmp_path / "results"
    call_orfbounder(check_config_sheet(config), output)

    folder = output / "matched/fiveprime/raw"
    candidate_path = folder / "matched_statistics/treated_vs_control_TIS.tsv"
    strata_path = folder / "matched_statistics/treated_vs_control_TIS_strata.tsv"
    candidates = pd.read_csv(candidate_path, sep="\t")
    strata = pd.read_csv(strata_path, sep="\t")
    assert len(candidates) > 0
    assert set(candidates.matched_pair_count) == {2}
    assert set(candidates.assay) == {"TIS"}
    assert set(candidates.numerator_condition) == {"treated"}
    assert set(candidates.denominator_condition) == {"control"}
    assert len(strata) == 2 * len(candidates)
    assert set(strata.replicate_id.astype(str)) == {"1", "2"}
    assert {
        "numerator_peak_count", "numerator_background_count",
        "denominator_peak_count", "denominator_background_count",
        "informative", "pair_direction",
    } <= set(strata.columns)

    merged = pd.read_csv(folder / "matched_final.csv", sep="\t")
    prefix = "treated_vs_control_TIS"
    assert f"{prefix}_common_odds_ratio" in merged
    assert f"{prefix}_q_value" in merged
    assert f"{prefix}_significant" in merged
    assert "TIS-treated-1_peak_fit_status" in merged

    complete = json.loads((folder / "complete.json").read_text())
    assert "matched_statistics/treated_vs_control_TIS.tsv" in complete["outputs"]
    manifest = json.loads((folder / "treated-1.run.json").read_text())
    assert manifest["matched_comparisons"][0]["comparison"] == "treated_vs_control"
    matched_metadata = manifest["matched_comparisons"][0]
    assert matched_metadata["model"] == "exact_stratified_fixed_margin_endpoint_count"
    assert matched_metadata["numeric_backend"] in {
        "extended_longdouble", "float64_with_underflow_guard",
    }
    assert matched_metadata["alternative"] == "greater"
    assert matched_metadata["replicate_ids"] == ["1", "2"]
    assert matched_metadata["minimum_informative_pairs"] == 2
    assert "2*concordant_pair_count>informative_pair_count" in matched_metadata["support_rule"]
    assert matched_metadata["maximum_supported_hypergeometric_total"] == (1 << 31) - 1
    assert matched_metadata["strata_table"].endswith("_strata.tsv")
    assert str(comparisons.resolve()) in {record["path"] for record in manifest["input_files"]}
    assert "matched_statistics/treated_vs_control_TIS.tsv" in (folder / "files.html").read_text()
    assert "matched_statistics/treated_vs_control_TIS_strata.tsv" in complete["outputs"]
    assert "matched_statistics/treated_vs_control_TIS_strata.tsv" in (folder / "files.html").read_text()


def test_batch_rejects_comparison_file_changed_while_it_is_parsed(
        matched_example, tmp_path, monkeypatch):
    _, _, comparisons, config = matched_example
    import lib.matched_statistics as matched_statistics

    original = matched_statistics.read_matched_comparisons

    def parse_then_mutate(path):
        result = original(path)
        before = Path(path).stat()
        Path(path).write_text(Path(path).read_text() + "\n")
        os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
        return result

    monkeypatch.setattr(matched_statistics, "read_matched_comparisons", parse_then_mutate)
    output = tmp_path / "results"
    with pytest.raises(OSError, match="changed during analysis"):
        call_orfbounder(check_config_sheet(config), output)
    assert not list(output.rglob("complete.json"))


def test_batch_rejects_comparison_file_changed_after_parsing(
        matched_example, tmp_path, monkeypatch):
    _, _, comparisons, config = matched_example
    import orfbounder_batch

    original = orfbounder_batch.ob.run_orfbounder
    calls = 0

    def run_then_mutate(*args, **kwargs):
        nonlocal calls
        result = original(*args, **kwargs)
        calls += 1
        if calls == 1:
            comparisons.write_text(comparisons.read_text() + "\n")
        return result

    monkeypatch.setattr(orfbounder_batch.ob, "run_orfbounder", run_then_mutate)
    output = tmp_path / "results"
    with pytest.raises(OSError, match="changed during analysis"):
        call_orfbounder(check_config_sheet(config), output)
    assert not list(output.rglob("complete.json"))
