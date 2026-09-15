"""Exercise count recovery with the same sample selection and reads as peaks."""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from helpers import recover_total_counts
from helpers.recover_total_counts import check_config_sheet, count_mapped_reads, recover_read_information
from lib.alignment_reader import AlignmentPolicy, PositionReader
from lib import run_output
from tests.test_core_regressions import write_alignment


def write_config(tmp_path, folder, **extra):
    path = tmp_path / "config.tsv"
    pd.DataFrame([{"experiment_name": "counts", "TIS_folder_path": str(folder), **extra}]).to_csv(path, sep="\t", index=False)
    return path


def test_minimal_config_recovers_only_selected_samples_and_creates_directories(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    primary = write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    # These unrelated libraries must not influence the minimum count.
    (inputs / "RNA-control-1.sam").write_text(primary.read_text())
    (inputs / "RIBO-other-1.sam").write_text(primary.read_text())
    config = check_config_sheet(write_config(tmp_path, inputs, normalization_method="min", mapped_counts_file_path="not-created-yet.tsv"))
    assert config.loc[0, "read_length_json"] == ""
    result = tmp_path / "nested" / "output"
    recover_read_information(config, result)
    assert (result / "counts_mapped_reads.tsv").read_text() == "TIS-control-1\tchr1\t1\n"
    with pytest.raises(FileExistsError, match="Output already exists"):
        recover_read_information(config, result)


def test_count_recovery_uses_query_length_and_optional_nh_primary_filter(tmp_path):
    path = write_alignment(tmp_path, [
        ("query_length_28", 0, 20, "5S10M2D10M3S", 28, None),
        ("query_length_22", 0, 40, "22M", 22, None),
        ("secondary", 256, 60, "28M", 28, None),
        ("multimapped", 0, 90, "28M", 28, 2),
    ], suffix=".bam")
    result = count_mapped_reads(path, {}, {"default": ["28"]})
    assert result == {("TIS-control-1", "chr1"): 1}


def test_count_recovery_infers_query_length_when_sam_sequence_is_omitted(tmp_path):
    path = write_alignment(tmp_path, [
        ("sequence_omitted", 0, 20, "5S10M2D10M3S", None, None),
        ("different_length", 0, 60, "22M", None, None),
    ])

    result = count_mapped_reads(path, {}, {"default": ["28"]})

    assert result == {("TIS-control-1", "chr1"): 1}


def test_count_recovery_uses_the_same_alignment_policy_and_diagnostics_as_peaks(tmp_path):
    path = write_alignment(tmp_path, [
        ("low_mapq", 0, 20, "10M", 10, 1, 5),
        ("accepted", 0, 40, "10M", 10, 1, 30),
        ("duplicate", 1024, 60, "10M", 10, 1, 30),
        ("multimapped", 0, 80, "10M", 10, 2, 30),
        ("without_nh", 0, 100, "10M", 10, None, 30),
    ])
    policy = AlignmentPolicy(min_mapq=20, include_duplicates=False, multimappers="include")
    reader = PositionReader(
        path, None, "threeprime", {"default": {"default": 0}}, alignment_policy=policy)
    diagnostics = {}

    recovered = count_mapped_reads(
        path, {}, None, alignment_policy=policy, diagnostics=diagnostics)

    assert recovered == {("TIS-control-1", "chr1"): 3}
    assert reader.output()[1] == {"chr1": 3}
    assert diagnostics == reader.alignment_diagnostics


def test_optional_read_length_path_is_converted_to_path(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("ten", 0, 20, "10M", 10, None), ("twenty", 0, 40, "20M", 20, None)])
    lengths = tmp_path / "lengths.json"
    lengths.write_text(json.dumps({"default": 10}))
    config = check_config_sheet(write_config(tmp_path, inputs, read_length_json=str(lengths)))
    recover_read_information(config, tmp_path / "output")
    assert (tmp_path / "output" / "counts_mapped_reads.tsv").read_text().endswith("\t1\n")


def test_recovery_requires_lengths_for_every_counted_alignment_sample(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    tis = write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    (inputs / "TTS-control-1.sam").write_text(tis.read_text())
    (inputs / "RIBO-control-1.sam").write_text(tis.read_text())
    lengths = tmp_path / "lengths.json"
    lengths.write_text(json.dumps({"TIS-control-1": 10}))
    config = write_config(
        tmp_path,
        inputs,
        TTS_folder_path=str(inputs),
        RIBO_folder_path=str(inputs),
        read_length_json=str(lengths),
    )

    with pytest.raises(
        ValueError,
        match=(
            "Read-length JSON has no entry or default for: "
            "RIBO-control-1, TTS-control-1"
        ),
    ):
        check_config_sheet(config)


def test_empty_experiment_fails_before_writing_counts(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [])
    config = check_config_sheet(write_config(tmp_path, inputs))
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="no accepted mapped reads"):
        recover_read_information(config, output)
    assert not output.exists()


def test_multiple_experiments_do_not_leave_partial_count_files_on_write_failure(
        tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    config_path = write_config(tmp_path, inputs)
    config = pd.read_csv(config_path, sep="\t", dtype=str, keep_default_na=False)
    config = pd.concat([
        config.assign(experiment_name="first"),
        config.assign(experiment_name="second"),
    ], ignore_index=True)
    config.to_csv(config_path, sep="\t", index=False)
    validated = check_config_sheet(config_path)
    original_write = recover_total_counts.write_read_counts_to_file
    calls = 0

    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("disk full")
        return original_write(*args, **kwargs)

    monkeypatch.setattr(
        recover_total_counts, "write_read_counts_to_file", fail_second,
    )
    output = tmp_path / "output"

    with pytest.raises(OSError, match="disk full"):
        recover_read_information(validated, output)

    assert not list(output.glob("*_mapped_reads.tsv"))
    assert not list(output.glob(".orfbounder-stage-*"))


def test_count_recovery_retries_an_interrupted_owned_output(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    validated = check_config_sheet(write_config(tmp_path, inputs))
    output = tmp_path / "output"
    output.mkdir()
    stage = output / ".orfbounder-stage-interrupted"
    stage.mkdir()
    (stage / ".transaction.lock").touch()
    source = stage / "counts_mapped_reads.tsv"
    source.write_text("partial\n")
    os.link(source, output / source.name)
    run_output._write_transaction_journal(stage, [{
        "relative": source.name,
        "completion": False,
    }])

    recover_read_information(validated, output)

    assert (output / source.name).read_text() == "TIS-control-1\tchr1\t1\n"
    assert not stage.exists()


def test_count_recovery_rejects_broken_output_symlink_before_counting(
        tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    validated = check_config_sheet(write_config(tmp_path, inputs))
    output = tmp_path / "output"
    output.mkdir()
    target = output / "counts_mapped_reads.tsv"
    target.symlink_to(tmp_path / "missing.tsv")

    def should_not_count(*args, **kwargs):
        raise AssertionError("existing outputs must fail before counting")

    monkeypatch.setattr(recover_total_counts, "count_mapped_reads", should_not_count)

    with pytest.raises(FileExistsError, match="Output already exists"):
        recover_read_information(validated, output)
    assert target.is_symlink()


def test_count_recovery_discards_results_when_an_input_changes_during_counting(
        tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    alignment = write_alignment(
        inputs, [("accepted", 0, 20, "10M", 10, None)],
    )
    validated = check_config_sheet(write_config(tmp_path, inputs))
    original_count = recover_total_counts.count_mapped_reads

    def mutate_after_counting(*args, **kwargs):
        counts = original_count(*args, **kwargs)
        alignment.write_text(alignment.read_text() + "\n")
        return counts

    monkeypatch.setattr(
        recover_total_counts, "count_mapped_reads", mutate_after_counting,
    )
    output = tmp_path / "output"

    with pytest.raises(OSError, match="changed during analysis"):
        recover_read_information(validated, output)

    assert not output.exists()


def test_count_recovery_rejects_a_changed_validated_config_before_counting(
        tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    config_path = write_config(tmp_path, inputs)
    validated = check_config_sheet(config_path)
    config_path.write_text(config_path.read_text() + "\n")

    def should_not_count(*args, **kwargs):
        raise AssertionError("changed config must fail before counting")

    monkeypatch.setattr(recover_total_counts, "count_mapped_reads", should_not_count)

    with pytest.raises(OSError, match="changed during analysis"):
        recover_read_information(validated, tmp_path / "output")


def test_count_recovery_rejects_unknown_duplicate_policy(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_alignment(inputs, [("accepted", 0, 20, "10M", 10, None)])
    config = write_config(tmp_path, inputs, duplicates="sometimes")
    with pytest.raises(ValueError, match="duplicates must be 'include' or 'exclude'"):
        check_config_sheet(config)


def test_count_recovery_module_cli_reports_invalid_config_without_traceback(tmp_path):
    config = tmp_path / "missing-experiment.tsv"
    config.write_text("TIS_folder_path\n/tmp\n")
    result = subprocess.run([sys.executable, "-m", "helpers.recover_total_counts", "-c", str(config), "-r", str(tmp_path / "output")],
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert "count recovery error:" in result.stderr
    assert "Traceback" not in result.stderr
