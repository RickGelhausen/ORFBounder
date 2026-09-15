"""Verify the files users receive, including failures after analysis succeeds."""

import json
from pathlib import Path
import sys

import pandas as pd
import pytest

import orfbounder
from orfbounder_batch import call_orfbounder
from lib import config, io
from helpers.toy_example_generation import ToyExampleGenerator


@pytest.fixture
def example(tmp_path):
    generator = ToyExampleGenerator(tmp_path / "example")
    generator.generate_all()
    return generator


def cli_args(example, destination, basename="WT-1"):
    return ["orfbounder", "--alignment_file_tis", str(example.bam_dir / "TIS-WT-1.bam"),
            "--annotation_file", str(example.data_dir / "annotation.gff"),
            "--genome_file", str(example.data_dir / "genome.fa"),
            "--offset_json", str(example.config_dir / "offsets.json"),
            "--mapping_method", "fiveprime", "--min_peak_height", "1", "--statistics", "local",
            "--output_basename", basename, "--output_path", str(destination)]


def test_direct_cli_publishes_portable_report_sequences_and_complete_manifest(example, tmp_path, monkeypatch):
    destination = tmp_path / "results"
    monkeypatch.setattr(sys, "argv", cli_args(example, destination))
    orfbounder.main()
    assert (destination / "WT-1.report.html").is_file()
    complete = json.loads((destination / "WT-1.complete.json").read_text())
    for output in complete["outputs"]:
        assert (destination / output).is_file()
        assert ".orfbounder-stage-" not in output
    manifest = json.loads((destination / "WT-1.run.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["sample_names"] == {"TIS": "TIS-WT-1"}
    assert all(len(item["sha256"]) == 64 for item in manifest["input_files"])
    assert (destination / "sequence_per_sample/WT-1.fna").read_text().count(">") == 3
    assert (destination / "sequence_per_sample/WT-1.faa").read_text().count(">") == 3
    # A second basename can share the folder without changing the first run.
    original = (destination / "WT-1.run.json").read_bytes()
    monkeypatch.setattr(sys, "argv", cli_args(example, destination, "second"))
    orfbounder.main()
    assert (destination / "second.complete.json").is_file()
    assert (destination / "WT-1.run.json").read_bytes() == original


def test_direct_cli_cleans_up_when_spreadsheet_export_fails(example, tmp_path, monkeypatch):
    destination = tmp_path / "results"
    destination.mkdir()
    (destination / "notes.txt").write_text("my notes")
    monkeypatch.setattr(sys, "argv", cli_args(example, destination))

    def fail_export(*args, **kwargs):
        raise OSError("spreadsheet write failed")

    monkeypatch.setattr(io, "write_results_to_table", fail_export)
    with pytest.raises(SystemExit) as exc:
        orfbounder.main()
    assert exc.value.code == 2
    assert {path.name for path in destination.iterdir()} == {"notes.txt"}
    assert (destination / "notes.txt").read_text() == "my notes"


def test_batch_report_and_fasta_match_final_table(example, tmp_path):
    config_path = example.config_dir / "config.tsv"
    frame = pd.read_csv(config_path, sep="\t", keep_default_na=False)
    frame["statistics"] = "local"
    frame.to_csv(config_path, sep="\t", index=False)
    call_orfbounder(config.check_config_sheet(config_path), tmp_path / "results")
    folder = tmp_path / "results/toy_example/fiveprime/raw"
    table = pd.read_csv(folder / "toy_example_final.csv", sep="\t")
    assert len(table) == 3
    assert (folder / "report.html").is_file()
    assert (folder / "toy_example_final.fna").read_text().count(">") == len(table)
    complete = json.loads((folder / "complete.json").read_text())
    assert complete["samples"] == ["WT-1", "WT-2"]
    assert all((folder / file).is_file() for file in complete["outputs"])
    for name in complete["samples"]:
        record = json.loads((folder / f"{name}.run.json").read_text())
        assert record["status"] == "complete"
        assert str(config_path) in {path["path"] for path in record["input_files"]}


def test_batch_config_mutation_is_rejected_before_output(example, tmp_path):
    config_path = example.config_dir / "config.tsv"
    frame = config.check_config_sheet(config_path)
    config_path.write_text(config_path.read_text().replace("toy_example\t", "changed\t"))
    with pytest.raises(OSError, match="changed during analysis"):
        call_orfbounder(frame, tmp_path / "results")
    assert not (tmp_path / "results").exists()


def test_changed_offsets_cannot_be_recorded_with_old_parsed_values(example, tmp_path, monkeypatch):
    original = io.parse_offset_json

    def mutate_after_parsing(path):
        parsed = original(path)
        Path(path).write_text('{"default":{"default":-1}}')
        return parsed

    monkeypatch.setattr(io, "parse_offset_json", mutate_after_parsing)
    destination = tmp_path / "results"
    monkeypatch.setattr(sys, "argv", cli_args(example, destination))
    with pytest.raises(SystemExit) as exc:
        orfbounder.main()
    assert exc.value.code == 2
    assert not destination.exists()


def test_batch_failure_after_first_sample_does_not_publish_partial_experiment(example, tmp_path, monkeypatch):
    original = orfbounder.run_orfbounder
    calls = 0

    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("second sample failed")
        return original(*args, **kwargs)

    monkeypatch.setattr(orfbounder, "run_orfbounder", fail_second)
    destination = tmp_path / "results"
    with pytest.raises(OSError, match="second sample failed"):
        call_orfbounder(config.check_config_sheet(example.config_dir / "config.tsv"), destination)
    assert not [path for path in destination.rglob("*") if path.is_file()]
    assert not list(destination.rglob(".orfbounder-stage-*"))
