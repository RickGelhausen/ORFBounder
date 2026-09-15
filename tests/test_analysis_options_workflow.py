"""Sequence/code and scaling options through real input validation and ORF search."""

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pysam
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from lib.config import check_config_sheet
from orfbounder import run_orfbounder


@pytest.fixture
def mixed_library(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    sequences = {"chr": list("C" * 200), "plasmid": list("C" * 200)}
    sequences["chr"][30:42] = "GTGAAATGATAA"
    sequences["chr"][100:112] = str(Seq("TTGAAATGATAA").reverse_complement())
    sequences["plasmid"][60:69] = "ATGCCCTAA"
    (inputs / "genome.fa").write_text("".join(f">{name}\n{''.join(seq)}\n" for name, seq in sequences.items()))
    (inputs / "annotation.gff").write_text(
        "##gff-version 3\n"
        "chr\ttest\tCDS\t31\t42\t.\t+\t0\tID=plus;locus_tag=plus\n"
        "chr\ttest\tCDS\t101\t112\t.\t-\t0\tID=minus;locus_tag=minus\n"
        "plasmid\ttest\tCDS\t61\t69\t.\t+\t0\tID=plasmid;locus_tag=plasmid\n"
    )
    (inputs / "offsets.json").write_text('{"default": {"default": 0}}\n')
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": name, "LN": 200} for name in sequences]}
    with pysam.AlignmentFile(inputs / "TIS-WT-1.sam", "w", header=header) as output:
        for index, (reference, position, flag) in enumerate([(0, 30, 0)] * 3 + [(0, 102, 16)] * 2 + [(1, 60, 0)]):
            read = pysam.AlignedSegment(output.header)
            read.query_name = f"read{index}"
            read.reference_id = reference
            read.reference_start = position
            read.flag = flag
            read.mapping_quality = 30
            read.cigarstring = "10M"
            read.query_sequence = "A" * 10
            output.write(read)
    return inputs


def analyze(inputs, output, **settings):
    total_read_file_path = settings.pop("total_read_file_path", None)
    min_peak_height = settings.pop("min_peak_height", 0)
    return run_orfbounder(
        inputs / "TIS-WT-1.sam", None, None, None, settings.pop("normalization", "raw"),
        "fiveprime", inputs / "annotation.gff", inputs / "genome.fa", None, None,
        output, "WT-1", inputs / "offsets.json", "furthest_inframe", min_peak_height, "sum",
        True, inputs, total_read_file_path, **settings,
    )[0]


@pytest.mark.parametrize("invalid_fdr", [True, float("nan"), "0.05"])
def test_direct_run_rejects_invalid_fdr_before_output(
        mixed_library, tmp_path, invalid_fdr):
    output = tmp_path / "invalid-fdr"
    with pytest.raises(ValueError, match="Statistical settings"):
        analyze(mixed_library, output, fdr=invalid_fdr)
    assert not output.exists()


@pytest.mark.parametrize("setting,value,message", [
    ("min_peak_height", True, "nonnegative integer"),
    ("min_peak_height", 1.5, "nonnegative integer"),
    ("background_width", True, "Statistical settings"),
    ("background_width", 1.5, "Statistical settings"),
])
def test_direct_run_rejects_noninteger_count_settings_before_output(
        mixed_library, tmp_path, setting, value, message):
    output = tmp_path / "invalid-count-setting"
    with pytest.raises(ValueError, match=message):
        analyze(mixed_library, output, **{setting: value})
    assert not output.exists()


def test_genetic_code_changes_boundaries_and_translation_on_both_strands(mixed_library, tmp_path):
    bacterial = analyze(mixed_library, tmp_path / "code11")
    alternative = analyze(mixed_library, tmp_path / "code4", genetic_code=4)
    assert list(zip(bacterial.Start, bacterial.Stop)) == [(31, 39), (104, 112), (61, 69)]
    assert list(zip(alternative.Start, alternative.Stop)) == [(31, 42), (101, 112), (61, 69)]
    assert bacterial.Amino_Acid_Seq.tolist() == ["MK*", "MK*", "MP*"]
    assert alternative.Amino_Acid_Seq.tolist() == ["MKW*", "MKW*", "MP*"]
    metadata = json.loads((tmp_path / "code4/WT-1.run.json").read_text())
    assert metadata["genetic_code"] == 4
    assert "ATA" in metadata["start_codons"]
    assert set(metadata["stop_codons"]) == {"TAA", "TAG"}


def test_local_statistics_invariant_across_library_and_contig_scaling(mixed_library, tmp_path):
    outputs = {}
    for normalization, scope in [("raw", "contig"), ("raw", "library"), ("mil", "contig"), ("mil", "library")]:
        directory = tmp_path / f"{normalization}_{scope}"
        frame = analyze(mixed_library, directory, genetic_code=4, normalization=normalization,
                        normalization_scope=scope, statistics="local")
        outputs[normalization, scope] = frame
        statistics = pd.read_csv(directory / "peak_statistics/WT-1_TIS-WT-1.tsv", sep="\t")
        if (normalization, scope) == ("raw", "contig"):
            expected_statistics = statistics
        else:
            pd.testing.assert_frame_equal(statistics, expected_statistics)
    assert outputs["raw", "library"]["TIS-WT-1_peak_height"].tolist() == [3, 2, 1]
    assert outputs["mil", "library"]["TIS-WT-1_peak_height"].tolist() == pytest.approx([500000, 1000000 / 3, 1000000 / 6])
    assert outputs["mil", "contig"]["TIS-WT-1_peak_height"].tolist() == [600000, 400000, 1000000]
    assert outputs["raw", "library"]["TIS-WT-1_rpkm"].tolist() == [round(3e9 / 72, 2), round(2e9 / 72, 2), round(1e9 / 54, 2)]
    for scope in ("library", "contig"):
        pd.testing.assert_series_equal(outputs["raw", scope]["TIS-WT-1_rpkm"], outputs["mil", scope]["TIS-WT-1_rpkm"])


def test_batch_config_resolves_code_specific_codons_and_library_min(mixed_library, tmp_path):
    config = tmp_path / "runs.tsv"
    row = {"experiment_name": "test", "annotation_file_path": str(mixed_library / "annotation.gff"),
           "genome_file_path": str(mixed_library / "genome.fa"), "offset_file_path": str(mixed_library / "offsets.json"),
           "TIS_folder_path": str(mixed_library), "genetic_code": "4", "normalization_scope": "library"}
    pd.DataFrame([row]).to_csv(config, sep="\t", index=False)
    resolved = check_config_sheet(config).iloc[0]
    assert set(resolved.stop_codons.split(",")) == {"TAA", "TAG"}
    assert "ATA" in resolved.start_codons.split(",")
    row["normalization_method"] = "min"
    counts = tmp_path / "counts.tsv"
    counts.write_text("TIS-WT-1\tchr1\t6\n")
    row["mapped_counts_file_path"] = str(counts)
    pd.DataFrame([row]).to_csv(config, sep="\t", index=False)
    resolved = check_config_sheet(config).iloc[0]
    assert resolved.normalization_scope == "library"


def test_library_min_uses_the_smallest_whole_sample_total(mixed_library, tmp_path):
    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "TIS-WT-1\tchr\t5\n"
        "TIS-WT-1\tplasmid\t1\n"
        "TIS-other-1\tchr\t3\n"
        "TIS-other-1\tplasmid\t1\n"
    )
    output = tmp_path / "library_min"
    frame = analyze(
        mixed_library, output, genetic_code=4, normalization="min",
        normalization_scope="library", total_read_file_path=counts,
    )
    assert frame["TIS-WT-1_peak_height"].tolist() == pytest.approx([2, 4 / 3, 2 / 3])
    metadata = json.loads((output / "WT-1.run.json").read_text())
    assert metadata["minimum_counts"] == {"__library_minimum__": 4}


def test_direct_min_requires_current_calling_sample_in_count_file(
        mixed_library, tmp_path):
    counts = tmp_path / "counts.tsv"
    counts.write_text("TIS-other-1\tchr\t4\n")
    output = tmp_path / "missing_min_sample"

    with pytest.raises(ValueError, match="missing samples: TIS-WT-1"):
        analyze(
            mixed_library, output, normalization="min",
            total_read_file_path=counts,
        )
    assert not output.exists()


def test_batch_normalization_hint_cannot_replace_current_calling_sample(
        mixed_library, tmp_path):
    counts = tmp_path / "reference-only-counts.tsv"
    counts.write_text(
        "TIS-other-1\tchr\t3\n"
        "TIS-other-1\tplasmid\t1\n"
    )
    output = tmp_path / "incomplete-batch-hint"

    with pytest.raises(ValueError, match="missing samples: TIS-WT-1"):
        analyze(
            mixed_library, output, normalization="min",
            total_read_file_path=counts,
            _normalization_samples=("TIS-other-1",),
        )
    assert not output.exists()


def test_direct_cli_exports_selected_code_and_library_scope(mixed_library, tmp_path):
    output = tmp_path / "cli"
    command = [sys.executable, str(Path(__file__).resolve().parents[1] / "orfbounder.py"),
               "--alignment_file_tis", str(mixed_library / "TIS-WT-1.sam"),
               "--annotation_file", str(mixed_library / "annotation.gff"), "--genome_file", str(mixed_library / "genome.fa"),
               "--offset_json", str(mixed_library / "offsets.json"), "--mapping_method", "fiveprime",
               "--min_peak_height", "0", "--peak_height_operator", "sum", "--output_path", str(output),
               "--output_basename", "WT-1", "--genetic-code", "4", "--normalization-scope", "library",
               "--normalization_method", "mil", "--min-mapq", "10",
               "--duplicates", "exclude", "--multimappers", "include"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    frame = pd.read_csv(output / "table_per_sample/WT-1.csv", sep="\t")
    assert frame.Amino_Acid_Seq.tolist() == ["MKW*", "MKW*", "MP*"]
    metadata = json.loads((output / "WT-1.run.json").read_text())
    assert metadata["genetic_code"] == 4
    assert metadata["normalization_scope"] == "library"
    assert metadata["alignment_policy"] == {
        "min_mapq": 10, "duplicates": "exclude", "multimappers": "include",
        "count_unit": "alignment_record",
    }
    assert metadata["status"] == "complete"
    assert (output / "WT-1.complete.json").is_file()
    assert (output / "WT-1.report.html").is_file()
    with (output / "sequence_per_sample/WT-1.faa").open(encoding="utf-8") as handle:
        assert [str(record.seq) for record in SeqIO.parse(handle, "fasta")] == [
            "MKW", "MKW", "MP",
        ]
