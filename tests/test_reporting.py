"""Offline result presentation preserves sample evidence and untrusted labels."""

from html.parser import HTMLParser
from io import StringIO
import json
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

import lib.reporting as reporting
from lib.reporting import write_report


class ReportParser(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.scripts = []
        self.links = []
        self.tags = []
        self._script = None
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        self.tags.append((tag, attrs))
        if tag == "script":
            self._script = {"attrs": attrs, "text": ""}
            self.scripts.append(self._script)
        if tag == "a":
            self.links.append(attrs.get("href"))

    def handle_endtag(self, tag):
        if tag == "script":
            self._script = None

    def handle_data(self, data):
        if self._script is not None:
            self._script["text"] += data

    @property
    def data(self):
        return json.loads(next(script["text"] for script in self.scripts
                               if script["attrs"].get("id") == "report-data"))


@pytest.fixture
def candidates():
    return pd.DataFrame({
        "Type": ["Annotated", "Intergenic", "Annotated"],
        "Identifier": ["chr:1-63:+", "chr:100-162:-", "chr:200-262:+"],
        "Genome": ["chr"] * 3,
        "Start": [1, 100, 200], "Stop": [63, 162, 262], "Strand": ["+", "-", "+"],
        "Locus_tag": ["gene1", pd.NA, "gene2"], "Codon_count": [21] * 3,
        "Nucleotide_Seq": ["ATG" + "GCT" * 19 + "TAA"] * 3,
        "Amino_Acid_Seq": ["M" + "A" * 19 + "*"] * 3,
        "15nt_window": ["AAAAAAATGGCTGCTG"] * 3,
        "TIS-WT-1_peak_height": pd.array([8, None, 6], dtype="Float64"),
        "TTS-WT-1_peak_height": pd.array([None, 7, 0], dtype="Float64"),
        "RIBO-WT-1_peak_height": [100, 100, 100],
        "TIS-WT-1_peak_significant": pd.array([False, True, None], dtype="boolean"),
        "TIS-WT-1_peak_q_value": pd.array([0.8, 0.001, None], dtype="Float64"),
        "TIS-WT-1_peak_count": pd.array([8, 4, None], dtype="Int64"),
        "TIS-WT-1_peak_fold_enrichment": pd.array([4.0, np.inf, None], dtype="Float64"),
        "TIS-WT-1_peak_fit_status": pd.array(
            ["ok", "p_value_floored", None], dtype="string",
        ),
        "RIBO-WT-1_peak_significant": [True] * 3,
        "RNA-WT-1_rpkm": [0, np.nan, 10], "WT-1_TE": [np.inf, np.nan, 1.5],
    })


def parse_report(candidates, tmp_path, **kwargs):
    path = write_report(candidates, tmp_path / "output" / "report.html", "Example results", **kwargs)
    return ReportParser(path.read_text(encoding="utf-8"))


def test_report_preserves_all_measurements_and_missing_values(candidates, tmp_path):
    before = candidates.copy(deep=True)
    report = parse_report(candidates, tmp_path)
    assert report.data["columns"] == candidates.columns.tolist()
    frame = report.data
    columns = {name: index for index, name in enumerate(frame["columns"])}
    assert frame["rows"][0]["values"][columns["RNA-WT-1_rpkm"]] == 0
    assert frame["rows"][1]["values"][columns["RNA-WT-1_rpkm"]] is None
    assert frame["rows"][0]["values"][columns["WT-1_TE"]] == "Infinity"
    assert frame["rows"][0]["values"][columns["TIS-WT-1_peak_count"]] == 8
    assert frame["rows"][0]["values"][columns["TIS-WT-1_peak_fold_enrichment"]] == 4
    assert frame["rows"][1]["values"][columns["TIS-WT-1_peak_fold_enrichment"]] == "Infinity"
    assert frame["rows"][2]["values"][columns["TIS-WT-1_peak_fold_enrichment"]] is None
    assert frame["rows"][1]["values"][columns["TIS-WT-1_peak_fit_status"]] == "p_value_floored"
    assert frame["rows"][1]["values"][columns["Locus_tag"]] is None
    pd.testing.assert_frame_equal(candidates, before)


def test_report_preserves_integers_beyond_javascript_number_precision(tmp_path):
    count = (1 << 63) + 1
    frame = pd.DataFrame({
        "Identifier": ["chr:1-3:+"],
        "TIS-WT-1_peak_height": pd.Series([count], dtype=object),
        "TIS-WT-1_peak_count": pd.Series([count], dtype=object),
    })
    report = parse_report(frame, tmp_path, run_metadata=[{
        "alignment_diagnostics": {
            "TIS-WT-1": {"accepted_reads_by_contig": {"chr": count}},
        },
    }])
    columns = {name: index for index, name in enumerate(report.data["columns"])}

    assert report.data["rows"][0]["values"][
        columns["TIS-WT-1_peak_height"]
    ] == str(count)
    assert report.data["rows"][0]["values"][
        columns["TIS-WT-1_peak_count"]
    ] == str(count)
    assert report.data["rows"][0]["called"] == ["TIS-WT-1"]
    assert report.data["runs"][0]["alignment_diagnostics"]["TIS-WT-1"][
        "accepted_reads_by_contig"
    ]["chr"] == str(count)
    assert "BigInt(value).toLocaleString()" in report.scripts[1]["text"]


def test_report_json_integer_precision_boundary_includes_numpy_nullable_scalars():
    safe = (1 << 53) - 1

    assert reporting._json_value(np.int64(safe)) == safe
    assert reporting._json_value(np.int64(-safe)) == -safe
    assert reporting._json_value(np.int64(safe + 1)) == str(safe + 1)
    assert reporting._json_value(np.int64(-(safe + 1))) == str(-(safe + 1))
    assert reporting._json_value(pd.array([safe + 1], dtype="Int64")[0]) == str(
        safe + 1
    )


def test_called_peaks_and_local_support_remain_independent(candidates, tmp_path):
    report = parse_report(candidates, tmp_path)
    rows = report.data["rows"]
    assert [row["called"] for row in rows] == [["TIS-WT-1"], ["TTS-WT-1"], ["TIS-WT-1"]]
    assert [row["significant"] for row in rows] == [[], ["TIS-WT-1"], []]
    assert [row["tested"] for row in rows] == [["TIS-WT-1"], ["TIS-WT-1"], []]
    assert "WT-1" not in [sample["name"] for sample in report.data["samples"]]
    source = (tmp_path / "output/report.html").read_text()
    assert 'With a called TIS peak</span><strong>2</strong>' in source
    assert 'With a called TTS peak</span><strong>1</strong>' in source
    assert 'With local statistical support</span><strong>1</strong>' in source
    assert "not a combined test" in source
    assert "Local fold enrichment" in source
    assert "Fit status" in source
    assert "floored p-value remains positive" in source
    assert "numerical limit means p=1" in source
    assert "correction limit means the valid raw p-value was retained but q=1" in source
    assert "Infinity means a positive peak with zero background endpoints and no pseudocount" in source
    assert "TE uses the unrounded rates" in source
    assert report.data["matched_comparisons"] == []


def test_sample_assay_from_run_record_handles_direct_custom_filenames(tmp_path):
    frame = pd.DataFrame({"start_library_peak_height": [8], "start_library_peak_significant": [True],
                          "stop_peak_height": [7], "TTS_peak_height": [np.nan], "RIBO_peak_height": [np.nan]})
    # Existing output columns take the basename before its first underscore.
    frame = frame.rename(columns={"start_library_peak_height": "start_peak_height",
                                  "start_library_peak_significant": "start_peak_significant"})
    report = parse_report(frame, tmp_path, run_metadata=[{
        "inputs": {"TIS": "/reads/start_library.bam", "TTS": "/reads/stop.bam", "RIBO": None},
        "genetic_code": np.int64(11), "normalization_scope": "library",
    }])
    assert report.data["rows"][0]["called"] == ["start", "stop"]
    assert report.data["rows"][0]["significant"] == ["start"]
    assert report.data["samples"] == [{"name": "start", "assay": "TIS"}, {"name": "stop", "assay": "TTS"}]
    assert report.data["runs"][0]["genetic_code"] == 11


def test_symlink_sample_label_takes_precedence_over_resolved_input_filename(tmp_path):
    report = parse_report(pd.DataFrame({"custom_peak_height": [8]}), tmp_path,
                          run_metadata=[{"inputs": {"TIS": "/alignments/mapped.bam"},
                                         "sample_names": {"TIS": "custom"}}])
    assert report.data["samples"] == [{"name": "custom", "assay": "TIS"}]
    assert report.data["rows"][0]["called"] == ["custom"]


def test_empty_report_is_valid_and_explains_next_steps(candidates, tmp_path):
    report = parse_report(candidates.iloc[:0], tmp_path)
    assert report.data["rows"] == []
    source = (tmp_path / "output/report.html").read_text()
    assert "ORF candidates</span><strong>0</strong>" in source
    assert "No candidate categories to display" in source
    assert "No ORF candidates were reported" in source
    assert "Check the read counts and settings" in source


def test_statistics_not_run_is_not_presented_as_negative_support(tmp_path):
    report = parse_report(pd.DataFrame({"TIS-WT-1_peak_height": [8]}), tmp_path)
    assert report.data["rows"][0]["tested"] == []
    source = (tmp_path / "output/report.html").read_text()
    assert "Local statistics were not included in this table" in source
    assert "Local statistics unavailable" in source


def test_empty_merged_table_uses_run_metadata_to_report_enabled_statistics(tmp_path):
    # Empty batch merges retain only identity fields, while the sample runs
    # still write their complete statistical candidate families.
    report = parse_report(pd.DataFrame(columns=["Type", "Identifier", "Genome", "Start", "Stop"]),
                          tmp_path, run_metadata=[{"statistics": "local", "orf_count": 0}])
    assert report.data["rows"] == []
    source = (tmp_path / "output/report.html").read_text()
    assert "Local statistics were enabled, but this result contains no ORF candidates" in source
    assert "statistical candidate TSV files" in source
    assert "With local statistical support</span><strong>0</strong>" in source
    assert "Local statistics were not included in this table" not in source
    assert "Local statistics unavailable" not in source


def test_text_false_is_not_truthy_statistical_support(tmp_path):
    report = parse_report(pd.DataFrame({"TIS-WT-1_peak_significant": ["False", "True", None]}), tmp_path)
    assert [row["significant"] for row in report.data["rows"]] == [[], ["TIS-WT-1"], []]


def test_matched_comparison_metadata_is_deduplicated_and_separate_from_samples(tmp_path):
    prefix = "treated_vs_control_TIS"
    frame = pd.DataFrame({
        "Identifier": ["chr:1-30:+", "chr:40-69:+", "chr:80-109:+"],
        "TIS-treated-1_peak_height": [8, 6, 4],
        f"{prefix}_matched_pair_count": pd.array([2, 2, 2], dtype="Int64"),
        f"{prefix}_informative_pair_count": pd.array([2, 1, 0], dtype="Int64"),
        f"{prefix}_concordant_pair_count": pd.array([2, 0, 0], dtype="Int64"),
        f"{prefix}_discordant_pair_count": pd.array([0, 1, 0], dtype="Int64"),
        f"{prefix}_tied_pair_count": pd.array([0, 0, 0], dtype="Int64"),
        f"{prefix}_numerator_peak_count": pd.array([20, 4, 0], dtype="Int64"),
        f"{prefix}_numerator_background_count": pd.array([3, 8, 0], dtype="Int64"),
        f"{prefix}_denominator_peak_count": pd.array([2, 5, 0], dtype="Int64"),
        f"{prefix}_denominator_background_count": pd.array([18, 7, 0], dtype="Int64"),
        f"{prefix}_common_odds_ratio": pd.array([30, 0.7, None], dtype="Float64"),
        f"{prefix}_log2_common_odds_ratio": pd.array([4.9069, -0.5146, None], dtype="Float64"),
        f"{prefix}_p_value": pd.array([0.0001, None, None], dtype="Float64"),
        f"{prefix}_q_value": pd.array([0.001, None, None], dtype="Float64"),
        f"{prefix}_significant": pd.array([True, None, None], dtype="boolean"),
        f"{prefix}_correction": pd.array(["by", "by", "by"], dtype="string"),
        f"{prefix}_tests_in_family": pd.array([2, 2, 2], dtype="Int64"),
        f"{prefix}_fit_status": pd.array(
            ["ok", "insufficient_informative_pairs", "no_informative_pairs"], dtype="string"
        ),
    })
    comparison = {
        "comparison": "treated_vs_control", "assay": "TIS",
        "numerator_condition": "treated", "denominator_condition": "control",
        "fdr": 0.05, "fdr_method": "by",
        "model": "exact_stratified_fixed_margin_endpoint_count",
        "numeric_backend": "extended_longdouble",
        "alternative": "greater", "replicate_ids": ["1", "2"],
        "minimum_informative_pairs": 2,
        "support_rule": "q_value<=fdr and common_odds_ratio>1 and strict majority",
        "max_exact_convolution_work": 2_000_000,
        "maximum_supported_hypergeometric_total": (1 << 31) - 1,
        "minimum_positive_p_value": 5e-324,
        "candidate_table": "matched_statistics/treated_vs_control_TIS.tsv",
        "strata_table": "matched_statistics/treated_vs_control_TIS_strata.tsv",
    }
    report = parse_report(frame, tmp_path, run_metadata=[
        {"output_basename": "treated-1", "matched_comparisons": [comparison]},
        {"output_basename": "treated-2", "matched_comparisons": [dict(comparison)]},
    ])

    assert report.data["matched_comparisons"] == [{"prefix": prefix, **comparison}]
    assert report.data["samples"] == [{"name": "TIS-treated-1", "assay": "TIS"}]
    assert [row["matched_tested"] for row in report.data["rows"]] == [[prefix], [], []]
    assert [row["matched_significant"] for row in report.data["rows"]] == [[prefix], [], []]
    source = (tmp_path / "output/report.html").read_text()
    assert 'With exploratory fixed-margin matched support</span><strong>1</strong>' in source
    assert "Exploratory matched support" in source
    assert "Exploratory fixed-margin matched evidence" in source
    for label in ("Conditions", "Matched pairs", "Informative pairs", "Numerator-favoring pairs",
                  "Denominator-favoring pairs", "Tied pairs", "Common OR", "log2 OR", "q-value",
                  "Fit status"):
        assert label in source
    assert "at least two informative pairs" in source
    assert "strict majority" in source
    assert "q-value at or below the configured FDR" in source
    assert "common odds ratio above 1" in source
    assert "does not estimate biological dispersion or a population-level condition effect" in source
    assert "predeclared independent experimental blocks" in source
    assert "Separate comparisons and assays are corrected separately" in source
    assert "companion _strata.tsv table, joined by candidate_id" in source


def test_conflicting_repeated_matched_metadata_is_rejected(tmp_path):
    prefix = "treated_vs_control_TIS"
    frame = pd.DataFrame({f"{prefix}_matched_pair_count": pd.array([2], dtype="Int64")})
    comparison = {
        "comparison": "treated_vs_control", "assay": "TIS",
        "numerator_condition": "treated", "denominator_condition": "control",
        "model": "exact_stratified_fixed_margin_endpoint_count",
    }
    conflicting = {**comparison, "denominator_condition": "unrelated"}
    with pytest.raises(ValueError, match="Conflicting matched-comparison metadata"):
        parse_report(frame, tmp_path, run_metadata=[
            {"matched_comparisons": [comparison]},
            {"matched_comparisons": [conflicting]},
        ])


def test_matched_comparison_is_safely_inferred_from_columns_without_run_metadata(tmp_path):
    prefix = "drug_vs_vehicle_TTS"
    frame = pd.DataFrame({
        f"{prefix}_matched_pair_count": pd.array([3], dtype="Int64"),
        f"{prefix}_numerator_peak_count": pd.array([9], dtype="Int64"),
        f"{prefix}_q_value": pd.array([0.4], dtype="Float64"),
        f"{prefix}_significant": pd.array([False], dtype="boolean"),
        f"{prefix}_fit_status": pd.array(["ok"], dtype="string"),
    })
    report = parse_report(frame, tmp_path)
    assert report.data["matched_comparisons"] == [{
        "prefix": prefix, "comparison": "drug_vs_vehicle", "assay": "TTS",
        "numerator_condition": None, "denominator_condition": None,
    }]
    assert report.data["samples"] == []
    assert report.data["rows"][0]["matched_tested"] == [prefix]
    assert report.data["rows"][0]["matched_significant"] == []
    source = (tmp_path / "output/report.html").read_text()
    assert 'With exploratory fixed-margin matched support</span><strong>0</strong>' in source
    assert "Not recorded" in source


def test_labels_metadata_sequences_and_template_markers_cannot_inject_html(tmp_path):
    hostile = '</script><script>alert("injected")</script><img src=x onerror=alert(1)> & __DATA__ \u2028'
    frame = pd.DataFrame({"Type": [hostile], "Locus_tag": [hostile], "Nucleotide_Seq": [hostile]})
    path = write_report(frame, tmp_path / "report.html", hostile, run_metadata=[{"output_basename": hostile}],
                        links={hostile: 'table/space%20and%20%22quote%22.tsv'})
    source = path.read_text()
    report = ReportParser(source)
    assert len(report.scripts) == 2
    assert not any(tag == "img" for tag, _ in report.tags)
    assert not any(key.startswith("on") for _, attrs in report.tags for key in attrs)
    assert report.data["rows"][0]["values"] == [hostile] * 3
    assert report.data["runs"][0]["output_basename"] == hostile
    assert "&lt;/script&gt;" in source
    assert "__DATA__" in source
    assert 'table/space%20and%20%22quote%22.tsv' in report.links
    assert "innerHTML" not in report.scripts[1]["text"]


@pytest.mark.parametrize("target", ["javascript:alert(1)", "data:text/html,hi", "https://example.com/table.tsv",
                                  "//example.com/table.tsv", "/absolute.tsv", "  /absolute.tsv", "\\\\host\\file", "table\n.tsv", ""])
def test_only_relative_result_links_are_accepted(tmp_path, target):
    with pytest.raises(ValueError, match="relative paths"):
        write_report(pd.DataFrame(), tmp_path / "report.html", "Results", links={"Table": target})
    assert not (tmp_path / "report.html").exists()


def test_report_is_self_contained_and_links_local_exports(candidates, tmp_path):
    report = parse_report(candidates, tmp_path, links={"TSV": "tables/result.csv", "GFF": "../features.gff"})
    assert {"tables/result.csv", "../features.gff"}.issubset(report.links)
    assert not any("src" in script["attrs"] for script in report.scripts)
    assert not any(tag == "link" and attrs.get("rel") == "stylesheet" for tag, attrs in report.tags)
    assert any(tag == "label" for tag, _ in report.tags)
    assert any(attrs.get("aria-live") == "polite" for _, attrs in report.tags)
    assert 'setAttribute("aria-expanded"' in report.scripts[1]["text"]


def test_report_streams_json_and_does_not_build_an_eager_search_copy(
        candidates, tmp_path, monkeypatch):
    def reject_full_json_copy(*args, **kwargs):
        raise AssertionError("write_report must stream JSON instead of using json.dumps")

    def reject_eager_table_copy(*args, **kwargs):
        raise AssertionError("write_report must not materialize every prepared report row")

    monkeypatch.setattr(reporting.json, "dumps", reject_full_json_copy)
    monkeypatch.setattr(reporting, "_prepare_data", reject_eager_table_copy)
    report = parse_report(candidates, tmp_path)

    assert len(report.data["rows"]) == len(candidates)
    browser_script = report.scripts[1]["text"]
    assert "const searchText" not in browser_script
    assert "row.values.some" in browser_script


def test_incremental_report_json_matches_eager_representation(candidates):
    metadata = [{"genetic_code": np.int64(11), "output_basename": "sample</script>"}]
    schema = reporting._prepare_schema(candidates, metadata)
    expected = StringIO()
    reporting._write_script_safe_json(reporting._prepare_data(candidates, metadata), expected)
    incremental = StringIO()
    reporting._write_report_data(candidates, schema, incremental)

    assert incremental.getvalue() == expected.getvalue()


def test_duplicate_column_names_fail_before_writing(tmp_path):
    with pytest.raises(ValueError, match="unique column"):
        write_report(pd.DataFrame([[1, 2]], columns=["Type", "Type"]), tmp_path / "report.html", "Results")


@pytest.mark.skipif(not shutil.which("node"), reason="Node is optional for checking the embedded browser script")
def test_embedded_browser_javascript_is_valid(candidates, tmp_path):
    report = parse_report(candidates, tmp_path)
    script = tmp_path / "report.js"
    script.write_text(report.scripts[1]["text"], encoding="utf-8")
    result = subprocess.run([shutil.which("node"), "--check", str(script)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
