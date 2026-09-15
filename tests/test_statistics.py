"""Statistical calibration, testing-family and coordinate regression checks."""

from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
from math import comb

import numpy as np
import pandas as pd
import pytest

import lib.statistics as local_statistics
from lib.statistics import (
    MAX_SUPPORTED_TOTAL_COUNT,
    MIN_POSITIVE_P_VALUE,
    _binomial_upper_tail,
    analyze_peak_enrichment,
    attach_to_orfs,
)


SEQUENCE = "AAAAAATGAAAAAAA"  # One + ATG, length 15, peak at 0-based 3..7.


@pytest.mark.parametrize("successes", range(9))
def test_exact_tail_matches_enumerated_conditional_null(successes):
    counts = {("chr", "+"): {5: successes, 0: 8 - successes}}
    result = analyze_peak_enrichment({"chr": SEQUENCE}, counts, ["ATG"], flank_width=50)
    expected = sum(comb(8, x) * (1 / 3) ** x * (2 / 3) ** (8 - x) for x in range(successes, 9))
    assert result.iloc[0].p_value == pytest.approx(expected)
    assert result.iloc[0].peak_count == successes
    assert result.iloc[0].background_count == 8 - successes


def test_direct_tail_guard_prevents_half_probability_threshold_crossing():
    value, status = _binomial_upper_tail(8, 15, 0.5)

    assert status == "ok"
    assert Fraction.from_float(value) >= Fraction(1, 2)

    # One five-base window and five background bases give the exact null p=1/2.
    # scipy's unguarded sf is one ULP below 1/2 for 8 successes in 15 trials.
    threshold = np.nextafter(0.5, 0.0)
    row = analyze_peak_enrichment(
        {"chr": "AAAAAATGAA"},
        {("chr", "+"): {5: 8, 0: 7}},
        ["ATG"], flank_width=50, correction="bh", fdr=threshold,
    ).iloc[0]

    assert row.peak_width == row.background_width == 5
    assert row.q_value > threshold
    assert not row.significant


def test_direct_tail_guard_does_not_undershoot_bounded_exact_width_fractions():
    # Peak windows are three to five bases after linear-reference clipping.
    # Sweep representative background geometries, counts and every possible
    # observed success count against an independent exact-rational sum.
    for peak_width in range(3, 6):
        for background_width in range(1, 13):
            width_total = peak_width + background_width
            probability = peak_width / width_total
            for trials in range(1, 31):
                denominator = width_total ** trials
                tail_numerator = 0
                for successes in range(trials, -1, -1):
                    tail_numerator += (
                        comb(trials, successes)
                        * peak_width ** successes
                        * background_width ** (trials - successes)
                    )
                    exact = Fraction(tail_numerator, denominator)
                    value, status = _binomial_upper_tail(
                        successes, trials, probability,
                    )
                    assert status == "ok"
                    assert Fraction.from_float(value) >= exact, (
                        peak_width, background_width, trials, successes,
                    )


def test_family_includes_uncovered_contigs_and_strands_and_by_is_conservative():
    genome = {"covered": SEQUENCE, "uncovered": SEQUENCE, "reverse": "AAAAACATAAAAAAA"}
    counts = {("covered", "+"): {5: 3}, ("covered", "-"): {5: 1000}}
    by = analyze_peak_enrichment(genome, counts, ["ATG"], correction="by")
    bh = analyze_peak_enrichment(genome, counts, ["ATG"], correction="bh")
    assert len(by) == 3
    assert set(by.tests_in_family) == {3}
    covered = by.loc[by.Genome == "covered"].iloc[0]
    assert covered.peak_count == 3
    assert covered.background_count == 0
    assert covered.p_value == pytest.approx(1 / 27)
    assert bh.loc[bh.Genome == "covered", "q_value"].iloc[0] == pytest.approx(3 / 27)
    assert covered.q_value == pytest.approx((3 / 27) * (1 + 1 / 2 + 1 / 3))
    assert (by.loc[by.Genome != "covered", "p_value"] == 1).all()


def test_premature_scipy_underflow_is_recomputed_in_log_space():
    # A centered five-base window with two 100-base flanks has p0=1/41.
    # scipy 1.15 reports zero for this tail even though its exact value is
    # 1.102832261181031756...e-288.
    sequence = "A" * 102 + "ATG" + "A" * 100
    counts = {("chr", "+"): {102: 206, 0: 38}}
    row = analyze_peak_enrichment(
        {"chr": sequence}, counts, ["ATG"], flank_width=100,
    ).iloc[0]

    exact = Fraction(
        sum(comb(244, x) * 40 ** (244 - x) for x in range(206, 245)),
        41 ** 244,
    )
    assert row.p_value == pytest.approx(float(exact), rel=1e-12)
    assert row.q_value == pytest.approx(row.p_value)
    assert row.p_value > 0
    assert row.fit_status == "ok"


def test_representable_subnormal_tail_is_rounded_up_not_blindly_floored():
    # Independent Decimal recurrence gives 6.29758547940374686...e-324 for
    # p0=1/21. It lies between float64's first and second positive values.
    value, status = _binomial_upper_tail(7564, 100_000, 1 / 21)
    reference = Decimal("6.29758547940374686e-324")

    assert status == "ok"
    assert value == np.nextafter(MIN_POSITIVE_P_VALUE, np.inf)
    assert Decimal.from_float(value) >= reference


def test_subnormal_fallback_does_not_undershoot_exact_decimal_tail():
    # An independent 100-digit Decimal binomial recurrence gives this tail.
    # scipy.sf is subnormal here, so the log fallback is used; float64 logpmf
    # previously undershot the true tail by about 2.4e-10 relative.
    reference = Decimal("2.3542044042689476614790010445984618687081926483419685e-310")
    value, status = _binomial_upper_tail(7_500, 100_000, 1 / 21)

    assert status == "ok"
    assert Decimal.from_float(value) >= reference
    assert Decimal.from_float(value) / reference < Decimal("1.00000000001")


def test_extreme_tail_combination_work_limit_fails_closed(monkeypatch):
    monkeypatch.setattr(local_statistics.binom, "sf", lambda *args: 0.0)
    value, status = _binomial_upper_tail(100_001, 300_000, 1 / 21)

    assert value == 1.0
    assert status == "numerical_limit"


def test_underflow_fallback_bounds_exact_window_width_ratio(monkeypatch):
    observed_probabilities = []
    original = local_statistics._log_binomial_upper_tail

    def capture_probability(successes, trials, probability):
        observed_probabilities.append(probability)
        return original(successes, trials, probability)

    monkeypatch.setattr(local_statistics.binom, "sf", lambda *args: 0.0)
    monkeypatch.setattr(local_statistics, "_log_binomial_upper_tail", capture_probability)
    row = analyze_peak_enrichment(
        {"chr": SEQUENCE}, {("chr", "+"): {5: 8}}, ["ATG"], flank_width=50,
    ).iloc[0]

    assert row.fit_status == "ok"
    assert observed_probabilities
    assert Fraction.from_float(observed_probabilities[0]) >= Fraction(1, 3)
    assert observed_probabilities[0] == np.nextafter(1 / 3, np.inf)


def test_log_tail_recurrence_matches_enumeration_when_direct_backend_fails(monkeypatch):
    monkeypatch.setattr(local_statistics.binom, "sf", lambda *args, **kwargs: 0.0)
    value, status = _binomial_upper_tail(12, 20, 1 / 3)
    expected = sum(comb(20, x) * (1 / 3) ** x * (2 / 3) ** (20 - x)
                   for x in range(12, 21))

    assert status == "ok"
    assert value == pytest.approx(expected, rel=1e-13)


def test_direct_binomial_backend_exception_uses_safe_log_tail(monkeypatch):
    def fail_sf(*args):
        raise FloatingPointError("simulated direct backend failure")

    monkeypatch.setattr(local_statistics.binom, "sf", fail_sf)
    value, status = _binomial_upper_tail(10, 10, 1 / 21)

    assert status == "ok"
    assert value == pytest.approx((1 / 21) ** 10, rel=1e-13)


def test_true_binomial_underflow_is_floored_and_never_exports_zero():
    # Here p=(1/21)^245 is positive but below the float64 output range.
    sequence = "A" * 52 + "ATG" + "A" * 50
    row = analyze_peak_enrichment(
        {"chr": sequence}, {("chr", "+"): {52: 245}}, ["ATG"], flank_width=50,
    ).iloc[0]

    assert row.p_value == MIN_POSITIVE_P_VALUE
    assert row.q_value == MIN_POSITIVE_P_VALUE
    assert row.fit_status == "p_value_floored"
    assert row.significant


def test_counts_beyond_exact_float_range_stay_in_family_conservatively():
    sequence = "A" * 52 + "ATG" + "A" * 50
    row = analyze_peak_enrichment(
        {"chr": sequence},
        {("chr", "+"): {52: MAX_SUPPORTED_TOTAL_COUNT + 1}},
        ["ATG"],
        flank_width=50,
        fdr=1,
    ).iloc[0]

    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert row.fit_status == "numerical_limit"
    assert not row.significant


def test_local_backend_exception_stays_in_family_conservatively(monkeypatch):
    sequence = "A" * 52 + "ATG" + "A" * 50
    monkeypatch.setattr(local_statistics.binom, "sf", lambda *args: 0.0)

    def fail_combination(*args):
        raise FloatingPointError("simulated backend failure")

    monkeypatch.setattr(local_statistics.math, "comb", fail_combination)
    row = analyze_peak_enrichment(
        {"chr": sequence}, {("chr", "+"): {52: 10}}, ["ATG"], fdr=1,
    ).iloc[0]

    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert row.fit_status == "numerical_limit"
    assert not row.significant


@pytest.mark.parametrize("failure", ["nonfinite", "exception"])
def test_local_correction_backend_failure_preserves_raw_p_and_fails_closed(
    monkeypatch, failure,
):
    sequence = "A" * 52 + "ATG" + "A" * 50
    if failure == "nonfinite":
        def backend(p_values, method):
            return np.full(np.shape(p_values), np.nan)
    else:
        def backend(*args, **kwargs):
            raise FloatingPointError("simulated correction backend failure")

    monkeypatch.setattr(local_statistics, "false_discovery_control", backend)
    row = analyze_peak_enrichment(
        {"chr": sequence}, {("chr", "+"): {52: 10}}, ["ATG"], fdr=1,
    ).iloc[0]

    assert 0 < row.p_value < 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert row.fit_status == "correction_limit"
    assert not row.significant


def test_local_materially_underadjusted_family_fails_closed(monkeypatch):
    sequence = "A" * 52 + "ATG" + "A" * 20 + "ATG" + "A" * 50
    monkeypatch.setattr(
        local_statistics, "false_discovery_control",
        lambda p_values, method: np.zeros(np.shape(p_values)),
    )
    rows = analyze_peak_enrichment(
        {"chr": sequence}, {("chr", "+"): {52: 10, 75: 10}}, ["ATG"], fdr=1,
    )

    assert len(rows) == 2
    assert rows.p_value.between(0, 1, inclusive="neither").all()
    assert rows.q_value.eq(1).all()
    assert rows.fit_status.eq("correction_limit").all()
    assert not rows.significant.any()


def test_bh_adjustment_does_not_round_below_exact_input_float_result():
    p_values = np.array([0.001, 0.005, 0.01])
    q_values = local_statistics._adjust_p_values(p_values, "bh")
    exact_middle = Fraction.from_float(float(p_values[1])) * 3 / 2

    assert Fraction.from_float(float(q_values[1])) >= exact_middle


def test_upward_adjustment_guard_prevents_threshold_support_edge(monkeypatch):
    p_values = np.array([0.001, 0.005, 0.01])
    unguarded_threshold = float(
        local_statistics.false_discovery_control(p_values, method="bh")[1]
    )
    generated = iter(p_values)
    monkeypatch.setattr(
        local_statistics, "_binomial_upper_tail",
        lambda *args: (float(next(generated)), "ok"),
    )
    rows = analyze_peak_enrichment(
        {"chr": "ATGAAAATGAAAATG"},
        {("chr", "+"): {1: 1, 6: 1, 12: 1}},
        ["ATG"], flank_width=1, correction="bh", fdr=unguarded_threshold,
    )

    assert rows.loc[1, "q_value"] > unguarded_threshold
    assert not rows.loc[1, "significant"]


def test_oversized_local_count_remains_exact_when_attached_to_orf(tmp_path):
    count = 1 << 63
    statistics = analyze_peak_enrichment(
        {"chr": SEQUENCE}, {("chr", "+"): {5: count}}, ["ATG"],
    )
    boundary = int(statistics.iloc[0].Codon_start)
    table = pd.DataFrame([{
        "Genome": "chr", "Start": boundary, "Stop": 15, "Strand": "+",
    }])

    result = attach_to_orfs(table, statistics, "TIS", "sample")

    assert statistics.iloc[0].fit_status == "numerical_limit"
    assert result.loc[0, "sample_peak_count"] == count
    assert result.sample_peak_count.dtype == object
    assert str(result.sample_peak_tests_in_family.dtype) == "Int64"
    output = tmp_path / "attached.tsv"
    result.to_csv(output, sep="\t", index=False)
    assert str(count) in output.read_text()


def test_raw_counts_are_preserved_and_background_excludes_the_peak():
    counts = {("chr", "+"): {0: 2, 3: 3, 5: 4, 7: 5, 8: 6, 14: 7, 15: 200}}
    original = deepcopy(counts)
    row = analyze_peak_enrichment({"chr": SEQUENCE}, counts, ["ATG"], flank_width=2).iloc[0]
    assert row.peak_count == 12
    assert row.background_count == 6
    assert row.background_width == 4
    assert row.fold_enrichment == pytest.approx((12 / 5) / (6 / 4))
    assert counts == original


@pytest.mark.parametrize("count", [1.5, 3.0, -1, np.nan, True])
def test_normalized_or_invalid_counts_are_rejected(count):
    with pytest.raises(ValueError, match="raw integer counts"):
        analyze_peak_enrichment({"chr": SEQUENCE}, {("chr", "+"): {5: count}}, ["ATG"])


@pytest.mark.parametrize("kwargs", [{"flank_width": 0}, {"flank_width": 2.5}, {"fdr": 0}, {"fdr": 1.1}, {"fdr": np.nan}, {"fdr": True}, {"fdr": "0.05"}, {"correction": "wrong"}])
def test_invalid_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"], **kwargs)


def test_boundary_windows_clip_and_use_actual_widths():
    genome = {"left": "ATGAAAAAAA", "right": "AAAAAAACAT"}
    counts = {("left", "+"): {0: 5, 3: 2}, ("right", "-"): {9: 5, 6: 2}}
    result = analyze_peak_enrichment(genome, counts, ["ATG"], flank_width=2)
    assert list(result.peak_width) == [3, 3]
    assert list(result.background_width) == [2, 2]
    assert list(result.peak_count) == [5, 5]
    assert list(result.background_count) == [2, 2]
    assert list(result.Window_start) == [1, 8]
    assert list(result.Codon_start) == [1, 8]
    assert result.p_value.iloc[0] == pytest.approx(result.p_value.iloc[1])


def test_exactly_neutral_fold_enrichment_cannot_round_into_support():
    row = analyze_peak_enrichment(
        {"chr": "ATG" + "A" * 17},
        {("chr", "+"): {1: 3, 3: 17}},
        ["ATG"], flank_width=50, correction="bh", fdr=0.9,
    ).iloc[0]

    assert (row.peak_count, row.background_count) == (3, 17)
    assert (row.peak_width, row.background_width) == (3, 17)
    assert row.q_value < 0.9
    assert row.fold_enrichment == 1
    assert not row.significant


def test_no_background_or_no_reads_never_reports_support():
    result = analyze_peak_enrichment({"tiny": "ATG", "empty": SEQUENCE}, {("tiny", "+"): {0: 100}}, ["ATG"])
    assert (result.p_value == 1).all()
    assert (result.q_value == 1).all()
    assert not result.significant.any()
    assert result.fold_enrichment.isna().all()


def test_no_codons_returns_stable_empty_table():
    result = analyze_peak_enrichment({"chr": "AAAAAA"}, {}, ["ATG"])
    assert result.empty
    assert "q_value" in result.columns
    assert "Codon_start" in result.columns


@pytest.mark.parametrize("method, strand, codon_start, boundary", [
    ("TIS", "+", 1, 1), ("TIS", "-", 18, 20),
    ("TTS", "+", 18, 20), ("TTS", "-", 1, 1),
])
def test_statistics_attach_to_correct_orf_boundary(method, strand, codon_start, boundary):
    table = pd.DataFrame([{"Genome": "chr", "Start": 1, "Stop": 20, "Strand": strand}])
    statistics = pd.DataFrame([{
        "Genome": "chr", "Strand": strand, "Codon_start": codon_start,
        "peak_count": 20, "background_count": 1, "fold_enrichment": 8.0,
        "p_value": 0.001,
        "q_value": 0.02, "significant": True, "correction": "by", "tests_in_family": 100,
        "fit_status": "p_value_floored",
    }])
    result = attach_to_orfs(table, statistics, method, "sample")
    assert result.loc[0, "sample_peak_count"] == 20
    assert result.loc[0, "sample_peak_fold_enrichment"] == 8
    assert result.loc[0, "sample_peak_q_value"] == 0.02
    assert result.loc[0, "sample_peak_correction"] == "by"
    assert result.loc[0, "sample_peak_fit_status"] == "p_value_floored"
    assert list(table.columns) == ["Genome", "Start", "Stop", "Strand"]


def test_missing_boundary_evidence_stays_missing_and_does_not_drop_orfs():
    statistics = analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"])
    table = pd.DataFrame([{"Genome": "chr", "Start": 2, "Stop": 15, "Strand": "+"}])
    result = attach_to_orfs(table, statistics, "TIS", "sample")
    assert len(result) == 1
    assert pd.isna(result.loc[0, "sample_peak_fold_enrichment"])
    assert pd.isna(result.loc[0, "sample_peak_q_value"])
    assert pd.isna(result.loc[0, "sample_peak_significant"])


def test_duplicate_local_statistical_boundaries_cannot_overwrite_evidence():
    statistics = analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"])
    duplicate = statistics.copy()
    duplicate.loc[0, "peak_count"] = 100
    statistics = pd.concat([statistics, duplicate], ignore_index=True)
    table = pd.DataFrame([{
        "Genome": "chr", "Start": int(statistics.iloc[0].Codon_start),
        "Stop": len(SEQUENCE), "Strand": "+",
    }])

    with pytest.raises(ValueError, match="Duplicate local-statistics boundary"):
        attach_to_orfs(table, statistics, "TIS", "sample")


def test_attaching_local_family_cannot_overwrite_existing_evidence():
    statistics = analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"])
    table = pd.DataFrame([{
        "Genome": "chr", "Start": int(statistics.iloc[0].Codon_start),
        "Stop": len(SEQUENCE), "Strand": "+", "sample_peak_q_value": 0.01,
    }])
    original = table.copy(deep=True)

    with pytest.raises(ValueError, match="output columns already exist.*sample_peak_q_value"):
        attach_to_orfs(table, statistics, "TIS", "sample")

    pd.testing.assert_frame_equal(table, original)


def test_empty_orf_table_keeps_statistical_schema():
    statistics = analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"])
    table = pd.DataFrame(columns=["Genome", "Start", "Stop", "Strand"])
    result = attach_to_orfs(table, statistics, "TIS", "sample")
    assert result.empty
    assert "sample_peak_q_value" in result.columns
    assert str(result["sample_peak_count"].dtype) == "Int64"
    assert str(result["sample_peak_fold_enrichment"].dtype) == "Float64"
    assert str(result["sample_peak_significant"].dtype) == "boolean"
    assert str(result["sample_peak_fit_status"].dtype) == "string"


def test_attached_fold_enrichment_preserves_zero_infinity_and_missing():
    boundary = int(analyze_peak_enrichment({"chr": SEQUENCE}, {}, ["ATG"]).iloc[0].Codon_start)
    table = pd.DataFrame([{"Genome": "chr", "Start": boundary, "Stop": 15, "Strand": "+"}])
    result = table
    for sample, counts in (
        ("peak_only", {("chr", "+"): {5: 5}}),
        ("background_only", {("chr", "+"): {0: 5}}),
        ("uncovered", {}),
    ):
        statistics = analyze_peak_enrichment({"chr": SEQUENCE}, counts, ["ATG"])
        result = attach_to_orfs(result, statistics, "TIS", sample)

    assert result.loc[0, "peak_only_peak_fold_enrichment"] == np.inf
    assert result.loc[0, "background_only_peak_fold_enrichment"] == 0
    assert pd.isna(result.loc[0, "uncovered_peak_fold_enrichment"])
