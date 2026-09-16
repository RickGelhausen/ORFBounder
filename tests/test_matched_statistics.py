"""Exact matched-condition statistics and ORF attachment regression tests."""

from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import comb, log2

import numpy as np
import pandas as pd
import pytest

import lib.matched_statistics as matched_statistics
from lib.matched_statistics import (
    ATTACH_FIELDS, MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL,
    MIN_POSITIVE_P_VALUE, MatchedComparison,
    _exact_conditional_tail,
    analyze_matched_condition,
    attach_matched_to_orfs,
    read_matched_comparisons,
)


def candidate_frame(rows):
    records = []
    for values in rows:
        candidate_id, peak_count, background_count, *optional = values
        overrides = optional[0] if optional else {}
        record = {
            "candidate_id": candidate_id,
            "Genome": "chr",
            "Window_start": 8,
            "Window_stop": 12,
            "Strand": "+",
            "Codon": "ATG",
            "Codon_start": 10,
            "peak_count": peak_count,
            "background_count": background_count,
            "peak_width": 5,
            "background_width": 100,
        }
        record.update(overrides)
        records.append(record)
    return pd.DataFrame.from_records(records)


def matched_inputs(tables):
    return {
        ("TIS", condition, replicate): candidate_frame(rows)
        for condition, replicate, rows in tables
    }


def replace_second_control_table(data, replacement):
    key = ("TIS", "control", "2")
    data[key] = replacement(data[key])


def test_comparison_file_is_explicit_directional_and_validated(tmp_path):
    path = tmp_path / "comparisons.tsv"
    path.write_text(
        "comparison\tassay\tnumerator_condition\tdenominator_condition\n"
        "drug_vs_vehicle\ttis\tdrug\tvehicle\n"
    )
    assert read_matched_comparisons(path) == [
        MatchedComparison("drug_vs_vehicle", "TIS", "drug", "vehicle")
    ]

    path.write_text(
        "comparison\tassay\tnumerator_condition\tdenominator_condition\n"
        "bad name\tRNA\tdrug-a\tdrug-a\n"
    )
    with pytest.raises(ValueError, match="comparison name"):
        read_matched_comparisons(path)


def test_comparison_name_can_be_reused_across_assays_but_not_within_one(tmp_path):
    path = tmp_path / "comparisons.tsv"
    header = "comparison\tassay\tnumerator_condition\tdenominator_condition\n"
    path.write_text(
        header
        + "drug_vs_vehicle\tTIS\tdrug\tvehicle\n"
        + "drug_vs_vehicle\tTTS\tdrug\tvehicle\n"
    )
    assert [comparison.prefix for comparison in read_matched_comparisons(path)] == [
        "drug_vs_vehicle_TIS", "drug_vs_vehicle_TTS",
    ]

    path.write_text(
        header
        + "drug_vs_vehicle\tTIS\tdrug\tvehicle\n"
        + "drug_vs_vehicle\ttis\tdrug\tvehicle\n"
    )
    with pytest.raises(ValueError, match="output prefix must be unique"):
        read_matched_comparisons(path)


def test_exact_tail_is_convolution_of_fixed_margin_strata():
    # Stratum 1 PMF is [1,16,36,16,1]/70. Stratum 2 is [1,3,1]/5.
    # The probability that the observed sum, 3 + 2 = 5, is met or exceeded is
    # 2/35. The Mantel-Haenszel common odds ratio is 17.
    inputs = matched_inputs([
        ("treated", "1", [("site", 3, 1)]),
        ("control", "1", [("site", 1, 3)]),
        ("treated", "2", [("site", 2, 1)]),
        ("control", "2", [("site", 0, 3)]),
    ])
    row = analyze_matched_condition(inputs, "TIS", "treated", "control", correction="bh").iloc[0]
    assert row.p_value == pytest.approx(2 / 35)
    assert row.q_value == pytest.approx(2 / 35)
    assert row.common_odds_ratio == pytest.approx(17)
    assert row.log2_common_odds_ratio == pytest.approx(log2(17))
    assert row.matched_pair_count == 2
    assert row.informative_pair_count == 2
    assert row.concordant_pair_count == 2
    assert row.discordant_pair_count == 0
    assert row.tied_pair_count == 0
    assert row.fit_status == "ok"
    assert not row.significant


def test_three_unequal_strata_match_brute_force_fixed_margin_enumeration():
    tables = [
        (2, 3, 1, 4),
        (1, 2, 2, 1),
        (3, 1, 1, 5),
    ]
    inputs = matched_inputs([
        (condition, str(index), [("site", peak, background)])
        for index, (a, b, c, d) in enumerate(tables, 1)
        for condition, peak, background in (
            ("treated", a, b), ("control", c, d),
        )
    ])

    distribution = {0: Fraction(1)}
    for a, b, c, d in tables:
        numerator_total = a + b
        peak_total = a + c
        total = a + b + c + d
        lower = max(0, numerator_total - (total - peak_total))
        upper = min(numerator_total, peak_total)
        denominator = comb(total, numerator_total)
        stratum = {
            x: Fraction(comb(peak_total, x) * comb(total - peak_total, numerator_total - x), denominator)
            for x in range(lower, upper + 1)
        }
        convolved = defaultdict(Fraction)
        for left, left_probability in distribution.items():
            for right, right_probability in stratum.items():
                convolved[left + right] += left_probability * right_probability
        distribution = dict(convolved)
    expected = sum(
        probability for value, probability in distribution.items()
        if value >= sum(table[0] for table in tables)
    )

    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
    ).iloc[0]
    assert row.p_value == pytest.approx(float(expected), rel=1e-14)


def test_mode_centered_recurrence_corrects_logpmf_cancellation():
    strata = ((5, 6, 0, 6), (5, 6, 0, 7))
    distributions = []
    for a, b, c, d in strata:
        total = a + b + c + d
        peak_total = a + c
        numerator_total = a + b
        distributions.append({
            x: Fraction(
                comb(peak_total, x)
                * comb(total - peak_total, numerator_total - x),
                comb(total, numerator_total),
            )
            for x in range(
                max(0, numerator_total - (total - peak_total)),
                min(numerator_total, peak_total) + 1,
            )
        })
    exact = sum(
        (left_probability * right_probability
         for left, left_probability in distributions[0].items()
         for right, right_probability in distributions[1].items()
         if left + right >= 10),
        Fraction(0),
    )
    computed, floored = _exact_conditional_tail(strata)

    assert not floored
    assert Fraction.from_float(computed) >= exact
    assert computed == pytest.approx(float(exact), rel=1e-13)


def test_all_small_informative_two_stratum_tails_match_exact_fractions():
    strata = [
        table for table in product(range(3), repeat=4)
        if table[0] + table[1] > 0 and table[2] + table[3] > 0
        and table[0] + table[2] > 0 and table[1] + table[3] > 0
    ]
    distributions = {}
    for table in strata:
        a, b, c, d = table
        total = a + b + c + d
        peak_total = a + c
        numerator_total = a + b
        distributions[table] = {
            x: Fraction(
                comb(peak_total, x)
                * comb(total - peak_total, numerator_total - x),
                comb(total, numerator_total),
            )
            for x in range(
                max(0, numerator_total - (total - peak_total)),
                min(numerator_total, peak_total) + 1,
            )
        }

    for left in strata:
        for right in strata:
            exact = sum(
                (left_probability * right_probability
                 for left_value, left_probability in distributions[left].items()
                 for right_value, right_probability in distributions[right].items()
                 if left_value + right_value >= left[0] + right[0]),
                Fraction(0),
            )
            computed, floored = _exact_conditional_tail((left, right))
            assert not floored
            assert Fraction.from_float(computed) >= exact, (left, right)
            assert computed == pytest.approx(float(exact), rel=1e-12), (left, right)


def test_swapping_direction_reciprocates_descriptive_common_odds_ratio():
    inputs = matched_inputs([
        ("treated", "1", [("site", 3, 1)]),
        ("control", "1", [("site", 1, 3)]),
        ("treated", "2", [("site", 2, 1)]),
        ("control", "2", [("site", 0, 3)]),
    ])
    forward = analyze_matched_condition(inputs, "TIS", "treated", "control").iloc[0]
    reverse = analyze_matched_condition(inputs, "TIS", "control", "treated").iloc[0]
    assert reverse.common_odds_ratio == pytest.approx(1 / forward.common_odds_ratio)
    assert reverse.concordant_pair_count == forward.discordant_pair_count
    assert reverse.discordant_pair_count == forward.concordant_pair_count


def test_exactly_neutral_common_effect_cannot_round_into_support():
    # The two numerator-favoring contributions sum exactly to the one
    # denominator-favoring contribution: 1802/87 + 1736/87 = 7442/183 = 122/3.
    # Naive binary64 accumulation rounds the left side above the right side.
    inputs = matched_inputs([
        ("treated", "1", [("site", 53, 0)]),
        ("control", "1", [("site", 0, 34)]),
        ("treated", "2", [("site", 56, 0)]),
        ("control", "2", [("site", 0, 31)]),
        ("treated", "3", [("site", 0, 61)]),
        ("control", "3", [("site", 122, 0)]),
    ])

    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh", fdr=0.6,
    ).iloc[0]

    assert row.q_value < 0.6
    assert row.concordant_pair_count == 2
    assert row.discordant_pair_count == 1
    assert row.common_odds_ratio == 1
    assert row.log2_common_odds_ratio == 0
    assert not row.significant


def test_no_information_stays_missing_and_is_excluded_from_correction_family():
    strong = [("strong", 10, 0), ("empty", 0, 0)]
    control = [("strong", 0, 10), ("empty", 0, 0)]
    inputs = matched_inputs([
        ("treated", "1", strong), ("control", "1", control),
        ("treated", "2", strong), ("control", "2", control),
    ])
    result = analyze_matched_condition(inputs, "TIS", "treated", "control", correction="by")
    strong_row = result.loc[result.candidate_id == "strong"].iloc[0]
    empty_row = result.loc[result.candidate_id == "empty"].iloc[0]
    assert strong_row.common_odds_ratio == np.inf
    assert strong_row.log2_common_odds_ratio == np.inf
    assert strong_row.p_value < 1e-10
    assert bool(strong_row.significant)
    assert strong_row.tests_in_family == 1
    assert empty_row.fit_status == "no_informative_pairs"
    assert empty_row.informative_pair_count == 0
    assert pd.isna(empty_row.p_value)
    assert pd.isna(empty_row.q_value)
    assert pd.isna(empty_row.significant)
    assert str(result.q_value.dtype) == "Float64"
    assert str(result.significant.dtype) == "boolean"


def test_one_informative_pair_is_descriptive_only():
    inputs = matched_inputs([
        ("treated", "1", [("site", 100, 0)]),
        ("control", "1", [("site", 0, 100)]),
        ("treated", "2", [("site", 0, 0)]),
        ("control", "2", [("site", 0, 0)]),
    ])
    row = analyze_matched_condition(inputs, "TIS", "treated", "control").iloc[0]
    assert row.informative_pair_count == 1
    assert row.fit_status == "insufficient_informative_pairs"
    assert row.common_odds_ratio == np.inf
    assert pd.isna(row.p_value)
    assert pd.isna(row.q_value)
    assert pd.isna(row.significant)
    assert row.tests_in_family == 0


def test_support_requires_strict_majority_directional_agreement():
    inputs = matched_inputs([
        ("treated", "1", [("site", 100, 0)]),
        ("control", "1", [("site", 0, 100)]),
        ("treated", "2", [("site", 0, 2)]),
        ("control", "2", [("site", 1, 1)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
    ).iloc[0]
    assert row.p_value < 0.05
    assert row.common_odds_ratio > 1
    assert row.concordant_pair_count == 1
    assert row.discordant_pair_count == 1
    assert row.tied_pair_count == 0
    assert not row.significant


def test_tied_informative_pairs_are_reported_and_fail_majority_gate():
    inputs = matched_inputs([
        ("treated", "1", [("site", 20, 0)]),
        ("control", "1", [("site", 0, 20)]),
        ("treated", "2", [("site", 1, 1)]),
        ("control", "2", [("site", 1, 1)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
    ).iloc[0]
    assert row.concordant_pair_count == 1
    assert row.discordant_pair_count == 0
    assert row.tied_pair_count == 1
    assert not row.significant


def test_extreme_nonzero_tail_is_floored_instead_of_reported_as_zero():
    inputs = matched_inputs([
        ("treated", "1", [("site", 1000, 0)]),
        ("control", "1", [("site", 0, 1000)]),
        ("treated", "2", [("site", 1000, 0)]),
        ("control", "2", [("site", 0, 1000)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
    ).iloc[0]
    assert row.p_value == MIN_POSITIVE_P_VALUE
    assert row.q_value > 0
    assert row.fit_status == "p_value_floored"


def test_subnormal_rounding_boundary_is_marked_as_floored():
    inputs = matched_inputs([
        ("treated", "1", [("site", 271, 0)]),
        ("control", "1", [("site", 0, 271)]),
        ("treated", "2", [("site", 271, 0)]),
        ("control", "2", [("site", 0, 271)]),
    ])
    row = analyze_matched_condition(inputs, "TIS", "treated", "control").iloc[0]
    assert row.p_value == MIN_POSITIVE_P_VALUE
    assert row.fit_status == "p_value_floored"


def test_narrow_longdouble_platform_refuses_underflow_prone_convolution(monkeypatch):
    monkeypatch.setattr(matched_statistics, "_LONGDOUBLE_HAS_WIDER_EXPONENT", False)
    inputs = matched_inputs([
        ("treated", "1", [("site", 271, 0)]),
        ("control", "1", [("site", 0, 271)]),
        ("treated", "2", [("site", 271, 0)]),
        ("control", "2", [("site", 0, 271)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]
    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


def test_work_limit_is_independent_of_replicate_labels():
    inputs = matched_inputs([
        ("treated", "1", [("site", 998, 0)]),
        ("control", "1", [("site", 0, 998)]),
        ("treated", "2", [("site", 2000, 0)]),
        ("control", "2", [("site", 0, 2000)]),
    ])
    relabeled = {
        (assay, condition, "2" if replicate == "1" else "1"): frame
        for (assay, condition, replicate), frame in inputs.items()
    }
    first = analyze_matched_condition(inputs, "TIS", "treated", "control").iloc[0]
    second = analyze_matched_condition(relabeled, "TIS", "treated", "control").iloc[0]
    assert first.fit_status == second.fit_status == "p_value_floored"
    assert first.p_value == second.p_value == MIN_POSITIVE_P_VALUE


def test_unbounded_exact_work_is_refused_with_explicit_status():
    inputs = matched_inputs([
        ("treated", "1", [("site", 4000, 0)]),
        ("control", "1", [("site", 0, 4000)]),
        ("treated", "2", [("site", 4000, 0)]),
        ("control", "2", [("site", 0, 4000)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]
    assert row.fit_status == "computational_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


def test_hypergeometric_total_bound_is_inclusive():
    maximum = MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL
    inputs = matched_inputs([
        ("treated", "1", [("site", 1, 0)]),
        ("control", "1", [("site", 0, maximum - 1)]),
        ("treated", "2", [("site", 1, 0)]),
        ("control", "2", [("site", 0, maximum - 1)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
    ).iloc[0]

    assert row.fit_status == "ok"
    assert row.p_value == pytest.approx(1 / maximum**2)
    assert row.tests_in_family == 1


def test_large_sparse_margins_match_exact_fixed_margin_fraction():
    maximum = MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL
    table = (2, 1, 0, maximum - 3)
    exact_one = Fraction(6, maximum * (maximum - 1))

    computed, floored = _exact_conditional_tail((table, table))

    assert not floored
    assert Fraction.from_float(computed) >= exact_one**2
    assert computed == pytest.approx(float(exact_one**2), rel=1e-12)


def test_over_bound_hypergeometric_stratum_stays_in_family_independent_of_order():
    over_bound = MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL + 1
    inputs = matched_inputs([
        ("treated", "1", [("ordinary", 4, 0), ("limit", 1, 0)]),
        ("control", "1", [("ordinary", 0, 4), ("limit", 0, over_bound - 1)]),
        ("treated", "2", [("ordinary", 4, 0), ("limit", 2, 1)]),
        ("control", "2", [("ordinary", 0, 4), ("limit", 1, 2)]),
    ])
    relabeled = {
        (assay, condition, "2" if replicate == "1" else "1"): frame
        for (assay, condition, replicate), frame in reversed(list(inputs.items()))
    }

    first = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh", fdr=1,
    ).set_index("candidate_id")
    second = analyze_matched_condition(
        relabeled, "TIS", "treated", "control", correction="bh", fdr=1,
    ).set_index("candidate_id")

    assert first.at["limit", "fit_status"] == "numeric_range_limit"
    assert first.at["limit", "p_value"] == 1
    assert first.at["limit", "q_value"] == 1
    assert not first.at["limit", "significant"]
    assert set(first.tests_in_family) == {2}
    pd.testing.assert_frame_equal(first, second)


def test_huge_python_counts_fail_safe_and_remain_exact_in_outputs(tmp_path):
    maximum_int64 = (1 << 63) - 1
    inputs = matched_inputs([
        ("treated", "1", [("site", maximum_int64, 1)]),
        ("control", "1", [("site", 1, maximum_int64)]),
        ("treated", "2", [("site", maximum_int64, 1)]),
        ("control", "2", [("site", 1, maximum_int64)]),
    ])

    result, strata = analyze_matched_condition(
        inputs, "TIS", "treated", "control", correction="bh",
        return_strata=True,
    )
    row = result.iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert not row.significant
    assert pd.isna(row.common_odds_ratio)
    assert row.numerator_peak_count == 2 * maximum_int64
    assert row.denominator_background_count == 2 * maximum_int64
    assert result.numerator_peak_count.dtype == object
    assert result.denominator_background_count.dtype == object
    assert strata.numerator_peak_count.tolist() == [maximum_int64, maximum_int64]
    audit_path = tmp_path / "matched.tsv"
    result.to_csv(audit_path, sep="\t", index=False)
    assert str(2 * maximum_int64) in audit_path.read_text()

    attached = attach_matched_to_orfs(
        pd.DataFrame([{"Genome": "chr", "Start": 10, "Stop": 30, "Strand": "+"}]),
        result,
        "TIS",
        "huge",
    )
    assert attached.loc[0, "huge_numerator_peak_count"] == 2 * maximum_int64
    assert attached.huge_numerator_peak_count.dtype == object


def test_optional_long_table_preserves_each_pair_for_reconstruction():
    inputs = matched_inputs([
        ("treated", "2", [("site", 2, 1)]),
        ("control", "2", [("site", 0, 3)]),
        ("treated", "1", [("site", 3, 1)]),
        ("control", "1", [("site", 1, 3)]),
    ])
    result, strata = analyze_matched_condition(
        inputs, "TIS", "treated", "control", return_strata=True,
    )
    assert len(result) == 1
    assert set(strata.assay) == {"TIS"}
    assert set(strata.numerator_condition) == {"treated"}
    assert set(strata.denominator_condition) == {"control"}
    assert strata.replicate_id.tolist() == ["1", "2"]
    assert strata.numerator_peak_count.tolist() == [3, 2]
    assert strata.denominator_peak_count.tolist() == [1, 0]
    assert strata.informative.tolist() == [True, True]
    assert strata.pair_direction.tolist() == ["numerator", "numerator"]


def test_empty_but_schema_complete_candidate_universe_is_stable():
    empty = candidate_frame([("unused", 0, 0)]).iloc[:0]
    inputs = {
        ("TIS", condition, replicate): empty.copy()
        for condition in ("treated", "control")
        for replicate in ("1", "2")
    }
    result = analyze_matched_condition(inputs, "TIS", "treated", "control")
    assert result.empty
    assert {"candidate_id", "p_value", "q_value", "significant", "fit_status"} <= set(result.columns)
    assert str(result.tests_in_family.dtype) == "Int64"
    assert str(result.q_value.dtype) == "Float64"
    assert str(result.significant.dtype) == "boolean"


def test_candidate_and_sample_order_do_not_change_results():
    second = {"Window_start": 20, "Window_stop": 24, "Codon_start": 22}
    inputs = matched_inputs([
        ("treated", "2", [("second", 1, 4, second), ("first", 4, 1)]),
        ("control", "2", [("second", 2, 3, second), ("first", 1, 4)]),
        ("treated", "1", [("first", 3, 2), ("second", 2, 3, second)]),
        ("control", "1", [("first", 1, 4), ("second", 2, 3, second)]),
    ])
    forward = analyze_matched_condition(inputs, "TIS", "treated", "control")
    reverse = analyze_matched_condition(dict(reversed(list(inputs.items()))), "TIS", "treated", "control")
    pd.testing.assert_frame_equal(forward, reverse)
    assert forward.candidate_id.tolist() == ["first", "second"]


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda data: data.pop(("TIS", "control", "2")), "identical replicate IDs"),
        (lambda data: [data.pop(key) for key in list(data) if key[2] == "2"], "at least two"),
        (
            lambda data: replace_second_control_table(data, lambda frame: frame.drop(index=0)),
            "candidate universes",
        ),
        (
            lambda data: replace_second_control_table(
                data, lambda frame: frame.assign(Codon=["GTG"]),
            ),
            "metadata",
        ),
        (
            lambda data: replace_second_control_table(
                data, lambda frame: frame.assign(peak_count=[1.5]),
            ),
            "raw integer",
        ),
        (
            lambda data: replace_second_control_table(
                data, lambda frame: frame.assign(background_count=[-1]),
            ),
            "raw integer",
        ),
    ],
)
def test_invalid_matched_inputs_are_rejected(mutate, message):
    inputs = matched_inputs([
        ("treated", "1", [("site", 3, 1)]), ("control", "1", [("site", 1, 3)]),
        ("treated", "2", [("site", 2, 1)]), ("control", "2", [("site", 1, 2)]),
    ])
    mutate(inputs)
    with pytest.raises(ValueError, match=message):
        analyze_matched_condition(inputs, "TIS", "treated", "control")


def test_invalid_assay_and_duplicate_input_columns_are_rejected():
    inputs = matched_inputs([
        ("treated", "1", [("site", 3, 1)]), ("control", "1", [("site", 1, 3)]),
        ("treated", "2", [("site", 2, 1)]), ("control", "2", [("site", 1, 2)]),
    ])
    with pytest.raises(ValueError, match="assay must be TIS or TTS"):
        analyze_matched_condition(inputs, "RIBO", "treated", "control")

    duplicate = inputs[("TIS", "treated", "1")]
    duplicate.columns = [*duplicate.columns[:-1], duplicate.columns[-2]]
    with pytest.raises(ValueError, match="duplicate column names"):
        analyze_matched_condition(inputs, "TIS", "treated", "control")


def test_bh_and_by_adjust_only_informative_candidates():
    rows_treated = [("one", 4, 0), ("two", 3, 1), ("empty", 0, 0)]
    rows_control = [("one", 0, 4), ("two", 1, 3), ("empty", 0, 0)]
    inputs = matched_inputs([
        ("treated", "1", rows_treated), ("control", "1", rows_control),
        ("treated", "2", rows_treated), ("control", "2", rows_control),
    ])
    bh = analyze_matched_condition(inputs, "TIS", "treated", "control", correction="bh")
    by = analyze_matched_condition(inputs, "TIS", "treated", "control", correction="by")
    assert set(bh.tests_in_family) == {2}
    assert set(by.tests_in_family) == {2}
    finite = bh.q_value.notna()
    assert (by.loc[finite, "q_value"] >= bh.loc[finite, "q_value"]).all()
    assert pd.isna(by.loc[by.candidate_id == "empty", "q_value"]).all()


def test_limit_case_remains_in_mixed_family_but_ineligible_candidate_does_not(monkeypatch):
    original_tail = matched_statistics._exact_conditional_tail

    def controlled_tail(strata):
        if strata[0][0] == 99:
            raise matched_statistics._ExactComputationLimit
        return original_tail(strata)

    monkeypatch.setattr(matched_statistics, "_exact_conditional_tail", controlled_tail)
    first_treated = [("ordinary", 4, 0), ("limit", 99, 0), ("ineligible", 4, 0)]
    first_control = [("ordinary", 0, 4), ("limit", 0, 99), ("ineligible", 0, 4)]
    second_treated = [("ordinary", 4, 0), ("limit", 99, 0), ("ineligible", 0, 0)]
    second_control = [("ordinary", 0, 4), ("limit", 0, 99), ("ineligible", 0, 0)]
    inputs = matched_inputs([
        ("treated", "1", first_treated), ("control", "1", first_control),
        ("treated", "2", second_treated), ("control", "2", second_control),
    ])
    ordinary_inputs = matched_inputs([
        ("treated", "1", [("ordinary", 4, 0)]),
        ("control", "1", [("ordinary", 0, 4)]),
        ("treated", "2", [("ordinary", 4, 0)]),
        ("control", "2", [("ordinary", 0, 4)]),
    ])

    for correction in ("bh", "by"):
        result = analyze_matched_condition(
            inputs, "TIS", "treated", "control", correction=correction, fdr=1,
        ).set_index("candidate_id")
        alone = analyze_matched_condition(
            ordinary_inputs, "TIS", "treated", "control", correction=correction,
        ).iloc[0]

        assert set(result.tests_in_family) == {2}
        assert result.at["limit", "fit_status"] == "computational_limit"
        assert result.at["limit", "p_value"] == 1
        assert result.at["limit", "q_value"] == 1
        assert not result.at["limit", "significant"]
        assert pd.isna(result.at["ineligible", "p_value"])
        assert pd.isna(result.at["ineligible", "q_value"])
        assert pd.isna(result.at["ineligible", "significant"])
        assert result.at["ordinary", "q_value"] >= alone.q_value


@pytest.mark.parametrize("invalid_p", [0.0, np.nan, np.inf, -0.1, 1.1])
def test_invalid_exact_backend_probability_fails_safe_in_family(monkeypatch, invalid_p):
    monkeypatch.setattr(
        matched_statistics, "_exact_conditional_tail", lambda strata: (invalid_p, False),
    )
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


@pytest.mark.parametrize("failure", ["nonfinite", "exception"])
def test_hypergeometric_recurrence_failure_stays_in_family_conservatively(
    monkeypatch, failure):
    if failure == "nonfinite":
        def backend(ratios):
            return np.full(np.shape(ratios), np.nan)
    else:
        def backend(*args):
            raise FloatingPointError("simulated backend failure")

    monkeypatch.setattr(matched_statistics.np, "cumprod", backend)
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


@pytest.mark.parametrize(
    "invalid_support",
    [(0.5, 4), (0, 3), (4, 0), (np.nan, 4)],
)
def test_invalid_hypergeometric_support_stays_in_family_conservatively(
    monkeypatch, invalid_support,
):
    monkeypatch.setattr(
        matched_statistics.hypergeom, "support", lambda *args: invalid_support,
    )
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


@pytest.mark.parametrize("failure", ["partial_nonfinite", "wrong_shape"])
def test_invalid_hypergeometric_recurrence_stays_in_family(
    monkeypatch, failure,
):
    if failure == "partial_nonfinite":
        def backend(ratios):
            result = np.ones(np.shape(ratios))
            result[0] = np.nan
            return result
    else:
        def backend(*args):
            return np.array([0.0])

    monkeypatch.setattr(matched_statistics.np, "cumprod", backend)
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


def test_invalid_convolution_stays_in_family_conservatively(monkeypatch):
    def zero_convolution(left, right):
        return np.zeros(len(left) + len(right) - 1)

    monkeypatch.setattr(matched_statistics.np, "convolve", zero_convolution)
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "numeric_range_limit"
    assert row.p_value == 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


@pytest.mark.parametrize("failure", ["nonfinite", "exception"])
def test_matched_correction_backend_failure_preserves_raw_p_and_fails_closed(
    monkeypatch, failure,
):
    if failure == "nonfinite":
        def backend(p_values, method):
            return np.full(np.shape(p_values), np.nan)
    else:
        def backend(*args, **kwargs):
            raise FloatingPointError("simulated correction backend failure")

    monkeypatch.setattr(matched_statistics, "false_discovery_control", backend)
    inputs = matched_inputs([
        ("treated", "1", [("site", 4, 0)]),
        ("control", "1", [("site", 0, 4)]),
        ("treated", "2", [("site", 4, 0)]),
        ("control", "2", [("site", 0, 4)]),
    ])
    row = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    ).iloc[0]

    assert row.fit_status == "correction_limit"
    assert 0 < row.p_value < 1
    assert row.q_value == 1
    assert row.tests_in_family == 1
    assert not row.significant


def test_matched_materially_underadjusted_family_fails_closed(monkeypatch):
    monkeypatch.setattr(
        matched_statistics, "false_discovery_control",
        lambda p_values, method: np.zeros(np.shape(p_values)),
    )
    inputs = matched_inputs([
        ("treated", "1", [("site1", 4, 0), ("site2", 5, 0)]),
        ("control", "1", [("site1", 0, 4), ("site2", 0, 5)]),
        ("treated", "2", [("site1", 4, 0), ("site2", 5, 0)]),
        ("control", "2", [("site1", 0, 4), ("site2", 0, 5)]),
    ])
    rows = analyze_matched_condition(
        inputs, "TIS", "treated", "control", fdr=1,
    )

    assert len(rows) == 2
    assert rows.p_value.between(0, 1, inclusive="neither").all()
    assert rows.q_value.eq(1).all()
    assert rows.fit_status.eq("correction_limit").all()
    assert not rows.significant.any()


def test_by_adjustment_does_not_round_below_exact_input_float_result():
    p_values = np.array([0.001, 0.005])
    q_values = matched_statistics._adjust_p_values(p_values, "by")
    exact_larger = Fraction.from_float(float(p_values[1])) * 3 / 2

    assert Fraction.from_float(float(q_values[1])) >= exact_larger


@pytest.mark.parametrize(
    "method,strand,codon_start,start,stop",
    [
        ("TIS", "+", 10, 10, 30),
        ("TIS", "-", 28, 10, 30),
        ("TTS", "+", 28, 10, 30),
        ("TTS", "-", 10, 10, 30),
    ],
)
def test_attach_matches_strand_specific_tis_and_tts_boundaries(method, strand, codon_start, start, stop):
    matched = pd.DataFrame([{
        "Genome": "chr", "Codon_start": codon_start, "Strand": strand, "assay": method,
        "matched_pair_count": 2, "informative_pair_count": 2, "concordant_pair_count": 1,
        "discordant_pair_count": 1, "tied_pair_count": 0,
        "numerator_peak_count": 9, "numerator_background_count": 4,
        "denominator_peak_count": 2, "denominator_background_count": 8,
        "common_odds_ratio": 8.0, "log2_common_odds_ratio": 3.0,
        "p_value": 0.01, "q_value": 0.02, "significant": True,
        "correction": "by", "tests_in_family": 100, "fit_status": "ok",
    }])
    orfs = pd.DataFrame([
        {"Genome": "chr", "Start": start, "Stop": stop, "Strand": strand},
        {"Genome": "chr", "Start": 100, "Stop": 120, "Strand": strand},
    ])
    result = attach_matched_to_orfs(orfs, matched, method, "treated_vs_control")
    assert result.loc[0, "treated_vs_control_common_odds_ratio"] == 8
    assert bool(result.loc[0, "treated_vs_control_significant"])
    assert pd.isna(result.loc[1, "treated_vs_control_q_value"])
    assert pd.isna(result.loc[1, "treated_vs_control_significant"])
    assert str(result["treated_vs_control_matched_pair_count"].dtype) == "Int64"
    assert str(result["treated_vs_control_significant"].dtype) == "boolean"


def test_attach_empty_orf_table_has_stable_schema():
    matched_columns = ["Genome", "Codon_start", "Strand", "assay", *ATTACH_FIELDS]
    result = attach_matched_to_orfs(
        pd.DataFrame(columns=["Genome", "Start", "Stop", "Strand"]),
        pd.DataFrame(columns=matched_columns),
        "TIS",
        "comparison",
    )
    assert result.empty
    assert all(f"comparison_{field}" in result for field in ATTACH_FIELDS)
    assert str(result["comparison_q_value"].dtype) == "Float64"
    assert str(result["comparison_significant"].dtype) == "boolean"


@pytest.mark.parametrize("prefix", ["", "bad prefix", "../escape", "x/y"])
def test_attach_rejects_unsafe_prefix(prefix):
    with pytest.raises(ValueError, match="prefix"):
        attach_matched_to_orfs(
            pd.DataFrame(columns=["Genome", "Start", "Stop", "Strand"]),
            pd.DataFrame(),
            "TIS",
            prefix,
        )


def test_attach_rejects_assay_mismatch_and_column_overwrite():
    matched = pd.DataFrame(columns=["Genome", "Codon_start", "Strand", "assay", *ATTACH_FIELDS])
    matched.loc[0] = ["chr", 10, "+", "TTS", *([pd.NA] * len(ATTACH_FIELDS))]
    orfs = pd.DataFrame([{"Genome": "chr", "Start": 10, "Stop": 30, "Strand": "+"}])
    with pytest.raises(ValueError, match="does not match"):
        attach_matched_to_orfs(orfs, matched, "TIS", "comparison")

    matched.at[0, "assay"] = "TIS"
    orfs["comparison_q_value"] = 0.5
    with pytest.raises(ValueError, match="already exist"):
        attach_matched_to_orfs(orfs, matched, "TIS", "comparison")


def test_attach_rejects_missing_assay_and_noninteger_orf_coordinates():
    matched = pd.DataFrame(columns=["Genome", "Codon_start", "Strand", "assay", *ATTACH_FIELDS])
    matched.loc[0] = ["chr", 10, "+", pd.NA, *([pd.NA] * len(ATTACH_FIELDS))]
    orfs = pd.DataFrame([{"Genome": "chr", "Start": 10, "Stop": 30, "Strand": "+"}])
    with pytest.raises(ValueError, match="does not match"):
        attach_matched_to_orfs(orfs, matched, "TIS", "comparison")

    matched.at[0, "assay"] = "TIS"
    orfs = pd.DataFrame([{"Genome": "chr", "Start": 10.5, "Stop": 30, "Strand": "+"}])
    with pytest.raises(ValueError, match="must be integers"):
        attach_matched_to_orfs(orfs, matched, "TIS", "comparison")
