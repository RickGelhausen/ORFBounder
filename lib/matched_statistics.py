"""Exploratory fixed-margin matched inference for local boundary counts.

This module compares two conditions within one assay while retaining biological
replicates as matched strata.  Inputs are the candidate tables produced by
``lib.statistics.analyze_peak_enrichment`` (or equivalent data frames) in a
mapping keyed by ``(assay, condition, replicate)``.  Equivalent frames must
contain the columns in :data:`CANDIDATE_METADATA_FIELDS` plus nonnegative raw
integer ``peak_count`` and ``background_count`` columns.

For every candidate and replicate, the two conditions form a 2x2 table::

                         peak       background
        numerator          a              b
        denominator        c              d

Under a common odds ratio of one, and conditional on both row and column
margins, ``a`` is hypergeometric.  The null distributions from informative
replicate strata are convolved, giving an exact one-sided tail probability for
the sum of numerator peak counts.  This is not a pooled-replicate test: each
replicate keeps its own margins.  The reported Mantel-Haenszel odds ratio is a
descriptive common-effect estimate; it is not used to calculate the exact
p-value.  This endpoint-count model does not estimate biological dispersion or
a population condition effect.  A support flag additionally requires at least
two informative pairs and strict majority agreement in the predeclared
direction.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, hypergeom

from lib.multiple_testing import conservative_adjustment


CANDIDATE_METADATA_FIELDS = (
    "candidate_id",
    "Genome",
    "Window_start",
    "Window_stop",
    "Strand",
    "Codon",
    "Codon_start",
    "peak_width",
    "background_width",
)
TOPOLOGY_FIELDS = ("Reference_length", "Is_circular")

COUNT_FIELDS = ("peak_count", "background_count")

MATCHED_RESULT_FIELDS = (
    "assay",
    "numerator_condition",
    "denominator_condition",
    "matched_pair_count",
    "informative_pair_count",
    "concordant_pair_count",
    "discordant_pair_count",
    "tied_pair_count",
    "numerator_peak_count",
    "numerator_background_count",
    "denominator_peak_count",
    "denominator_background_count",
    "common_odds_ratio",
    "log2_common_odds_ratio",
    "p_value",
    "q_value",
    "significant",
    "correction",
    "tests_in_family",
    "fit_status",
)

MATCHED_STRATUM_FIELDS = (
    "assay",
    "numerator_condition",
    "denominator_condition",
    "replicate_id",
    "numerator_peak_count",
    "numerator_background_count",
    "denominator_peak_count",
    "denominator_background_count",
    "informative",
    "pair_direction",
)

MATCHED_MODEL_NAME = "exact_stratified_fixed_margin_endpoint_count"
MATCHED_ALTERNATIVE = "greater"
MIN_INFORMATIVE_PAIRS = 2
MATCHED_SUPPORT_RULE = (
    "q_value<=fdr and common_odds_ratio>1 and "
    "2*concordant_pair_count>informative_pair_count"
)
# Bound direct polynomial convolution rather than allowing a single very deep
# candidate to consume unbounded CPU and memory.
MAX_EXACT_CONVOLUTION_WORK = 2_000_000
# A deliberately conservative, platform-independent ceiling for every
# hypergeometric population (the grand total of one matched 2x2 stratum).
# scipy's discrete-distribution backend is integer/float backed and can fail or
# lose accuracy for otherwise legal, arbitrarily large Python integers.
MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL = (1 << 31) - 1
MIN_POSITIVE_P_VALUE = float(np.nextafter(0.0, 1.0))
_LONGDOUBLE_HAS_WIDER_EXPONENT = (
    np.finfo(np.longdouble).minexp < np.finfo(np.float64).minexp
)
MATCHED_NUMERIC_BACKEND = (
    "extended_longdouble"
    if _LONGDOUBLE_HAS_WIDER_EXPONENT
    else "float64_with_underflow_guard"
)
_MIN_POSITIVE_FLOAT_LOG = math.log(MIN_POSITIVE_P_VALUE)

# Fields copied to ORF tables. Conditions and assay remain in the complete
# candidate table and are represented by the caller-supplied column prefix.
ATTACH_FIELDS = MATCHED_RESULT_FIELDS[3:]

_INTEGER_RESULT_FIELDS = {
    "matched_pair_count",
    "informative_pair_count",
    "concordant_pair_count",
    "discordant_pair_count",
    "tied_pair_count",
    "numerator_peak_count",
    "numerator_background_count",
    "denominator_peak_count",
    "denominator_background_count",
    "tests_in_family",
}
_FLOAT_RESULT_FIELDS = {
    "common_odds_ratio",
    "log2_common_odds_ratio",
    "p_value",
    "q_value",
}
_STRING_RESULT_FIELDS = {
    "assay",
    "numerator_condition",
    "denominator_condition",
    "correction",
    "fit_status",
}
_PREFIX_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_CONDITION_PATTERN = re.compile(r"[A-Za-z0-9]+\Z")
_SIGNED_INT64_MIN = -(1 << 63)
_SIGNED_INT64_MAX = (1 << 63) - 1


@dataclass(frozen=True, slots=True)
class MatchedComparison:
    """One predeclared, directional, same-assay condition comparison."""

    name: str
    assay: str
    numerator_condition: str
    denominator_condition: str

    @property
    def prefix(self):
        return f"{self.name}_{self.assay}"


def read_matched_comparisons(path) -> list[MatchedComparison]:
    """Read a small TSV declaring matched comparisons for a batch analysis."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Matched-comparison file does not exist: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if frame.columns.has_duplicates:
        raise ValueError("Matched-comparison file has duplicate column names.")
    frame = frame.apply(lambda column: column.str.strip())
    required = {"comparison", "assay", "numerator_condition", "denominator_condition"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Matched-comparison file is missing columns: {', '.join(missing)}.")
    unknown = sorted(set(frame.columns) - required)
    if unknown:
        raise ValueError(f"Unknown matched-comparison columns: {', '.join(unknown)}.")
    if frame.empty:
        raise ValueError("Matched-comparison file contains no comparisons.")

    comparisons = []
    seen_prefixes = set()
    for index, row in frame.iterrows():
        name = row["comparison"]
        assay = row["assay"].upper()
        numerator = row["numerator_condition"]
        denominator = row["denominator_condition"]
        if _PREFIX_PATTERN.fullmatch(name) is None:
            raise ValueError(f"Matched-comparison row {index + 2}: invalid comparison name {name!r}.")
        if assay not in {"TIS", "TTS"}:
            raise ValueError(f"Matched-comparison row {index + 2}: assay must be TIS or TTS.")
        prefix = (name, assay)
        if prefix in seen_prefixes:
            raise ValueError(
                "Matched-comparison output prefix must be unique: "
                f"{name!r} with assay {assay!r}."
            )
        if (_CONDITION_PATTERN.fullmatch(numerator) is None
                or _CONDITION_PATTERN.fullmatch(denominator) is None):
            raise ValueError(
                f"Matched-comparison row {index + 2}: conditions must match sample filename labels "
                "and contain only letters and digits."
            )
        if numerator == denominator:
            raise ValueError(f"Matched-comparison row {index + 2}: numerator and denominator must differ.")
        seen_prefixes.add(prefix)
        comparisons.append(MatchedComparison(name, assay, numerator, denominator))
    return comparisons


def _require_label(value, name):
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty string.")


def _validate_raw_count(value, sample_key, candidate_id, field):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise ValueError(
            "Matched statistics require nonnegative raw integer counts; "
            f"got {value!r} for {sample_key!r}, candidate {candidate_id!r}, {field}."
        )


def _validate_candidate_frame(frame, sample_key):
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"Matched-statistics input for {sample_key!r} must be a pandas DataFrame.")
    if frame.columns.has_duplicates:
        raise ValueError(f"Matched-statistics input for {sample_key!r} has duplicate column names.")
    required = set(CANDIDATE_METADATA_FIELDS) | set(COUNT_FIELDS)
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            f"Matched-statistics input for {sample_key!r} is missing columns: {', '.join(missing)}."
        )
    if frame["candidate_id"].isna().any() or not frame["candidate_id"].map(lambda value: isinstance(value, str) and bool(value)).all():
        raise ValueError(f"candidate_id values for {sample_key!r} must be nonempty strings.")
    if frame["candidate_id"].duplicated().any():
        duplicate = frame.loc[frame["candidate_id"].duplicated(), "candidate_id"].iloc[0]
        raise ValueError(f"Duplicate candidate_id {duplicate!r} in matched-statistics input {sample_key!r}.")

    topology_columns = set(TOPOLOGY_FIELDS).intersection(frame.columns)
    if topology_columns and topology_columns != set(TOPOLOGY_FIELDS):
        raise ValueError(
            f"Matched-statistics input for {sample_key!r} must provide both "
            "Reference_length and Is_circular."
        )

    for row in frame.itertuples(index=False):
        values = row._asdict()
        candidate_id = values["candidate_id"]
        for field in CANDIDATE_METADATA_FIELDS[1:]:
            value = values[field]
            if pd.isna(value):
                raise ValueError(
                    f"Candidate metadata cannot be missing: {sample_key!r}, "
                    f"candidate {candidate_id!r}, {field}."
                )
        if values["Strand"] not in {"+", "-"}:
            raise ValueError(f"Invalid strand for candidate {candidate_id!r} in {sample_key!r}.")
        for field in ("Window_start", "Window_stop", "Codon_start", "peak_width", "background_width"):
            value = values[field]
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise ValueError(
                    f"Candidate metadata {field} must be an integer for {sample_key!r}, "
                    f"candidate {candidate_id!r}."
                )
        if values["Window_start"] < 1 or values["Window_stop"] < values["Window_start"]:
            raise ValueError(f"Invalid candidate window for {candidate_id!r} in {sample_key!r}.")
        if values["Codon_start"] < 1 or values["peak_width"] < 1 or values["background_width"] < 0:
            raise ValueError(f"Invalid candidate widths or codon coordinate for {candidate_id!r} in {sample_key!r}.")
        if topology_columns:
            reference_length = values["Reference_length"]
            is_circular = values["Is_circular"]
            if (isinstance(reference_length, bool) or not isinstance(reference_length, Integral)
                    or reference_length < 1):
                raise ValueError(f"Invalid Reference_length for candidate {candidate_id!r} in {sample_key!r}.")
            if not isinstance(is_circular, (bool, np.bool_)):
                raise ValueError(f"Is_circular must be boolean for candidate {candidate_id!r} in {sample_key!r}.")
            if values["Window_start"] > reference_length or values["Codon_start"] > reference_length:
                raise ValueError(f"Candidate coordinates exceed the reference for {candidate_id!r} in {sample_key!r}.")
            if is_circular:
                if values["Window_stop"] - values["Window_start"] + 1 > reference_length:
                    raise ValueError(f"Circular candidate spans more than one revolution: {candidate_id!r}.")
            elif values["Window_stop"] > reference_length:
                raise ValueError(f"Linear candidate exceeds the reference: {candidate_id!r}.")
        for field in COUNT_FIELDS:
            _validate_raw_count(values[field], sample_key, candidate_id, field)


def _metadata_equal(left, right):
    if pd.isna(left) and pd.isna(right):
        return True
    return left == right


def _validate_identical_universe(indexed_frames):
    reference_key, reference = next(iter(indexed_frames.items()))
    topology_present = set(TOPOLOGY_FIELDS).issubset(reference.columns)
    metadata_fields = (*CANDIDATE_METADATA_FIELDS, *TOPOLOGY_FIELDS) if topology_present else CANDIDATE_METADATA_FIELDS
    reference_ids = set(reference.index)
    for sample_key, frame in indexed_frames.items():
        if set(TOPOLOGY_FIELDS).issubset(frame.columns) != topology_present:
            raise ValueError(
                "Matched samples must either all provide reference topology or all omit it."
            )
        candidate_ids = set(frame.index)
        if candidate_ids != reference_ids:
            missing = sorted(reference_ids - candidate_ids)
            extra = sorted(candidate_ids - reference_ids)
            detail = []
            if missing:
                detail.append(f"missing {missing[:3]!r}")
            if extra:
                detail.append(f"extra {extra[:3]!r}")
            raise ValueError(
                "Matched samples must have identical candidate universes; "
                f"{sample_key!r} differs from {reference_key!r} ({'; '.join(detail)})."
            )
        for candidate_id in reference_ids:
            for field in metadata_fields[1:]:
                left = reference.at[candidate_id, field]
                right = frame.at[candidate_id, field]
                if not _metadata_equal(left, right):
                    raise ValueError(
                        "Matched samples must have identical candidate metadata; "
                        f"candidate {candidate_id!r}, field {field!r} differs between "
                        f"{reference_key!r} and {sample_key!r}."
                    )
    return reference, metadata_fields


class _ExactComputationLimit(RuntimeError):
    """Raised when exact convolution would exceed the declared work bound."""


class _ExactNumericRangeLimit(RuntimeError):
    """Raised rather than export a potentially anti-conservative underflowed tail."""


def _adjust_p_values(p_values, correction):
    """Return finite in-range adjusted values, or ``None`` on backend failure."""
    return conservative_adjustment(p_values, correction, false_discovery_control)


def _stratum_support(a, b, c, d):
    """Return hypergeometric parameters and integer support without allocating it."""
    numerator_total = a + b
    peak_total = a + c
    total = numerator_total + c + d
    if total > MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL:
        raise _ExactNumericRangeLimit(
            "Matched stratum total exceeds the supported hypergeometric bound "
            f"({MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL:,})."
        )
    expected_lower = max(0, numerator_total - (total - peak_total))
    expected_upper = min(numerator_total, peak_total)
    try:
        raw_lower, raw_upper = hypergeom.support(total, peak_total, numerator_total)
        if isinstance(raw_lower, (bool, np.bool_)) or isinstance(raw_upper, (bool, np.bool_)):
            raise TypeError("Boolean support endpoint")
        lower_float, upper_float = float(raw_lower), float(raw_upper)
        if not math.isfinite(lower_float) or not math.isfinite(upper_float):
            raise ValueError("Non-finite support endpoint")
        lower, upper = int(raw_lower), int(raw_upper)
        if lower != raw_lower or upper != raw_upper:
            raise ValueError("Non-integral support endpoint")
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise _ExactNumericRangeLimit(
            "Could not evaluate finite integer support for a matched stratum."
        ) from exc
    if (lower, upper) != (expected_lower, expected_upper):
        raise _ExactNumericRangeLimit(
            "The hypergeometric backend returned an invalid matched-stratum support."
        )
    return total, peak_total, numerator_total, lower, upper


def _stratum_distribution(a, b, c, d):
    """Return the support origin and PMF for one fixed-margin 2x2 table."""
    total, peak_total, numerator_total, lower, upper = _stratum_support(a, b, c, d)
    # A mode-centered ratio recurrence avoids the cancellation in scipy's
    # logpmf for large margins. Its successive ratios are rational functions
    # of the exact integer margins, so only their evaluation and accumulation
    # incur floating-point error. Keep the array and temporary ratio chunks
    # bounded by the caller's convolution-work limit.
    mode = min(upper, max(lower, (numerator_total + 1) * (peak_total + 1) // (total + 2)))
    probabilities = np.empty(upper - lower + 1, dtype=np.longdouble)
    probabilities[mode - lower] = 1
    chunk_size = 16_384
    try:
        for first in range(mode, upper, chunk_size):
            last = min(upper, first + chunk_size)
            x = np.arange(first, last, dtype=np.longdouble)
            ratios = ((peak_total - x) * (numerator_total - x)) / (
                (x + 1) * (total - peak_total - numerator_total + x + 1)
            )
            if not np.all(np.isfinite(ratios)) or np.any(ratios <= 0):
                raise _ExactNumericRangeLimit("Invalid upper hypergeometric recurrence ratio.")
            segment = probabilities[first - lower] * np.cumprod(ratios)
            if segment.shape != ratios.shape:
                raise _ExactNumericRangeLimit("Invalid upper hypergeometric recurrence shape.")
            probabilities[first + 1 - lower:last + 1 - lower] = segment
        for first in range(mode, lower, -chunk_size):
            last = max(lower, first - chunk_size)
            x = np.arange(first, last, -1, dtype=np.longdouble)
            ratios = (x * (total - peak_total - numerator_total + x)) / (
                (peak_total - x + 1) * (numerator_total - x + 1)
            )
            if not np.all(np.isfinite(ratios)) or np.any(ratios <= 0):
                raise _ExactNumericRangeLimit("Invalid lower hypergeometric recurrence ratio.")
            stop_index = last - 1 - lower if last > lower else None
            segment = probabilities[first - lower] * np.cumprod(ratios)
            if segment.shape != ratios.shape:
                raise _ExactNumericRangeLimit("Invalid lower hypergeometric recurrence shape.")
            probabilities[first - 1 - lower:stop_index:-1] = segment
        probability_sum = probabilities.sum(dtype=np.longdouble)
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise _ExactNumericRangeLimit(
            "Could not safely construct a stratum's hypergeometric distribution."
        ) from exc
    if (probabilities.shape != (upper - lower + 1,)
            or not np.all(np.isfinite(probabilities))
            or np.any(probabilities <= 0)
            or not np.isfinite(probability_sum) or probability_sum <= 0):
        raise _ExactNumericRangeLimit(
            "Could not safely normalize a stratum's hypergeometric distribution."
        )
    try:
        normalized = np.asarray(
            probabilities / probability_sum, dtype=np.longdouble,
        )
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise _ExactNumericRangeLimit(
            "Could not safely normalize a stratum's hypergeometric distribution."
        ) from exc
    if not np.all(np.isfinite(normalized)) or np.any(normalized <= 0):
        raise _ExactNumericRangeLimit(
            "Could not safely normalize a stratum's hypergeometric distribution."
        )
    return lower, normalized


def _exact_conditional_tail(strata, *, max_work=MAX_EXACT_CONVOLUTION_WORK):
    """Return ``(P(sum(A) >= observed), floored)`` for fixed margins.

    Direct convolution is deterministic and avoids FFT tail artefacts. Work is
    explicitly bounded. A positive probability below float64 range is
    conservatively exported as the smallest positive float, with ``floored``
    set to true.
    """
    distribution = np.array([1.0], dtype=np.longdouble)
    support_start = 0
    observed = 0
    work = 0
    prepared = []
    for table in strata:
        a, b, c, d = table
        _, _, _, lower, upper = _stratum_support(*table)
        prepared.append((upper - lower + 1, table))
    # Smallest supports first minimizes sequential direct-convolution work and,
    # crucially, makes the work-limit decision independent of replicate labels.
    prepared.sort(key=lambda item: (item[0], item[1]))
    for support_length, (a, b, c, d) in prepared:
        work += len(distribution) * support_length
        if work > max_work:
            raise _ExactComputationLimit(
                f"Exact matched convolution exceeds the work limit ({max_work:,})."
            )
        lower, probabilities = _stratum_distribution(a, b, c, d)
        if not _LONGDOUBLE_HAS_WIDER_EXPONENT:
            smallest_left = float(np.min(distribution))
            smallest_right = float(np.min(probabilities))
            if (smallest_left <= 0 or smallest_right <= 0
                    or math.log(smallest_left) + math.log(smallest_right)
                    <= _MIN_POSITIVE_FLOAT_LOG):
                raise _ExactNumericRangeLimit(
                    "This platform lacks the exponent range needed for exact convolution."
                )
        expected_length = len(distribution) + len(probabilities) - 1
        try:
            distribution = np.asarray(
                np.convolve(distribution, probabilities), dtype=np.longdouble,
            )
        except (ArithmeticError, TypeError, ValueError) as exc:
            raise _ExactNumericRangeLimit(
                "Could not safely convolve the matched conditional distribution."
            ) from exc
        if (distribution.shape != (expected_length,)
                or not np.all(np.isfinite(distribution))
                or np.any(distribution <= 0)):
            raise _ExactNumericRangeLimit(
                "The exact matched distribution underflowed or became invalid."
            )
        probability_sum = distribution.sum(dtype=np.longdouble)
        if not np.isfinite(probability_sum) or probability_sum <= 0:
            raise _ExactNumericRangeLimit(
                "Could not safely normalize the matched conditional distribution."
            )
        try:
            distribution /= probability_sum
        except (ArithmeticError, TypeError, ValueError) as exc:
            raise _ExactNumericRangeLimit(
                "Could not safely normalize the matched conditional distribution."
            ) from exc
        support_start += lower
        observed += a
    tail_index = observed - support_start
    if tail_index <= 0:
        return 1.0, False
    if tail_index >= len(distribution):
        raise _ExactNumericRangeLimit(
            "Observed matched statistic lies outside its conditional support."
        )
    tail = distribution[tail_index:].sum(dtype=np.longdouble)
    if not np.isfinite(tail) or tail <= 0 or tail > 1:
        raise _ExactNumericRangeLimit(
            "Could not evaluate a valid matched conditional tail."
        )
    if tail <= np.longdouble(MIN_POSITIVE_P_VALUE):
        return MIN_POSITIVE_P_VALUE, True
    # The fixed-margin law is exact, but its recurrence, normalization and
    # direct convolution are floating-point calculations. Allow a deliberately
    # generous operation-scaled upward tolerance before exporting the tail;
    # this is an engineering rounding guard, not a claim of exact arithmetic.
    relative_guard = np.longdouble(64) * np.finfo(np.longdouble).eps * (
        work + len(distribution) + len(strata)
    )
    tail = min(np.longdouble(1), tail * (1 + relative_guard))
    value = float(tail)
    if value == 0.0:
        return MIN_POSITIVE_P_VALUE, True
    if np.longdouble(value) < tail:
        value = float(np.nextafter(value, math.inf))
    return float(np.clip(value, MIN_POSITIVE_P_VALUE, 1.0)), False


def _common_odds_ratio(strata):
    """Return the common odds ratio and its exact direction around one.

    The ratio is descriptive floating-point output.  Its support-rule direction
    is evaluated as the exact sign of ``sum((a*d - b*c) / total)`` so rounding
    at the neutral boundary cannot turn a null or denominator-favoring common
    effect into numerator support.
    """
    numerator_terms = []
    denominator_terms = []
    difference_numerator = 0
    difference_denominator = 1
    for a, b, c, d in strata:
        total = a + b + c + d
        numerator_terms.append((a * d) / total)
        denominator_terms.append((b * c) / total)
        difference_numerator = (
            difference_numerator * total
            + (a * d - b * c) * difference_denominator
        )
        difference_denominator *= total
        divisor = math.gcd(abs(difference_numerator), difference_denominator)
        difference_numerator //= divisor
        difference_denominator //= divisor

    numerator = math.fsum(numerator_terms)
    denominator = math.fsum(denominator_terms)
    if numerator == 0 and denominator == 0:
        return math.nan, 0
    if denominator == 0:
        return math.inf, 1
    ratio = numerator / denominator
    if difference_numerator == 0:
        return 1.0, 0
    if difference_numerator > 0:
        # Keep the exported value consistent with the exact support-rule test
        # even if a very small positive effect rounds to one in binary64.
        return max(ratio, math.nextafter(1.0, math.inf)), 1
    return min(ratio, math.nextafter(1.0, 0.0)), -1


def _nullable_integer_array(values):
    """Use nullable Int64 when exact; retain larger Python integers as objects."""
    values = list(values)
    if any(
        value is not pd.NA and pd.notna(value)
        and not _SIGNED_INT64_MIN <= int(value) <= _SIGNED_INT64_MAX
        for value in values
    ):
        # Aggregate raw counts can exceed signed int64 even when every input
        # cell is a valid Python integer. Object dtype preserves those counts
        # exactly in the audit TSV instead of wrapping or converting to float.
        return pd.array(
            [pd.NA if value is pd.NA or pd.isna(value) else int(value) for value in values],
            dtype=object,
        )
    return pd.array(values, dtype="Int64")


def _cast_result_columns(result):
    for field in _INTEGER_RESULT_FIELDS:
        result[field] = _nullable_integer_array(result[field])
    for field in _FLOAT_RESULT_FIELDS:
        result[field] = pd.array(result[field], dtype="Float64")
    for field in _STRING_RESULT_FIELDS:
        result[field] = pd.array(result[field], dtype="string")
    result["significant"] = pd.array(result["significant"], dtype="boolean")
    return result


def _cast_stratum_columns(result):
    for field in (
        "numerator_peak_count", "numerator_background_count",
        "denominator_peak_count", "denominator_background_count",
    ):
        result[field] = _nullable_integer_array(result[field])
    for field in ("assay", "numerator_condition", "denominator_condition", "replicate_id"):
        result[field] = pd.array(result[field], dtype="string")
    result["informative"] = pd.array(result["informative"], dtype="boolean")
    result["pair_direction"] = pd.array(result["pair_direction"], dtype="string")
    return result


def analyze_matched_condition(
    sample_statistics: Mapping[tuple[str, str, str], pd.DataFrame],
    assay: str,
    numerator_condition: str,
    denominator_condition: str,
    *,
    fdr: float = 0.05,
    correction: str = "by",
    return_strata: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame]:
    """Compare two conditions within an assay using matched replicate strata.

    Replicate sets must match exactly and contain at least two identifiers.  A
    stratum is informative when both conditions have at least one local read and
    the combined table has at least one peak and one background read.  Empty or
    otherwise fixed strata are retained in aggregate counts but do not change
    the conditional distribution. Candidates need at least two informative
    strata for p/q/support values. With ``return_strata=True``, return the
    complete candidate result and a long, per-pair count table that can
    reconstruct the conditional test.

    The alternative is predeclared and one-sided: the numerator condition has
    greater peak-versus-background odds than the denominator.  ``concordant``
    counts informative pairs whose observed cross product satisfies ``a*d >
    b*c``. Every candidate with at least two informative pairs remains in the
    correction family. A numerical or work-limit failure contributes the
    conservative value p=1 rather than shrinking that family.
    """
    if not isinstance(sample_statistics, Mapping):
        raise TypeError("sample_statistics must be a mapping keyed by (assay, condition, replicate).")
    _require_label(assay, "assay")
    if assay not in {"TIS", "TTS"}:
        raise ValueError("Matched-statistics assay must be TIS or TTS.")
    _require_label(numerator_condition, "numerator_condition")
    _require_label(denominator_condition, "denominator_condition")
    if numerator_condition == denominator_condition:
        raise ValueError("Numerator and denominator conditions must be different.")
    if isinstance(fdr, bool) or not isinstance(fdr, Real) or not math.isfinite(fdr) or not 0 < fdr <= 1:
        raise ValueError("Matched-statistics FDR must be greater than 0 and at most 1.")
    if correction not in {"bh", "by"}:
        raise ValueError("Matched-statistics correction must be 'bh' or 'by'.")
    if not isinstance(return_strata, bool):
        raise TypeError("return_strata must be boolean.")

    selected = {}
    for key, frame in sample_statistics.items():
        if (not isinstance(key, tuple) or len(key) != 3
                or any(not isinstance(part, str) or not part for part in key)):
            raise ValueError("Matched-statistics keys must be nonempty (assay, condition, replicate) string tuples.")
        key_assay, condition, replicate = key
        if key_assay == assay and condition in {numerator_condition, denominator_condition}:
            _validate_candidate_frame(frame, key)
            selected[(condition, replicate)] = frame

    numerator_replicates = {replicate for condition, replicate in selected if condition == numerator_condition}
    denominator_replicates = {replicate for condition, replicate in selected if condition == denominator_condition}
    if numerator_replicates != denominator_replicates:
        missing_numerator = sorted(denominator_replicates - numerator_replicates)
        missing_denominator = sorted(numerator_replicates - denominator_replicates)
        raise ValueError(
            "Matched conditions require identical replicate IDs; "
            f"missing numerator replicates={missing_numerator}, "
            f"missing denominator replicates={missing_denominator}."
        )
    if len(numerator_replicates) < 2:
        raise ValueError("Matched-condition statistics require at least two matched replicate IDs.")

    replicates = sorted(numerator_replicates)
    indexed_frames = {
        (condition, replicate): selected[(condition, replicate)].set_index("candidate_id", drop=False)
        for condition in (numerator_condition, denominator_condition)
        for replicate in replicates
    }
    reference, metadata_fields = _validate_identical_universe(indexed_frames)
    candidate_ids = sorted(
        reference.index,
        key=lambda candidate_id: (
            str(reference.at[candidate_id, "Genome"]),
            int(reference.at[candidate_id, "Window_start"]),
            int(reference.at[candidate_id, "Window_stop"]),
            str(reference.at[candidate_id, "Strand"]),
            candidate_id,
        ),
    )

    records = []
    common_effect_directions = []
    stratum_records = [] if return_strata else None
    for candidate_id in candidate_ids:
        all_strata = []
        informative_strata = []
        concordant = 0
        discordant = 0
        tied = 0
        for replicate in replicates:
            numerator = indexed_frames[(numerator_condition, replicate)].loc[candidate_id]
            denominator = indexed_frames[(denominator_condition, replicate)].loc[candidate_id]
            a = int(numerator["peak_count"])
            b = int(numerator["background_count"])
            c = int(denominator["peak_count"])
            d = int(denominator["background_count"])
            all_strata.append((a, b, c, d))
            informative = (
                a + b > 0 and c + d > 0
                and a + c > 0 and b + d > 0
            )
            cross_product_difference = a * d - b * c
            direction = "uninformative"
            if informative:
                informative_strata.append((a, b, c, d))
                if cross_product_difference > 0:
                    concordant += 1
                    direction = "numerator"
                elif cross_product_difference < 0:
                    discordant += 1
                    direction = "denominator"
                else:
                    tied += 1
                    direction = "tie"
            if stratum_records is not None:
                stratum_records.append({
                    **{
                        field: reference.at[candidate_id, field]
                        for field in metadata_fields
                    },
                    "assay": assay,
                    "numerator_condition": numerator_condition,
                    "denominator_condition": denominator_condition,
                    "replicate_id": replicate,
                    "numerator_peak_count": a,
                    "numerator_background_count": b,
                    "denominator_peak_count": c,
                    "denominator_background_count": d,
                    "informative": informative,
                    "pair_direction": direction,
                })

        metadata = {field: reference.at[candidate_id, field] for field in metadata_fields}
        record = {
            **metadata,
            "assay": assay,
            "numerator_condition": numerator_condition,
            "denominator_condition": denominator_condition,
            "matched_pair_count": len(replicates),
            "informative_pair_count": len(informative_strata),
            "concordant_pair_count": concordant,
            "discordant_pair_count": discordant,
            "tied_pair_count": tied,
            "numerator_peak_count": sum(table[0] for table in all_strata),
            "numerator_background_count": sum(table[1] for table in all_strata),
            "denominator_peak_count": sum(table[2] for table in all_strata),
            "denominator_background_count": sum(table[3] for table in all_strata),
            "common_odds_ratio": pd.NA,
            "log2_common_odds_ratio": pd.NA,
            "p_value": pd.NA,
            "q_value": pd.NA,
            "significant": pd.NA,
            "correction": correction,
            "tests_in_family": 0,
            "fit_status": "no_informative_pairs",
        }
        common_effect_direction = 0
        if informative_strata:
            exceeds_numeric_range = any(
                sum(table) > MAX_SUPPORTED_HYPERGEOMETRIC_TOTAL
                for table in informative_strata
            )
            if not exceeds_numeric_range:
                odds_ratio, common_effect_direction = _common_odds_ratio(informative_strata)
                record["common_odds_ratio"] = odds_ratio
                record["log2_common_odds_ratio"] = (
                    math.log2(odds_ratio) if odds_ratio > 0 and math.isfinite(odds_ratio)
                    else math.inf if odds_ratio == math.inf
                    else -math.inf if odds_ratio == 0
                    else pd.NA
                )
            if len(informative_strata) < MIN_INFORMATIVE_PAIRS:
                record["fit_status"] = "insufficient_informative_pairs"
            elif exceeds_numeric_range:
                # Check before both the float-backed descriptive effect and the
                # scipy hypergeometric backend. The candidate is still eligible
                # for multiple testing, so retain it conservatively with p=1.
                record["fit_status"] = "numeric_range_limit"
                record["p_value"] = 1.0
            else:
                try:
                    p_value, floored = _exact_conditional_tail(informative_strata)
                except _ExactComputationLimit:
                    record["fit_status"] = "computational_limit"
                    record["p_value"] = 1.0
                except _ExactNumericRangeLimit:
                    record["fit_status"] = "numeric_range_limit"
                    record["p_value"] = 1.0
                else:
                    if (not isinstance(p_value, Real) or not math.isfinite(float(p_value))
                            or not 0 < float(p_value) <= 1):
                        record["p_value"] = 1.0
                        record["fit_status"] = "numeric_range_limit"
                    else:
                        record["p_value"] = float(p_value)
                        record["fit_status"] = "p_value_floored" if floored else "ok"
        records.append(record)
        common_effect_directions.append(common_effect_direction)

    columns = [*metadata_fields, *MATCHED_RESULT_FIELDS]
    result = pd.DataFrame.from_records(records, columns=columns)
    finite_indices = [
        index for index, value in enumerate(result["p_value"])
        if value is not pd.NA and pd.notna(value) and math.isfinite(float(value))
    ]
    tests_in_family = len(finite_indices)
    if finite_indices:
        p_values = result.loc[finite_indices, "p_value"].to_numpy(dtype=float)
        q_values = _adjust_p_values(p_values, correction)
        if q_values is None:
            for index in finite_indices:
                result.at[index, "q_value"] = 1.0
                result.at[index, "significant"] = False
                if result.at[index, "fit_status"] not in {
                        "computational_limit", "numeric_range_limit"}:
                    result.at[index, "fit_status"] = "correction_limit"
        else:
            for index, q_value in zip(finite_indices, q_values):
                q_value = max(float(q_value), MIN_POSITIVE_P_VALUE)
                result.at[index, "q_value"] = q_value
                concordant = int(result.at[index, "concordant_pair_count"])
                informative = int(result.at[index, "informative_pair_count"])
                result.at[index, "significant"] = bool(
                    result.at[index, "fit_status"] in {"ok", "p_value_floored"}
                    and q_value <= fdr and common_effect_directions[index] > 0
                    and 2 * concordant > informative
                )
    result["tests_in_family"] = tests_in_family
    result = _cast_result_columns(result)
    if not return_strata:
        return result
    stratum_columns = [*metadata_fields, *MATCHED_STRATUM_FIELDS]
    strata = pd.DataFrame.from_records(stratum_records, columns=stratum_columns)
    return result, _cast_stratum_columns(strata)


def attach_matched_to_orfs(
    result_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    method: str,
    prefix: str,
) -> pd.DataFrame:
    """Attach matched-condition fields to the appropriate ORF boundary.

    ``method`` is ``TIS`` or ``TTS``. Candidate coordinates are interpreted in
    the same 1-based, strand-aware way as :func:`lib.statistics.attach_to_orfs`.
    Output columns are named ``<prefix>_<field>``. Invalid prefixes and existing
    destination columns are rejected so evidence cannot be silently overwritten.
    Missing candidate boundaries remain nullable rather than becoming zero or
    false.
    """
    if not isinstance(result_df, pd.DataFrame) or not isinstance(matched_df, pd.DataFrame):
        raise TypeError("ORF and matched-statistics inputs must be pandas DataFrames.")
    if result_df.columns.has_duplicates:
        raise ValueError("ORF table has duplicate column names.")
    if matched_df.columns.has_duplicates:
        raise ValueError("Matched-statistics table has duplicate column names.")
    if method not in {"TIS", "TTS"}:
        raise ValueError("Matched-statistics method must be TIS or TTS.")
    if not isinstance(prefix, str) or _PREFIX_PATTERN.fullmatch(prefix) is None:
        raise ValueError(
            "Matched-statistics prefix must start with a letter or digit and contain "
            "only letters, digits, _, - or ."
        )
    required_orf = {"Genome", "Start", "Stop", "Strand"}
    missing_orf = sorted(required_orf - set(result_df.columns))
    if missing_orf:
        raise ValueError(f"ORF table is missing columns: {', '.join(missing_orf)}.")
    required_matched = {"Genome", "Codon_start", "Strand", "assay", *ATTACH_FIELDS}
    missing_matched = sorted(required_matched - set(matched_df.columns))
    if missing_matched:
        raise ValueError(f"Matched-statistics table is missing columns: {', '.join(missing_matched)}.")
    invalid_assays = [
        value for value in matched_df["assay"]
        if not isinstance(value, str) or value != method
    ]
    if invalid_assays:
        raise ValueError(
            f"Matched-statistics assay does not match method {method!r}: "
            f"{invalid_assays[:3]!r}."
        )

    matched_topology = set(TOPOLOGY_FIELDS).intersection(matched_df.columns)
    result_topology = set(TOPOLOGY_FIELDS).intersection(result_df.columns)
    if matched_topology not in (set(), set(TOPOLOGY_FIELDS)):
        raise ValueError("Matched-statistics table must provide both topology columns.")
    if result_topology not in (set(), set(TOPOLOGY_FIELDS)):
        raise ValueError("ORF table must provide both topology columns.")
    if bool(matched_topology) != bool(result_topology):
        raise ValueError(
            "ORF and matched-statistics tables must either both provide reference topology or both omit it."
        )

    topology_by_genome = {}

    def topology(row, source):
        if not matched_topology:
            return None
        length = getattr(row, "Reference_length")
        circular = getattr(row, "Is_circular")
        if (isinstance(length, bool) or not isinstance(length, Integral) or length < 1
                or not isinstance(circular, (bool, np.bool_))):
            raise ValueError(f"Invalid reference topology in {source} for {row.Genome!r}.")
        value = (int(length), bool(circular))
        previous = topology_by_genome.setdefault(row.Genome, value)
        if previous != value:
            raise ValueError(f"Conflicting reference topology for {row.Genome!r}.")
        return int(length) if circular else None

    destination_columns = [f"{prefix}_{field}" for field in ATTACH_FIELDS]
    collisions = sorted(set(destination_columns) & set(result_df.columns))
    if collisions:
        raise ValueError(f"Matched-statistics output columns already exist: {', '.join(collisions)}.")

    lookup = {}
    for row in matched_df.itertuples(index=False):
        if row.Strand not in {"+", "-"}:
            raise ValueError(f"Invalid strand in matched-statistics table: {row.Strand!r}.")
        if isinstance(row.Codon_start, bool) or not isinstance(row.Codon_start, Integral):
            raise ValueError("Matched-statistics Codon_start values must be integers.")
        if row.Codon_start < 1:
            raise ValueError("Matched-statistics Codon_start values must be positive.")
        if method == "TIS":
            boundary = row.Codon_start + (0 if row.Strand == "+" else 2)
        else:
            boundary = row.Codon_start + (2 if row.Strand == "+" else 0)
        reference_length = topology(row, "matched-statistics table")
        if matched_topology and row.Codon_start > int(row.Reference_length):
            raise ValueError("Matched-statistics Codon_start exceeds its reference length.")
        if reference_length is not None:
            boundary = (boundary - 1) % reference_length + 1
        key = (row.Genome, boundary, row.Strand)
        if key in lookup:
            raise ValueError(f"Duplicate matched-statistics boundary: {key!r}.")
        lookup[key] = row

    output = {field: [] for field in ATTACH_FIELDS}
    for row in result_df.itertuples(index=False):
        if row.Strand not in {"+", "-"}:
            raise ValueError(f"Invalid strand in ORF table: {row.Strand!r}.")
        if (isinstance(row.Start, bool) or not isinstance(row.Start, Integral)
                or isinstance(row.Stop, bool) or not isinstance(row.Stop, Integral)):
            raise ValueError("ORF Start and Stop values must be integers.")
        if row.Start < 1 or row.Stop < row.Start:
            raise ValueError("ORF Start and Stop coordinates are invalid.")
        if method == "TIS":
            boundary = row.Start if row.Strand == "+" else row.Stop
        else:
            boundary = row.Stop if row.Strand == "+" else row.Start
        reference_length = topology(row, "ORF table")
        if result_topology:
            physical_length = int(row.Reference_length)
            if row.Start > physical_length:
                raise ValueError("ORF Start exceeds its reference length.")
            if bool(row.Is_circular):
                if row.Stop - row.Start + 1 > physical_length:
                    raise ValueError("Circular ORF spans more than one revolution.")
            elif row.Stop > physical_length:
                raise ValueError("Linear ORF exceeds its reference length.")
        if reference_length is not None:
            boundary = (boundary - 1) % reference_length + 1
        matched = lookup.get((row.Genome, boundary, row.Strand))
        for field in ATTACH_FIELDS:
            output[field].append(getattr(matched, field) if matched is not None else pd.NA)

    result = result_df.copy()
    for field, values in output.items():
        dtype = (
            "Float64" if field in _FLOAT_RESULT_FIELDS
            else "boolean" if field == "significant"
            else "string"
        )
        result[f"{prefix}_{field}"] = (
            _nullable_integer_array(values)
            if field in _INTEGER_RESULT_FIELDS else pd.array(values, dtype=dtype)
        )
    return result
