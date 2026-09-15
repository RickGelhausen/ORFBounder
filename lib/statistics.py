"""Optional, model-conditional local enrichment of raw read endpoints.

The testing family is every eligible codon window in the reference for one
sample and method. It is defined without consulting peak height or ORF calls.
"""

from bisect import bisect_left, bisect_right
from itertools import accumulate
import math
from numbers import Integral, Real

import numpy as np
import pandas as pd
from scipy.stats import binom, false_discovery_control

from lib import misc
from lib.coordinates import format_coordinate_id, parse_coordinate_id, unique_modulo_positions
from lib.multiple_testing import conservative_adjustment


STATISTIC_FIELDS = (
    "peak_count", "background_count", "peak_fold_enrichment",
    "peak_p_value", "peak_q_value", "peak_significant", "peak_correction",
    "peak_tests_in_family", "peak_fit_status",
)

LOCAL_MODEL_NAME = "conditional_binomial_local_endpoint_count"
LOCAL_ALTERNATIVE = "greater"
MIN_POSITIVE_P_VALUE = float(np.nextafter(0.0, 1.0))
# A deliberately conservative engineering limit for scipy's float-backed
# binomial calculation. It is far above realistic local endpoint depth while
# avoiding observed loss of accuracy at much larger, otherwise legal integers.
MAX_SUPPORTED_TOTAL_COUNT = (1 << 31) - 1
# Bound the rare log-tail fallback so an adversarial count cannot consume
# unbounded CPU. A failure remains in the correction family with p=1.
MAX_TAIL_RECURRENCE_WORK = 1_000_000
# Constructing an exact binomial coefficient also has a work cost. Very large
# coefficients fail closed instead of relying on cancellation-prone log-gamma
# calculations in the rare extreme-tail fallback.
MAX_EXACT_COMBINATION_WORK = 100_000
_MIN_NORMAL_P_VALUE = float(np.finfo(np.float64).tiny)
_LOG_MIN_POSITIVE_P_VALUE = math.log(MIN_POSITIVE_P_VALUE)
_DIRECT_TAIL_ROUNDING_GUARD_PER_TRIAL = 64 * np.finfo(np.float64).eps
_SIGNED_INT64_MIN = -(1 << 63)
_SIGNED_INT64_MAX = (1 << 63) - 1

_COLUMNS = (
    "candidate_id", "Genome", "Window_start", "Window_stop", "Strand", "Codon", "Codon_start",
    "peak_count", "background_count", "peak_width", "background_width",
    "fold_enrichment", "p_value", "q_value", "significant", "correction",
    "tests_in_family", "fit_status", "Reference_length", "Is_circular",
)


class _BinomialNumericLimit(RuntimeError):
    """Raised when a safe finite local binomial tail cannot be produced."""


def _adjust_p_values(p_values, correction):
    """Return finite in-range adjusted values, or ``None`` on backend failure."""
    return conservative_adjustment(p_values, correction, false_discovery_control)


def _log_binomial_upper_tail(successes, trials, probability):
    """Evaluate an extreme upper tail in log space with a bounded recurrence.

    scipy's survival function can return zero well before the mathematical tail
    leaves float64 range. Starting at the observed probability and applying the
    successive-PMF ratio avoids that premature underflow. The final omitted
    mass is bounded by a geometric series because the ratio decreases above the
    binomial mode.
    """
    extended = np.longdouble
    if min(successes, trials - successes) > MAX_EXACT_COMBINATION_WORK:
        raise _BinomialNumericLimit("Exact binomial coefficient exceeded its work limit.")
    coefficient = math.comb(trials, successes)
    # Keep the integer coefficient exact until its leading bits are converted
    # to longdouble. scipy's float64 logpmf can lose many tail ULPs through
    # cancellation of large log-gamma terms, even for supported read counts.
    shift = max(0, coefficient.bit_length() - 64)
    leading_bits = coefficient >> shift
    log_first = (
        np.log(extended(leading_bits)) + extended(shift) * np.log(extended(2))
        + extended(successes) * np.log(extended(probability))
        + extended(trials - successes) * np.log1p(-extended(probability))
    )
    if not np.isfinite(log_first):
        raise _BinomialNumericLimit("Could not evaluate the observed binomial mass.")

    odds = extended(probability) / (extended(1.0) - extended(probability))
    relative_term = extended(1.0)
    relative_sum = extended(1.0)
    epsilon = np.finfo(np.longdouble).eps
    current = successes
    work = 0
    while current < trials:
        if work >= MAX_TAIL_RECURRENCE_WORK:
            raise _BinomialNumericLimit("Binomial tail recurrence exceeded its work limit.")
        ratio = (extended(trials - current) / extended(current + 1)) * odds
        # A survival probability small enough to need this fallback must start
        # above the mode, where successive masses decrease monotonically.
        if not np.isfinite(ratio) or ratio < 0 or ratio >= 1:
            raise _BinomialNumericLimit("Binomial tail recurrence is outside its safe range.")
        relative_term *= ratio
        relative_sum += relative_term
        current += 1
        work += 1
        if current == trials or relative_term == 0:
            break

        next_ratio = (extended(trials - current) / extended(current + 1)) * odds
        if not np.isfinite(next_ratio) or next_ratio < 0 or next_ratio >= 1:
            raise _BinomialNumericLimit("Binomial tail recurrence is outside its safe range.")
        remainder_bound = relative_term * next_ratio / (extended(1.0) - next_ratio)
        if remainder_bound <= extended(8.0) * epsilon * relative_sum:
            # Add the recurrence's geometric remainder bound rather than
            # silently dropping the final small terms.
            relative_sum += remainder_bound
            break

    log_tail = log_first + np.log(relative_sum)
    # The coefficient is exact, but logarithms and summation are not. Bound
    # their accumulated rounding at the scale of the input and recurrence,
    # then round the final float upward as before. This is a conservative
    # engineering guard, not a formal interval-arithmetic proof.
    log_tail += extended(64) * epsilon * extended(trials + work + 1)
    if not np.isfinite(log_tail) or log_tail > 0:
        raise _BinomialNumericLimit("Could not evaluate a valid binomial tail.")
    return log_tail


def _binomial_upper_tail(successes, trials, probability):
    """Return ``(p_value, fit_status)`` without ever exporting a zero tail."""
    if trials > MAX_SUPPORTED_TOTAL_COUNT:
        return 1.0, "numerical_limit"
    # The caller's peak-width fraction may round just below its exact rational
    # value. The upper tail is monotone in this probability, so one upward step
    # places the binary64 input on the conservative side of that conversion.
    guarded_probability = min(1.0, math.nextafter(probability, math.inf))
    try:
        direct = float(binom.sf(successes - 1, trials, guarded_probability))
    except (ArithmeticError, TypeError, ValueError):
        direct = math.nan
    if math.isfinite(direct) and _MIN_NORMAL_P_VALUE <= direct <= 1.0:
        # scipy's otherwise accurate direct tail can still undershoot by a few
        # output ULPs. Apply a deliberately small operation-scaled engineering
        # tolerance and round outward. This is not a formal interval bound.
        relative_guard = _DIRECT_TAIL_ROUNDING_GUARD_PER_TRIAL * (trials + 1)
        direct = min(1.0, direct * (1.0 + relative_guard))
        if 0 < direct < 1:
            direct = math.nextafter(direct, math.inf)
        return direct, "ok"

    try:
        log_tail = _log_binomial_upper_tail(successes, trials, guarded_probability)
    except (_BinomialNumericLimit, ArithmeticError, TypeError, ValueError):
        # p=1 is deliberately conservative and keeps this eligible hypothesis
        # in the predeclared multiple-testing family. The standard-library and
        # scipy exceptions cover backend failures before our explicit finite
        # checks can translate them to _BinomialNumericLimit.
        return 1.0, "numerical_limit"
    if log_tail <= _LOG_MIN_POSITIVE_P_VALUE:
        return MIN_POSITIVE_P_VALUE, "p_value_floored"

    value = float(np.exp(log_tail))
    if value <= 0 or not math.isfinite(value):
        return 1.0, "numerical_limit"
    # The log calculation is approximate. Move one representable value upward
    # so conversion cannot make an extreme p-value anti-conservatively small.
    value = float(np.nextafter(value, math.inf))
    return min(value, 1.0), "ok"


def _fold_enrichment(peak_count, background_count, peak_width, background_width):
    """Return the descriptive density ratio and its exact direction around one."""
    if not background_width or not peak_count + background_count:
        return np.nan, 0
    if not background_count:
        return np.inf, 1
    scaled_peak = peak_count * background_width
    scaled_background = background_count * peak_width
    try:
        ratio = scaled_peak / scaled_background
    except OverflowError:
        ratio = np.inf
    if scaled_peak == scaled_background:
        return 1.0, 0
    if scaled_peak > scaled_background:
        return max(ratio, math.nextafter(1.0, math.inf)), 1
    return min(ratio, math.nextafter(1.0, 0.0)), -1


def _count_index(position_counts, genome_length):
    """Index sparse raw counts without allocating an array the genome's size."""
    for position, count in position_counts.items():
        if isinstance(position, bool) or not isinstance(position, Integral):
            raise ValueError("Statistical testing requires integer endpoint positions.")
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 0:
            raise ValueError("Statistical testing requires nonnegative raw integer counts.")
    positions = sorted(pos for pos in position_counts if 0 <= pos < genome_length)
    totals = [0, *accumulate(int(position_counts[pos]) for pos in positions)]

    def count(start, stop):
        return totals[bisect_right(positions, stop)] - totals[bisect_left(positions, start)]

    return count


def _nullable_integer_array(values):
    """Use nullable Int64 when exact; retain larger Python integers as objects."""
    values = list(values)
    if any(
        value is not pd.NA and pd.notna(value)
        and not _SIGNED_INT64_MIN <= int(value) <= _SIGNED_INT64_MAX
        for value in values
    ):
        return pd.array(
            [pd.NA if value is pd.NA or pd.isna(value) else int(value) for value in values],
            dtype=object,
        )
    return pd.array(values, dtype="Int64")


def analyze_peak_enrichment(
    genome_dict,
    raw_position_dict,
    search_codons,
    flank_width=50,
    fdr=0.05,
    correction="by",
    circular_contigs=None,
):
    """Return all candidate windows with a one-sided conditional binomial test.

    Inputs use 0-based positions; output window coordinates are 1-based,
    inclusive. ``raw_position_dict`` must come from a fiveprime/threeprime
    PositionReader *before* normalization. Mode validation belongs to the
    caller because integer counts alone cannot distinguish global coverage.

    Under a uniform local endpoint-rate null, the peak-window count conditional
    on peak plus flank counts is binomial, with probability equal to the peak
    window's fraction of the tested bases. Reference edges clip windows and flanks;
    all eligible codon windows, including those with zero reads, enter FDR
    adjustment jointly across chromosomes and strands for this sample/method.
    """
    if isinstance(flank_width, bool) or not isinstance(flank_width, Integral) or flank_width < 1:
        raise ValueError("Statistical flank width must be a positive integer.")
    if isinstance(fdr, bool) or not isinstance(fdr, Real) or not math.isfinite(fdr) or not 0 < fdr <= 1:
        raise ValueError("Statistical FDR threshold must be greater than 0 and at most 1.")
    if correction not in {"bh", "by"}:
        raise ValueError("Statistical correction must be 'bh' or 'by'.")

    circular_contigs = set(circular_contigs or ())
    unknown_circular = circular_contigs - set(genome_dict)
    if unknown_circular:
        raise ValueError(
            f"Circular sequences are missing from the genome: {', '.join(sorted(unknown_circular))}."
        )
    records = []
    enrichment_directions = []
    for chrom, sequence in genome_dict.items():
        is_circular = chrom in circular_contigs
        _, codons = misc.create_codon_interlaps(
            chrom, sequence.upper(), search_codons, circular=is_circular,
        )
        counters = {
            strand: _count_index(raw_position_dict.get((chrom, strand), {}), len(sequence))
            for strand in ("+", "-")
        }
        for candidate_id, values in codons.items():
            _, theoretical_start, theoretical_stop, strand = parse_coordinate_id(candidate_id)
            if is_circular:
                peak_positions = unique_modulo_positions(
                    theoretical_start, theoretical_stop, len(sequence),
                )
                seen = set(peak_positions)
                background_positions = []
                for flank_start, flank_stop in (
                    (theoretical_start - flank_width, theoretical_start - 1),
                    (theoretical_stop + 1, theoretical_stop + flank_width),
                ):
                    for position in unique_modulo_positions(
                            flank_start, flank_stop, len(sequence)):
                        if position not in seen:
                            seen.add(position)
                            background_positions.append(position)
                position_counts = raw_position_dict.get((chrom, strand), {})
                peak_count = sum(int(position_counts.get(position, 0)) for position in peak_positions)
                background_count = sum(
                    int(position_counts.get(position, 0)) for position in background_positions
                )
                peak_width = len(peak_positions)
                background_width = len(background_positions)
                start, stop = theoretical_start, theoretical_stop
                codon_start = (theoretical_start + (2 if strand == "+" else 0)) % len(sequence)
            else:
                codon_start = theoretical_start + (2 if strand == "+" else 0)
                start, stop = max(0, theoretical_start), min(len(sequence) - 1, theoretical_stop)
                left, right = max(0, start - flank_width), min(len(sequence) - 1, stop + flank_width)
                count = counters[strand]
                peak_count = count(start, stop)
                background_count = count(left, start - 1) + count(stop + 1, right)
                peak_width = stop - start + 1
                background_width = (start - left) + (right - stop)
            total_count = peak_count + background_count
            p_value = 1.0
            fit_status = "ok"
            if total_count and background_width:
                null_probability = peak_width / (peak_width + background_width)
                p_value, fit_status = _binomial_upper_tail(
                    peak_count, total_count, null_probability,
                )
            fold_enrichment, enrichment_direction = _fold_enrichment(
                peak_count, background_count, peak_width, background_width,
            )
            records.append({
                "candidate_id": format_coordinate_id(chrom, start + 1, stop + 1, strand),
                "Genome": chrom, "Window_start": start + 1, "Window_stop": stop + 1,
                "Strand": strand, "Codon": values[0], "Codon_start": codon_start + 1, "peak_count": peak_count,
                "background_count": background_count, "peak_width": peak_width,
                "background_width": background_width, "fold_enrichment": fold_enrichment,
                "p_value": p_value, "fit_status": fit_status,
                "Reference_length": len(sequence), "Is_circular": is_circular,
            })
            enrichment_directions.append(enrichment_direction)

    result = pd.DataFrame.from_records(records, columns=_COLUMNS)
    if not result.empty:
        p_values = result["p_value"].to_numpy(dtype=float)
        q_values = _adjust_p_values(p_values, correction)
        if q_values is None:
            # The entire predeclared family remains present. Discard its
            # adjusted evidence rather than exporting unsafe q-values. Valid
            # raw p-values remain available for audit.
            result["q_value"] = 1.0
            result["significant"] = False
            result.loc[
                result["fit_status"] != "numerical_limit", "fit_status",
            ] = "correction_limit"
        else:
            result["q_value"] = np.maximum(q_values, MIN_POSITIVE_P_VALUE)
            result["significant"] = (
                (result["q_value"] <= fdr)
                & (np.asarray(enrichment_directions) > 0)
                & (result["fit_status"] != "numerical_limit")
            )
        result["correction"] = correction
        result["tests_in_family"] = len(result)
        result = result.sort_values(["Genome", "Window_start", "Window_stop", "Strand"]).reset_index(drop=True)
    return result


def attach_to_orfs(
    result_df, statistics_df, method, sample_name,
    reference_lengths=None, circular_contigs=None,
):
    """Append evidence for an ORF's TIS or TTS boundary without filtering rows.

    Candidate windows match the same strand-dependent offsets used by
    ``create_codon_interlaps`` and ``detect_potential_orfs``. An inferred boundary
    can therefore receive statistical evidence even if this sample did not
    pass its peak-height threshold. Missing candidate boundaries remain NA.
    Existing destination fields are rejected so a previously attached testing
    family cannot be silently replaced.
    """
    if method not in {"TIS", "TTS", "RIBO"}:
        raise ValueError("Statistical method must be TIS, TTS, or RIBO.")
    destination_columns = [f"{sample_name}_{field}" for field in STATISTIC_FIELDS]
    collisions = sorted(set(destination_columns) & set(result_df.columns))
    if collisions:
        raise ValueError(
            "Local-statistics output columns already exist: "
            f"{', '.join(collisions)}."
        )
    reference_lengths = dict(reference_lengths or {})
    circular_contigs = set(circular_contigs or ())

    def topology(row):
        is_circular = bool(getattr(row, "Is_circular", False)) or row.Genome in circular_contigs
        length = getattr(row, "Reference_length", reference_lengths.get(row.Genome))
        if is_circular:
            if pd.isna(length) or int(length) < 1:
                raise ValueError(f"Missing reference length for circular sequence {row.Genome!r}.")
            return int(length)
        return None

    result = result_df.copy()
    lookup = {}
    for row in statistics_df.itertuples(index=False):
        if method in {"TIS", "RIBO"}:
            boundary = row.Codon_start + (0 if row.Strand == "+" else 2)
        else:
            boundary = row.Codon_start + (2 if row.Strand == "+" else 0)
        length = topology(row)
        if length is not None:
            boundary = (boundary - 1) % length + 1
        key = (row.Genome, boundary, row.Strand)
        if key in lookup:
            raise ValueError(f"Duplicate local-statistics boundary: {key!r}.")
        lookup[key] = row

    output = {field: [] for field in STATISTIC_FIELDS}
    for row in result.itertuples(index=False):
        use_start = (method in {"TIS", "RIBO"}) == (row.Strand == "+")
        boundary = row.Start if use_start else row.Stop
        length = topology(row)
        if length is not None:
            boundary = (boundary - 1) % length + 1
        statistics = lookup.get((row.Genome, boundary, row.Strand))
        values = [
            statistics.peak_count, statistics.background_count,
            getattr(statistics, "fold_enrichment", pd.NA), statistics.p_value,
            statistics.q_value, statistics.significant, statistics.correction,
            statistics.tests_in_family, getattr(statistics, "fit_status", pd.NA),
        ] if statistics is not None else [pd.NA] * len(STATISTIC_FIELDS)
        for field, value in zip(STATISTIC_FIELDS, values):
            output[field].append(value)

    for field, values in output.items():
        integer_field = field in {"peak_count", "background_count", "peak_tests_in_family"}
        dtype = "boolean" if field == "peak_significant" else (
            "string" if field in {"peak_correction", "peak_fit_status"} else "Float64"
        )
        result[f"{sample_name}_{field}"] = (
            _nullable_integer_array(values) if integer_field else pd.array(values, dtype=dtype)
        )
    return result
