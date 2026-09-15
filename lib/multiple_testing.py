"""Shared validation and conservative rounding for adjusted p-values."""

import math

import numpy as np


# Adjustment uses a constant number of arithmetic operations plus a monotone
# reduction whose depth grows with the family. This deliberately small upward
# engineering tolerance covers observed one-to-two-ULP backend undershoots. It
# is not a formal floating-point error proof.
_ROUNDING_GUARD_PER_STAGE = 8 * np.finfo(np.float64).eps


def conservative_adjustment(p_values, correction, backend):
    """Return validated, slightly upward-rounded q-values or ``None``."""
    p_values = np.asarray(p_values, dtype=float)
    try:
        q_values = np.asarray(
            backend(p_values.copy(), method=correction), dtype=float,
        )
    except (ArithmeticError, TypeError, ValueError):
        return None
    if (q_values.shape != p_values.shape
            or not np.all(np.isfinite(q_values))
            or np.any(q_values < 0)
            or np.any(q_values > 1)):
        return None
    # An adjusted p-value cannot be smaller than its raw p-value. Permit only a
    # one-ULP backend undershoot here; a material violation invalidates the
    # complete family rather than silently removing multiplicity correction.
    if np.any(q_values < np.nextafter(p_values, -math.inf)):
        return None
    q_values = np.maximum(q_values, p_values)
    if len(q_values) <= 1:
        return q_values

    reduction_stages = math.ceil(math.log2(len(q_values)))
    relative_guard = _ROUNDING_GUARD_PER_STAGE * (4 + reduction_stages)
    guarded = np.minimum(1.0, q_values * (1.0 + relative_guard))
    movable = (guarded > 0) & (guarded < 1)
    guarded[movable] = np.nextafter(guarded[movable], math.inf)
    return np.maximum(guarded, p_values)
