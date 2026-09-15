"""Coordinate identifier helpers shared by prediction and output code."""

import re


_INTERVAL_PATTERN = re.compile(r"(-?\d+)-(-?\d+)")


def format_coordinate_id(
    seqid: str,
    start: int,
    stop: int,
    strand: str,
) -> str:
    """Return ``seqid:start-stop:strand`` without restricting the seqid.

    Coordinates may be zero/negative for theoretical codon windows near a
    reference edge. Public ORF identifiers use the same format with positive,
    one-based coordinates.
    """
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid strand in coordinate identifier: {strand!r}")
    return f"{seqid}:{int(start)}-{int(stop)}:{strand}"


def parse_coordinate_id(identifier: str) -> tuple[str, int, int, str]:
    """Parse a coordinate identifier from the right.

    Splitting from the right is essential because GFF3/BAM sequence IDs may
    themselves contain colons and hyphens.
    """
    try:
        seqid, interval, strand = str(identifier).rsplit(":", 2)
        match = _INTERVAL_PATTERN.fullmatch(interval)
        if match is None:
            raise ValueError
        start, stop = map(int, match.groups())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid coordinate identifier: {identifier!r}") from exc
    if not seqid or strand not in {"+", "-"}:
        raise ValueError(f"Invalid coordinate identifier: {identifier!r}")
    return seqid, start, stop, strand


def normalize_circular_interval(start: int, stop: int, reference_length: int) -> tuple[int, int]:
    """Return a one-revolution interval with a canonical zero-based start.

    Circular ORFs use virtual coordinates: the start is on the physical
    reference and the inclusive stop may be greater than or equal to the
    reference length when the interval crosses the origin.
    """
    if reference_length <= 0:
        raise ValueError("Circular reference length must be positive.")
    span = int(stop) - int(start) + 1
    if not 1 <= span <= reference_length:
        raise ValueError("Circular intervals must span between one base and one revolution.")
    normalized_start = int(start) % reference_length
    return normalized_start, normalized_start + span - 1


def split_circular_interval(start: int, stop: int, reference_length: int) -> tuple[tuple[int, int], ...]:
    """Split a zero-based inclusive virtual interval into physical intervals."""
    start, stop = normalize_circular_interval(start, stop, reference_length)
    if stop < reference_length:
        return ((start, stop),)
    return ((start, reference_length - 1), (0, stop % reference_length))


def unique_modulo_positions(start: int, stop: int, reference_length: int) -> tuple[int, ...]:
    """Return physical positions for an inclusive range, without duplicates."""
    if reference_length <= 0:
        raise ValueError("Circular reference length must be positive.")
    if stop < start:
        return ()
    span = min(stop - start + 1, reference_length)
    return tuple((start + offset) % reference_length for offset in range(span))
