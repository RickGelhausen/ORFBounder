"""Resolve codon searches and protein translation against one NCBI table."""

from numbers import Integral

from Bio.Data import CodonTable
from Bio.Seq import Seq


def get_genetic_code(code: int = 11) -> CodonTable.CodonTable:
    """Return a DNA table with unambiguous stop behavior suitable for ORF search."""
    if isinstance(code, bool) or not isinstance(code, Integral):
        raise ValueError("genetic_code must be an integer NCBI translation table ID.")
    try:
        table = CodonTable.unambiguous_dna_by_id[int(code)]
    except KeyError as exc:
        supported = ", ".join(str(key) for key in sorted(CodonTable.unambiguous_dna_by_id))
        raise ValueError(f"Unknown genetic_code {code}. NCBI table IDs available: {supported}.") from exc
    dual_codons = set(table.stop_codons) & set(table.forward_table)
    if dual_codons:
        raise ValueError(
            f"genetic_code {code} has context-dependent stop/sense codons "
            f"({', '.join(sorted(dual_codons))}); ORFBounder does not support these tables."
        )
    return table


def _codon_list(values: list[str] | str, kind: str) -> list[str]:
    if isinstance(values, str):
        values = values.split(",")
    try:
        codons = list(dict.fromkeys(value.strip().upper() for value in values))
    except (TypeError, AttributeError) as exc:
        raise ValueError(f"{kind}_codons must contain DNA triplets.") from exc
    if not codons or any(len(value) != 3 or set(value) - set("ACGT") for value in codons):
        raise ValueError(f"{kind}_codons must contain nonempty DNA triplets using A, C, G and T.")
    return codons


def resolve_genetic_code(
    code: int = 11,
    start_codons: list[str] | str | None = None,
    stop_codons: list[str] | str | None = None,
) -> tuple[list[str], list[str]]:
    """Resolve validated start/stop lists, retaining the historical table-11 defaults.

    Other tables use their full NCBI start and stop lists when an override is
    omitted. Overrides select a subset of the codons supported by the table.
    """
    table = get_genetic_code(code)
    if start_codons is None:
        start_codons = ["ATG", "GTG", "TTG"] if code == 11 else table.start_codons
    if stop_codons is None:
        stop_codons = ["TAG", "TAA", "TGA"] if code == 11 else table.stop_codons
    starts = _codon_list(start_codons, "start")
    stops = _codon_list(stop_codons, "stop")
    if set(starts) & set(stops):
        raise ValueError("start_codons and stop_codons must not overlap.")
    for kind, codons, allowed in (("start", starts, table.start_codons), ("stop", stops, table.stop_codons)):
        invalid = set(codons) - set(allowed)
        if invalid:
            raise ValueError(
                f"{kind}_codons {', '.join(sorted(invalid))} are not {kind} codons in "
                f"genetic_code {code}. Choose from {', '.join(allowed)}."
            )
    return starts, stops


def translate_orf(nt_sequence: str, genetic_code: int = 11) -> str:
    """Translate with the selected code and retain the terminal stop marker.

    A complete CDS receives initiator methionine, including valid alternative
    starts. Other sequences retain ordinary translation so result generation can
    reject internal stops instead of silently truncating a candidate.
    """
    table = get_genetic_code(genetic_code)
    protein = str(Seq(nt_sequence).translate(table=int(genetic_code), to_stop=False))
    upper_sequence = nt_sequence.upper()
    if (len(upper_sequence) >= 6 and len(upper_sequence) % 3 == 0
            and upper_sequence[:3] in table.start_codons
            and upper_sequence[-3:] in table.stop_codons
            and protein.endswith("*") and "*" not in protein[:-1]):
        protein = "M" + protein[1:]
    return protein
