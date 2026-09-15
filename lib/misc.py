#!/usr/bin/env python

"""
Miscellaneous functions used throughout the ORFBounder pipeline.
"""

from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import collections
import numpy as np
import pandas as pd

from Bio.Seq import Seq

from interlap import InterLap
import lib.expression as expr
from lib.alignment_reader import validate_normalization_scope
from lib.coordinates import (
    format_coordinate_id, normalize_circular_interval, split_circular_interval,
    unique_modulo_positions,
)
from lib.genetic_code import get_genetic_code, translate_orf



CODON_LENGTH = 3
NT_WINDOW_SIZE = 15
NEAR_ANNOTATED_THRESHOLD = 10
CODON_INTERVAL_OFFSET = 2
CODON_INTERVAL_SIZE = 5


@dataclass(frozen=True)
class AnnotationFeature:
    """One logical CDS, possibly represented by multiple GFF3 rows."""

    feature_id: str
    chromosome: str
    start: int
    stop: int
    strand: str
    gene_id: str
    locus_tag: str
    name: str
    read_list: tuple
    gene_name: str
    old_locus_tag: str
    segments: tuple[tuple[int, int], ...]
    phases: tuple[str, ...]
    exception: str = ""
    part_count: int = 1
    reference_length: int | None = None
    is_circular: bool = False
    crosses_origin: bool = False

    @property
    def display_label(self) -> str:
        """Return a human-facing label; unlike feature identity it need not be unique."""
        return (self.locus_tag or self.gene_name or self.name or self.gene_id
                or format_coordinate_id(self.chromosome, self.start, self.stop, self.strand))

    @property
    def is_complex(self) -> bool:
        """Whether annotation semantics exceed a single ordinary CDS interval."""
        return self.part_count > 1 or "ribosomal slippage" in self.exception.lower()


class AnnotationGene(list):
    """Backward-compatible five-field density record with annotation metadata.

    Existing callers can continue to index or compare this value as
    ``[chromosome, start, stop, strand, density]``. New code can use the
    attributes to keep feature identity distinct from its display label and to
    inspect multipart/phase metadata.
    """

    def __init__(self, feature: AnnotationFeature):
        super().__init__([
            feature.chromosome, feature.start, feature.stop, feature.strand, 0,
        ])
        self.feature_id = feature.feature_id
        self.display_label = feature.display_label
        self.segments = feature.segments
        self.phases = feature.phases
        self.exception = feature.exception
        self.is_complex = feature.is_complex
        self.reference_length = feature.reference_length
        self.is_circular = feature.is_circular
        self.crosses_origin = feature.crosses_origin


def _parse_gff_attributes(attributes: str) -> dict[str, str]:
    """Parse GFF3 attributes while retaining the project's case-tolerant input."""
    parsed = {}
    for field in str(attributes).split(";"):
        if field.strip() in {"", "."}:
            continue
        if "=" not in field:
            raise ValueError(f"Error: invalid gff, wrongly formatted attribute fields.\n{field}")
        name, value = field.split("=", 1)
        normalized_name = name.strip().lower()
        if not normalized_name:
            raise ValueError("Error: GFF3 attribute names must not be empty.")
        if normalized_name in parsed:
            raise ValueError(
                f"Error: duplicate GFF3 attribute name {name.strip()!r}; "
                "use one comma-separated value list where the attribute permits multiple values."
            )
        parsed[normalized_name] = unquote(value.strip())
    return parsed


def _parse_gff_list_attribute(attributes: str, attribute_name: str) -> tuple[str, ...]:
    """Return one GFF3 list attribute without splitting escaped commas.

    Commas separate values in list-valued attributes such as ``Parent`` only
    before percent decoding.  Decoding the complete value first would turn an
    escaped literal comma (``%2C``) into a false list boundary.
    """
    values = []
    for field in str(attributes).split(";"):
        if field.strip() in {"", "."} or "=" not in field:
            continue
        name, value = field.split("=", 1)
        if name.strip().lower() != attribute_name.lower():
            continue
        values.extend(
            unquote(item.strip()) for item in value.split(",") if item.strip()
        )
    return tuple(values)


def parse_reference_topology(annotation_path: Path) -> tuple[dict[str, int], set[str]]:
    """Return declared GFF3 sequence lengths and circular sequence identifiers.

    Circularity is deliberately opt-in and is recognized only from a whole
    ``region`` feature carrying ``Is_circular=true``. Sequence identifiers are
    URL-decoded just like ordinary feature records.
    """
    annotation_path = Path(annotation_path)
    if not annotation_path.is_file():
        raise FileNotFoundError(f"Annotation file does not exist: {annotation_path}")

    reference_lengths: dict[str, int] = {}
    circular_contigs: set[str] = set()
    topology_rows: list[tuple[str, int, int, bool, int]] = []

    def register_length(seqid: str, length: int, context: str) -> None:
        previous = reference_lengths.get(seqid)
        if previous is not None and previous != length:
            raise ValueError(
                f"Conflicting GFF3 reference lengths for {seqid!r}: {previous} and {length} ({context})."
            )
        reference_lengths[seqid] = length

    with annotation_path.open(encoding="utf-8") as annotation_file:
        for line_number, line in enumerate(annotation_file, start=1):
            stripped = line.strip()
            if stripped == "##FASTA":
                break
            if stripped.startswith("##sequence-region"):
                fields = stripped.split()
                if len(fields) != 4:
                    raise ValueError(f"Invalid ##sequence-region directive on line {line_number}.")
                try:
                    start, stop = int(fields[2]), int(fields[3])
                except ValueError as exc:
                    raise ValueError(f"Invalid ##sequence-region coordinates on line {line_number}.") from exc
                if start != 1 or stop < 1:
                    raise ValueError("ORFBounder requires ##sequence-region coordinates to start at 1.")
                register_length(unquote(fields[1]), stop, f"line {line_number}")
                continue
            if not stripped or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9 or fields[2].lower() != "region":
                continue
            attributes = _parse_gff_attributes(fields[8])
            circular_value = attributes.get("is_circular")
            if circular_value is None:
                continue
            if circular_value.lower() not in {"true", "false"}:
                raise ValueError(
                    f"Invalid Is_circular value on GFF3 line {line_number}: {circular_value!r}."
                )
            try:
                start, stop = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid circular region coordinates on line {line_number}.") from exc
            is_circular = circular_value.lower() == "true"
            if is_circular and (start != 1 or stop < 1):
                raise ValueError(
                    f"Circular region for {unquote(fields[0])!r} must span 1 through its reference length."
                )
            seqid = unquote(fields[0])
            if is_circular:
                register_length(seqid, stop, f"line {line_number}")
            topology_rows.append((seqid, start, stop, is_circular, line_number))

    topology_declarations: dict[str, tuple[bool, int]] = {}
    for seqid, start, stop, is_circular, line_number in topology_rows:
        # A false value on a subregion says nothing about reference topology.
        # Once all directives/true declarations have supplied lengths, only a
        # whole-reference row participates in consistency checks.
        if start != 1 or stop != reference_lengths.get(seqid):
            continue
        previous = topology_declarations.get(seqid)
        if previous is not None and previous[0] != is_circular:
            raise ValueError(
                f"Conflicting Is_circular declarations for {seqid!r} "
                f"on GFF3 lines {previous[1]} and {line_number}."
            )
        topology_declarations[seqid] = (is_circular, line_number)

    circular_contigs.update(
        seqid for seqid, (is_circular, _) in topology_declarations.items()
        if is_circular
    )

    return reference_lengths, circular_contigs


def _circular_feature_bounds(group: list[dict], reference_length: int) -> tuple[int, int]:
    """Return the shortest forward virtual arc covering all CDS parts."""
    physical_segments = []
    for row in group:
        physical_segments.extend(split_circular_interval(
            row["start"] - 1, row["stop"] - 1, reference_length,
        ))
    merged = []
    for start, stop in sorted(physical_segments):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(stop, merged[-1][1]))
        else:
            merged.append((start, stop))

    largest_gap = -1
    arc_start = 0
    for index, (_, stop) in enumerate(merged):
        next_start = merged[(index + 1) % len(merged)][0]
        if index == len(merged) - 1:
            next_start += reference_length
        gap = next_start - stop - 1
        candidate_start = next_start % reference_length
        if gap > largest_gap or (gap == largest_gap and candidate_start < arc_start):
            largest_gap, arc_start = gap, candidate_start
    span = reference_length - largest_gap
    return arc_start, arc_start + span - 1


def parse_annotation_features(
    annotation_path: Path,
    *,
    reference_topology: tuple[dict[str, int], set[str]] | None = None,
) -> list[AnnotationFeature]:
    """Parse logical CDS features from GFF3.

    GFF3 permits a discontinuous feature to occupy multiple rows with the same
    ``ID``. Those rows are grouped here. ``locus_tag`` and ``Name`` remain
    display labels and are deliberately not used as unique identities.
    """
    annotation_path = Path(annotation_path)
    if not annotation_path.is_file():
        raise FileNotFoundError(f"Annotation file does not exist: {annotation_path}")

    if reference_topology is None:
        reference_lengths, circular_contigs = parse_reference_topology(annotation_path)
    else:
        reference_lengths, circular_contigs = reference_topology
        reference_lengths = dict(reference_lengths)
        circular_contigs = set(circular_contigs)
    rows = []
    feature_nodes: dict[str, list[dict]] = {}
    with annotation_path.open(encoding="utf-8") as annotation_file:
        for line_number, line in enumerate(annotation_file, start=1):
            if line.strip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    "Error: invalid GFF; feature records require at least nine tab-separated columns."
                )
            chromosome = unquote(fields[0])
            feature_type = fields[2].lower()
            attributes = _parse_gff_attributes(fields[8])
            parents = _parse_gff_list_attribute(fields[8], "parent")
            if attributes.get("id"):
                feature_nodes.setdefault(attributes["id"], []).append({
                    "feature_type": feature_type,
                    "attributes": attributes,
                    "parents": parents,
                })
            if feature_type not in {"cds", "gene", "pseudogene"}:
                continue
            try:
                start, stop = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Error: invalid GFF coordinates on line {line_number}.") from exc
            if start < 1 or stop < start:
                raise ValueError(f"Error: invalid GFF coordinates on line {line_number}.")
            if chromosome in reference_lengths and stop > reference_lengths[chromosome]:
                topology = "Circular" if chromosome in circular_contigs else "Declared-reference"
                raise ValueError(
                    f"{topology} GFF3 feature rows must use physical coordinates within "
                    f"1..{reference_lengths[chromosome]} (line {line_number})."
                )
            rows.append({
                "line_number": line_number,
                "chromosome": chromosome,
                "feature_type": feature_type,
                "start": start,
                "stop": stop,
                "strand": fields[6],
                "phase": fields[7],
                "attributes": attributes,
                "parents": parents,
                "read_list": tuple(fields[9:]),
            })

    gene_by_coordinates = {}

    def annotation_metadata(node: dict) -> dict[str, str]:
        attributes = node["attributes"]
        gene_identifier = attributes.get("gene_id", "")
        if not gene_identifier and node["feature_type"] in {"gene", "pseudogene"}:
            gene_identifier = attributes.get("id", "")
        return {
            "gene_id": gene_identifier,
            "name": attributes.get("name", attributes.get("gene_name", "")),
            "locus_tag": attributes.get("locus_tag", attributes.get("gene_id", "")),
            "old_locus_tag": attributes.get("old_locus_tag", ""),
        }

    for row in rows:
        if row["feature_type"] not in {"gene", "pseudogene"}:
            continue
        metadata = annotation_metadata(row)
        coordinates = (row["chromosome"], row["start"], row["stop"], row["strand"])
        gene_by_coordinates[coordinates] = metadata

    def inherited_metadata(parent_ids: tuple[str, ...]) -> dict[str, str]:
        """Resolve the first gene ancestor, following ordinary transcript nodes."""
        pending = deque(parent_ids)
        visited = set()
        fallback: dict[str, str] = {}
        while pending:
            parent_id = pending.popleft()
            if parent_id in visited:
                continue
            visited.add(parent_id)
            for node in feature_nodes.get(parent_id, ()):
                metadata = annotation_metadata(node)
                for name, value in metadata.items():
                    if value:
                        fallback.setdefault(name, value)
                if node["feature_type"] in {"gene", "pseudogene"}:
                    return {
                        name: value or fallback.get(name, "")
                        for name, value in metadata.items()
                    }
                pending.extend(node["parents"])
        return fallback

    cds_groups = {}
    for row in rows:
        if row["feature_type"] != "cds":
            continue
        gff_id = row["attributes"].get("id")
        group_key = ("id", gff_id) if gff_id else ("line", row["line_number"])
        group = cds_groups.setdefault(group_key, [])
        if group and (group[0]["chromosome"], group[0]["strand"]) != (row["chromosome"], row["strand"]):
            raise ValueError(
                f"Error: multipart CDS ID {gff_id!r} occurs on different sequences or strands."
            )
        group.append(row)

    features = []
    for group_key, group in cds_groups.items():
        strand = group[0]["strand"]
        chromosome = group[0]["chromosome"]
        reference_length = reference_lengths.get(chromosome)
        is_circular = chromosome in circular_contigs
        if is_circular:
            start0, stop0 = _circular_feature_bounds(group, reference_length)
            start, stop = start0 + 1, stop0 + 1
            ordered = sorted(
                group,
                key=lambda row: ((row["start"] - 1) % reference_length - start0) % reference_length,
                reverse=strand == "-",
            )
        else:
            ordered = sorted(group, key=lambda row: row["start"], reverse=strand == "-")
            start = min(row["start"] for row in group)
            stop = max(row["stop"] for row in group)

        def first_attribute(name: str, default: str = "") -> str:
            return next((row["attributes"].get(name, "") for row in ordered
                         if row["attributes"].get(name, "")), default)

        parent_ids = tuple(dict.fromkeys(
            parent_id for row in ordered for parent_id in row["parents"]
        ))
        parent = inherited_metadata(parent_ids)
        if not any(parent.values()):
            parent = gene_by_coordinates.get((chromosome, start, stop, strand), {})

        gff_id = group_key[1] if group_key[0] == "id" else ""
        gene_id = first_attribute("gene_id") or parent.get("gene_id", "") or gff_id
        locus_tag = first_attribute("locus_tag") or parent.get("locus_tag", "")
        name = first_attribute("name") or first_attribute("gene_name")
        gene_name = parent.get("name", "")
        old_locus_tag = first_attribute("old_locus_tag") or parent.get("old_locus_tag", "")
        exceptions = list(dict.fromkeys(
            row["attributes"].get("exception", "") for row in ordered
            if row["attributes"].get("exception", "")
        ))
        feature_id = gff_id or f"gff-row-{group[0]['line_number']}"
        segments = []
        for row in ordered:
            if is_circular:
                physical = split_circular_interval(
                    row["start"] - 1, row["stop"] - 1, reference_length,
                )
                segments.extend((segment_start + 1, segment_stop + 1)
                                for segment_start, segment_stop in physical)
            else:
                segments.append((row["start"], row["stop"]))
        features.append(AnnotationFeature(
            feature_id=feature_id,
            chromosome=chromosome,
            start=start,
            stop=stop,
            strand=strand,
            gene_id=gene_id,
            locus_tag=locus_tag,
            name=name,
            read_list=tuple(value for row in ordered for value in row["read_list"]),
            gene_name=gene_name,
            old_locus_tag=old_locus_tag,
            segments=tuple(dict.fromkeys(segments)),
            phases=tuple(row["phase"] for row in ordered),
            exception=", ".join(exceptions),
            part_count=len(group),
            reference_length=reference_length,
            is_circular=is_circular,
            crosses_origin=is_circular and stop > reference_length,
        ))
    return features


def generate_annotation_dict(annotation_path: Path) -> dict[str, tuple]:
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """

    annotation_dict = {}
    for feature in parse_annotation_features(annotation_path):
        key = format_coordinate_id(feature.chromosome, feature.start, feature.stop, feature.strand)
        annotation_dict[key] = (
            feature.gene_id, feature.locus_tag, feature.name, list(feature.read_list),
            feature.gene_name, feature.old_locus_tag,
        )
    return annotation_dict

def annotation_interlap(
    annotation_file: Path,
    *,
    annotation_features: tuple[AnnotationFeature, ...] | None = None,
) -> tuple[dict[tuple[str, str], InterLap], dict[str, list]]:
    """
    create an interlap object for the annotation
    """
    annotation_interlap_dict = {}
    gene_density_dict = {}
    features = (
        parse_annotation_features(annotation_file)
        if annotation_features is None else annotation_features
    )
    label_counts = Counter(feature.display_label for feature in features)
    reserved_labels = set(label_counts)
    used_keys = set()

    for feature in features:
        label = feature.display_label
        if label_counts[label] == 1:
            feature_key = label
        else:
            coordinate = format_coordinate_id(
                feature.chromosome, feature.start, feature.stop, feature.strand)
            base = f"{label} [{feature.feature_id};{coordinate}]"
            feature_key = base
            suffix = 2
            while feature_key in reserved_labels or feature_key in used_keys:
                feature_key = f"{base}#{suffix}"
                suffix += 1
        used_keys.add(feature_key)
        gene_density_dict[feature_key] = AnnotationGene(feature)

        interlap = annotation_interlap_dict.setdefault(
            (feature.chromosome, feature.strand), InterLap())
        for start, stop in feature.segments:
            interlap.add((start - 1, stop - 1, feature_key))

    return annotation_interlap_dict, gene_density_dict

def calculate_density(
    alignment_position_dict: dict[tuple[str, str], dict[int, int]],
    annotation_interlap_dict: dict[tuple[str, str], InterLap],
    gene_density_dict: dict[str, list]
) -> dict[str, list]:
    """
    calculate density for each annotated gene
    """
    for (chrom, strand) in alignment_position_dict:
        if (chrom, strand) not in annotation_interlap_dict:
            continue
        for position, read_count in alignment_position_dict[(chrom, strand)].items():
            matching_genes = list(annotation_interlap_dict[(chrom, strand)].find((position, position)))
            # Multipart annotations may contain overlapping rows (notably
            # ribosomal-slippage annotations). Count one genomic observation
            # once per logical feature, not once per overlapping part.
            for feature_key in dict.fromkeys(gene[2] for gene in matching_genes):
                gene_density_dict[feature_key][4] += read_count

    return gene_density_dict


def create_codon_interlaps(
    chrom: str,
    genome_seq: str,
    codons: list[str],
    circular: bool = False,
) -> tuple[dict[tuple[str, str], InterLap], dict[str, list]]:
    """
    create interlaps around each codon
    """
    reverse_codons = [str(Seq(codon).reverse_complement()) for codon in codons]
    interlap_dict = {}
    codon_dict = {}
    reference_length = len(genome_seq)
    positions = range(reference_length if circular else max(0, reference_length - CODON_LENGTH + 1))
    for pos in positions:
        codon = ("".join(genome_seq[index % reference_length]
                         for index in range(pos, pos + CODON_LENGTH))
                 if circular else genome_seq[pos:pos + CODON_LENGTH])
        for strand, accepted_codons in (("+", codons), ("-", reverse_codons)):
            if codon not in accepted_codons:
                continue
            interval_start = pos - CODON_INTERVAL_OFFSET if strand == "+" else pos
            interval_stop = interval_start + CODON_INTERVAL_SIZE - 1
            if circular:
                interval_start %= reference_length
                interval_stop = interval_start + min(CODON_INTERVAL_SIZE, reference_length) - 1
            # Preserve the theoretical window in the key so callers can recover
            # the actual codon coordinate; clip only the searchable interval.
            key = format_coordinate_id(chrom, interval_start, interval_stop, strand)
            if circular:
                physical_positions = sorted(set(unique_modulo_positions(
                    interval_start, interval_stop, reference_length,
                )))
                runs = []
                for physical_position in physical_positions:
                    if runs and physical_position == runs[-1][1] + 1:
                        runs[-1] = (runs[-1][0], physical_position)
                    else:
                        runs.append((physical_position, physical_position))
                interlap = interlap_dict.setdefault((chrom, strand), InterLap())
                for start, stop in runs:
                    interlap.add((start, stop, key))
            else:
                interval = (max(0, interval_start), min(reference_length - 1, interval_stop), key)
                interlap_dict.setdefault((chrom, strand), InterLap()).add(interval)
            codon_dict[key] = [codon if strand == "+" else str(Seq(codon).reverse_complement()), 0]
    return interlap_dict, codon_dict


def get_frame(position: int) -> int:
    """
    get the reading frame for the given position
    """
    return position % CODON_LENGTH


def _candidate_coordinate_copies(
    start: int,
    stop: int,
    gene_record: list,
) -> tuple[tuple[int, int], ...]:
    """Return candidate copies that can be compared on a circular reference.

    Circular result and annotation intervals each use a canonical virtual start,
    so two physically overlapping intervals can land one reference length apart
    on that virtual axis. The adjacent copies are sufficient because both
    intervals are bounded to one revolution. Ordinary/legacy annotation
    records retain the historical single-coordinate comparison.
    """
    if not getattr(gene_record, "is_circular", False):
        return ((start, stop),)
    reference_length = getattr(gene_record, "reference_length", None)
    if not isinstance(reference_length, (int, np.integer)) or reference_length < 1:
        return ((start, stop),)
    return tuple(
        (start + shift * int(reference_length), stop + shift * int(reference_length))
        for shift in (-1, 0, 1)
    )


def _candidate_overlaps_annotation_segments(
    start: int,
    stop: int,
    gene_record: list,
) -> bool:
    """Return whether a candidate touches an actual CDS segment.

    A multipart CDS can contain a genomic gap between its outer bounds.  Those
    bounds remain useful for describing the logical feature, but a candidate
    wholly inside the gap must not be classified as an internal CDS.  Legacy
    annotation records have no segment metadata and retain their historical
    bounding-interval behavior.
    """
    annotation_segments = getattr(gene_record, "segments", ())
    if not annotation_segments:
        return True

    candidate_segments = ((start, stop),)
    if getattr(gene_record, "is_circular", False):
        reference_length = getattr(gene_record, "reference_length", None)
        if isinstance(reference_length, (int, np.integer)) and reference_length > 0:
            try:
                candidate_segments = tuple(
                    (segment_start + 1, segment_stop + 1)
                    for segment_start, segment_stop in split_circular_interval(
                        start - 1, stop - 1, int(reference_length),
                    )
                )
            except ValueError:
                return False

    return any(
        candidate_start <= annotation_stop and annotation_start <= candidate_stop
        for candidate_start, candidate_stop in candidate_segments
        for annotation_start, annotation_stop in annotation_segments
    )


def _align_candidate_to_gene(
    start: int,
    stop: int,
    gene_record: list,
) -> tuple[int, int]:
    """Choose the circular candidate copy sharing the gene's local axis."""
    gene_start, gene_stop, strand = gene_record[1:4]

    def score(candidate: tuple[int, int]) -> tuple[int, int, int, int]:
        candidate_start, candidate_stop = candidate
        # Near-annotated and extension relationships share the translation-stop
        # boundary, so prefer that exact alignment before overlap/midpoint ties.
        shared_stop = (
            candidate_stop == gene_stop if strand == "+"
            else candidate_start == gene_start
        )
        overlap = max(
            0,
            min(candidate_stop, gene_stop) - max(candidate_start, gene_start) + 1,
        )
        midpoint_distance = abs(
            (candidate_start + candidate_stop) - (gene_start + gene_stop)
        )
        return (not shared_stop, -overlap, midpoint_distance, candidate_start)

    return min(_candidate_coordinate_copies(start, stop, gene_record), key=score)


def get_genome_information(
    start: int,
    stop: int,
    strand: str,
    genome_seq: str,
    genetic_code: int = 11,
    circular: bool = False,
) -> tuple[str, str, str, str, str]:
    """
    retrieve information from genome including nucleotide sequence, start_codon, stop_codon, amino acid sequence, 15nt window
    """
    if circular:
        positions = unique_modulo_positions(start, stop, len(genome_seq))
        genomic_sequence = "".join(genome_seq[position] for position in positions)
        window_size = min(NT_WINDOW_SIZE, len(genome_seq))
        if strand == "+":
            nt_seq = genomic_sequence
            nt_window = "".join(
                genome_seq[position]
                for position in unique_modulo_positions(
                    start - window_size, start - 1, len(genome_seq),
                )
            )
        else:
            nt_seq = str(Seq(genomic_sequence).reverse_complement())
            downstream = "".join(
                genome_seq[position]
                for position in unique_modulo_positions(
                    stop + 1, stop + window_size, len(genome_seq),
                )
            )
            nt_window = str(Seq(downstream).reverse_complement())
    elif strand == "+":
        nt_seq = genome_seq[start:stop + 1]
        nt_window = genome_seq[max(0, start - NT_WINDOW_SIZE):start]
    else:
        nt_seq = str(Seq(genome_seq[start:stop + 1]).reverse_complement())
        nt_window = str(Seq(genome_seq[stop + 1:stop + 1 + NT_WINDOW_SIZE]).reverse_complement())

    aa_seq = translate_orf(nt_seq, genetic_code)
    start_codon, stop_codon = nt_seq[:CODON_LENGTH], nt_seq[-CODON_LENGTH:]

    return nt_seq, aa_seq, nt_window, start_codon, stop_codon


def get_gene_information(
    chrom: str,
    start_position: int,
    stop_position: int,
    strand: str,
    gene_dict: dict[str, tuple[str, int, int, str, int]]
) -> tuple[str, str]:
    """
    determine the gene_name and gene_type
    """
    # ensure that the positions are 1-based
    start_position += 1
    stop_position += 1

    label = "Unannotated"
    gene_name_assigned = None

    for gene_name, gene_record in gene_dict.items():
        gene_chrom, gene_start, gene_stop, gene_strand, _ = gene_record
        if gene_chrom != chrom or gene_strand != strand:
            continue

        if not _candidate_overlaps_annotation_segments(
                start_position, stop_position, gene_record):
            continue

        candidate_copies = _candidate_coordinate_copies(
            start_position, stop_position, gene_record,
        )

        # Annotated - exact match (highest priority)
        if any(gene_start == candidate_start and gene_stop == candidate_stop
               for candidate_start, candidate_stop in candidate_copies):
            category = "Annotated_Complex" if getattr(gene_record, "is_complex", False) else "Annotated"
            return category, gene_name

        # Near-Annotated (second highest priority)
        if strand == "+":
            if any(abs(candidate_start - gene_start) < NEAR_ANNOTATED_THRESHOLD
                   and candidate_stop == gene_stop
                   for candidate_start, candidate_stop in candidate_copies):
                label = "Near_Annotated"
                gene_name_assigned = gene_name
                continue
        else:  # strand == "-"
            if any(abs(candidate_stop - gene_stop) < NEAR_ANNOTATED_THRESHOLD
                   and candidate_start == gene_start
                   for candidate_start, candidate_stop in candidate_copies):
                label = "Near_Annotated"
                gene_name_assigned = gene_name
                continue

        if label != "Near_Annotated":
            # Internal_Inframe (third priority)
            if strand == "+":
                if any(candidate_start > gene_start and candidate_stop <= gene_stop
                       and candidate_start % CODON_LENGTH == gene_start % CODON_LENGTH
                       for candidate_start, candidate_stop in candidate_copies):
                    label = "Internal_Inframe"
                    gene_name_assigned = gene_name
                    continue
            else:  # strand == "-"
                if any(candidate_stop < gene_stop and candidate_start >= gene_start
                       and candidate_stop % CODON_LENGTH == gene_stop % CODON_LENGTH
                       for candidate_start, candidate_stop in candidate_copies):
                    label = "Internal_Inframe"
                    gene_name_assigned = gene_name
                    continue

            # N-terminal extension (fourth priority)
            if strand == "+":
                if any(candidate_stop == gene_stop and candidate_start < gene_start
                       for candidate_start, candidate_stop in candidate_copies):
                    label = "N-terminal_extension"
                    gene_name_assigned = gene_name
            else:  # strand == "-"
                if any(candidate_start == gene_start and candidate_stop > gene_stop
                       for candidate_start, candidate_stop in candidate_copies):
                    label = "N-terminal_extension"
                    gene_name_assigned = gene_name

        # Internal-OutofFrame (lowest priority among overlapping categories)
        if any(candidate_start >= gene_start and candidate_stop <= gene_stop
               and (
                   candidate_start % CODON_LENGTH != gene_start % CODON_LENGTH
                   if strand == "+"
                   else candidate_stop % CODON_LENGTH != gene_stop % CODON_LENGTH
               )
               for candidate_start, candidate_stop in candidate_copies):
            # Only set if no higher priority label already assigned
            if label not in ["Near_Annotated", "Internal_Inframe", "N-terminal_extension"]:
                label = "Internal_OutofFrame"
                gene_name_assigned = gene_name

    return label, gene_name_assigned if gene_name_assigned else format_coordinate_id(
        chrom, start_position, stop_position, strand)


def annotation_display_label(gene_key: str, gene_dict: dict[str, list]) -> str:
    """Resolve an internal annotation key to its original human-facing label."""
    record = gene_dict.get(gene_key)
    return getattr(record, "display_label", gene_key) if record is not None else gene_key


def _build_annotation_candidate_index(gene_dict: dict[str, list]) -> dict:
    """Index genes that could have a supported relationship to an ORF.

    Every classification in :func:`get_gene_information` requires the ORF and
    annotation intervals to overlap (exact, near, internal, or extension).
    Circular annotations are indexed on adjacent virtual copies because result
    and annotation intervals can use virtual axes one revolution apart.
    """
    intervals = {}
    for ordinal, (gene_key, gene_record) in enumerate(gene_dict.items()):
        chrom, gene_start, gene_stop, strand = gene_record[:4]
        shifts = (0,)
        if getattr(gene_record, "is_circular", False):
            reference_length = getattr(gene_record, "reference_length", None)
            if isinstance(reference_length, (int, np.integer)) and reference_length > 0:
                shifts = (-int(reference_length), 0, int(reference_length))
        records = intervals.setdefault((chrom, strand), [])
        records.extend(
            (gene_start + shift, gene_stop + shift, ordinal, gene_key)
            for shift in shifts
        )
    return {key: InterLap(records) for key, records in intervals.items()}


def _overlapping_annotation_genes(
    annotation_index: dict,
    gene_dict: dict[str, list],
    chrom: str,
    start_position: int,
    stop_position: int,
    strand: str,
) -> dict[str, list]:
    """Return possible annotation matches in original dictionary order."""
    index = annotation_index.get((chrom, strand))
    if index is None:
        return {}
    matches = index.find((start_position + 1, stop_position + 1))
    ordered_keys = sorted({record[3]: record[2] for record in matches}.items(),
                          key=lambda item: item[1])
    return {gene_key: gene_dict[gene_key] for gene_key, _ in ordered_keys}


def calculate_utr_distance(
    start_position: int,
    stop_position: int,
    gene_name: str,
    gene_dict: dict[str, tuple[str, int, int, str, int]]
) -> tuple[int | float, int | float]:
    """
    Distances from the predicted translation start to annotated 5' and 3' ends.

    Inputs are zero-based inclusive ORF coordinates; gene metadata is one-based.
    These describe the start site's location, not inferred untranslated lengths.
    """
    start_position += 1
    stop_position += 1

    if gene_name in gene_dict:
        gene_record = gene_dict[gene_name]
        _, gene_start, gene_stop, strand, _ = gene_record
        start_position, stop_position = _align_candidate_to_gene(
            start_position, stop_position, gene_record,
        )
        if strand == "+":
            fiveprime_dist = start_position - gene_start
            threeprime_dist = gene_stop - start_position
        else:
            fiveprime_dist = gene_stop - stop_position
            threeprime_dist = stop_position - gene_start
    else:
        fiveprime_dist, threeprime_dist = np.nan, np.nan

    return fiveprime_dist, threeprime_dist


def calculate_relative_density(
    rpm: float,
    gene_name: str,
    gene_type: str,
    gene_dict: dict[str, tuple[str, int, int, str, int]]
) -> float:
    """
    calculate the relative density of the signal to the overall gene density
    """
    if gene_dict == {} or gene_type == "N-terminal_extension" or gene_name not in gene_dict:
        return np.nan

    gene_rpm = gene_dict[gene_name][4]
    if gene_rpm != 0:
        density = rpm / gene_rpm
    else:
        density = 0.0

    return density

def generate_result_dataframe(
    detected_orfs_dict: dict,
    gene_tis_dict: dict,
    gene_tts_dict: dict,
    gene_ribo_dict: dict,
    genome: dict[str, str],
    read_count_dict: dict,
    accepted_read_list: list[dict],
    wildcards: list[str],
    headers: tuple[str, str, str],
    genetic_code: int = 11,
    normalization_scope: str = "contig",
    circular_contigs: set[str] | None = None,
) -> pd.DataFrame:
    """
    Generate the final dataframe to be written to file.
    This contains RPKM, TE, nucleotide and aminoacid sequences and more.
    """
    validate_normalization_scope(normalization_scope)
    get_genetic_code(genetic_code)
    circular_contigs = set(circular_contigs or ())
    tis_header, tts_header, ribo_header = headers
    annotation_gene_dict = gene_tis_dict if gene_tis_dict else gene_tts_dict
    annotation_index = _build_annotation_candidate_index(annotation_gene_dict)

    te_header = expr.get_te_header(wildcards)

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count",
              f"{tis_header}_peak_height", f"{tts_header}_peak_height", f"{ribo_header}_peak_height",
              "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq",
              f"{tis_header}_relative_density", f"{tts_header}_relative_density", f"{ribo_header}_relative_density",
              "5'-distance", "3'-distance"] + [f"{card}_rpkm" for card in wildcards] + \
              [f"{cond}_TE" for cond in te_header] + ["Reference_length", "Is_circular"]

    name_list = [f"s{x}" for x in range(len(header))]
    n_tuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for (chrom, strand) in detected_orfs_dict.keys():
        for (raw_start, raw_stop), peak_values in detected_orfs_dict[(chrom, strand)].items():
            start, stop = raw_start, raw_stop
            is_circular = chrom in circular_contigs
            if is_circular:
                start, stop = normalize_circular_interval(start, stop, len(genome[chrom]))
            rpm_tis, rpm_tts, rpm_ribo = peak_values

            if np.isnan(rpm_tis) and np.isnan(rpm_tts):
                continue

            candidate_gene_dict = _overlapping_annotation_genes(
                annotation_index, annotation_gene_dict,
                chrom, start, stop, strand,
            )
            gene_type, gene_key = get_gene_information(
                chrom, start, stop, strand, candidate_gene_dict)
            gene_name = annotation_display_label(gene_key, annotation_gene_dict)

            nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(
                start, stop, strand, genome[chrom], genetic_code, circular=is_circular)

            if aa_seq.count("*") > 1:
                continue

            rpkm_list = []
            te_rate_list = []
            te_list = []
            if read_count_dict != {}:
                for idx, val in enumerate(read_count_dict[(chrom, start, stop, strand)]):
                    total_mapped = expr.normalization_read_count(accepted_read_list[idx], chrom, normalization_scope)
                    rpkm_list.append(expr.calculate_rpkm(total_mapped, val, len(nt_seq)))
                    # TE is a ratio of matched RPKMs for this same ORF. The
                    # shared length and 1e9 factors cancel, so use the exact
                    # pre-rounding rates and reserve two-decimal RPKM values
                    # for display only.
                    te_rate_list.append(val / total_mapped if total_mapped else 0.0)

                te_list = expr.calculate_te(te_rate_list, wildcards)

            identifier = format_coordinate_id(chrom, start + 1, stop + 1, strand)
            codon_count = int(len(nt_seq) / CODON_LENGTH)

            if gene_tis_dict != {}:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_key, gene_tis_dict)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_key, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_key, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_key, gene_type, gene_ribo_dict)
            else:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_key, gene_tts_dict)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_key, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_key, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_key, gene_type, gene_ribo_dict)

            result = [gene_type, identifier, chrom, start + 1, stop + 1, strand, gene_name, codon_count,
                      rpm_tis, rpm_tts, rpm_ribo, start_codon, stop_codon,
                      nt_window, nt_seq, aa_seq, relative_density_tis, relative_density_tts, relative_density_ribo,
                      fiveprime_dist, threeprime_dist] + rpkm_list + te_list + [len(genome[chrom]), is_circular]

            result_rows.append(n_tuple(*result))

    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])
    df_results = df_results.sort_values(by=["Genome", "Start", "Stop", "Strand"])
    df_results.attrs["reference_lengths"] = {
        chrom: len(sequence) for chrom, sequence in genome.items()
    }
    df_results.attrs["circular_contigs"] = sorted(circular_contigs)

    return df_results


def dictionary_depth(dic: dict) -> int:
    """
    get the depth of a dictionary
    """
    queue = deque([(id(dic), dic, 0)])
    already_visited = set()
    while queue:
        id_, o, level = queue.popleft()
        if id_ in already_visited:
            continue
        already_visited.add(id_)
        if isinstance(o, dict):
            queue += ((id(v), v, level + 1) for v in o.values())
    return level
