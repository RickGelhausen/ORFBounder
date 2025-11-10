#!/usr/bin/env python

"""
Miscellaneous functions used throughout the ORFBounder pipeline.
"""

from collections import deque
from pathlib import Path

import re
import collections
import numpy as np
import pandas as pd

from Bio.Seq import Seq

from interlap import InterLap
import lib.expression as expr
import lib.messaging as msg


CODON_LENGTH = 3
NT_WINDOW_SIZE = 15
NEAR_ANNOTATED_THRESHOLD = 10
CODON_INTERVAL_OFFSET = 2
CODON_INTERVAL_SIZE = 5

def generate_annotation_dict(annotation_path: Path) -> dict[str, tuple]:
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)

    ### MAJOR CHANGE 2: Added file validation and pathlib
    # Original code didn't check if annotation file exists
    # Now validates file exists before processing
    """

    if not annotation_path.is_file():
        raise FileNotFoundError(
            msg.error(f"Annotation file does not exist: {annotation_path}")
        )

    annotation_df = pd.read_csv(annotation_path, sep="\t", comment="#", header=None)
    annotation_dict = {}
    gene_dict = {}
    cds_dict = {}

    for row in annotation_df.itertuples(index=False, name='Pandas'):
        chromosome = getattr(row, "_0")
        feature = getattr(row, "_2")
        start = getattr(row, "_3")
        stop = getattr(row, "_4")
        strand = getattr(row, "_6")
        attributes = getattr(row, "_8")
        read_list = [getattr(row, f"_{x}") for x in range(9, len(row))]

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i, _ in enumerate(attribute_list):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            raise ValueError(msg.error(f"Error: invalid gff, wrongly formatted attribute fields.\n{attribute_list}"))

        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag") + 1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag") + 1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name") + 1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name") + 1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id") + 1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id") + 1]

            new_key = f"{chromosome}:{start}-{stop}:{strand}"
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list, old_locus_tag)

        elif feature.lower() in ["gene", "pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name") + 1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name") + 1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag") + 1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id") + 1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag") + 1]

            new_key = f"{chromosome}:{start}-{stop}:{strand}"
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)

    for key, value in cds_dict.items():
        gene_name = ""
        gene_id, locus_tag, name, read_list, old_locus_tag = value

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_dict[key] = (gene_id, locus_tag, name, read_list, gene_name, old_locus_tag)

    return annotation_dict

def annotation_interlap(
    annotation_file: Path
) -> tuple[dict[tuple[str, str], InterLap], dict[str, list]]:
    """
    create an interlap object for the annotation
    """
    annotation_dict = generate_annotation_dict(annotation_file)

    annotation_interlap_dict = {}
    gene_density_dict = {}

    for key, val in annotation_dict.items():
        chromosome, mid, strand = key.split(":")
        start, stop = mid.split("-")
        start, stop = int(start), int(stop)

        if val[1] != "":
            locus_tag = val[1]
        else:
            locus_tag = val[4]

        if (chromosome, strand) in annotation_interlap_dict:
            annotation_interlap_dict[(chromosome, strand)].add((start, stop, locus_tag))
        else:
            inter = InterLap()
            inter.add((start, stop, locus_tag))
            annotation_interlap_dict[(chromosome, strand)] = inter

        if locus_tag not in gene_density_dict:
            gene_density_dict[locus_tag] = [chromosome, start, stop, strand, 0]

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
            for gene in matching_genes:
                gene_density_dict[gene[2]][4] += read_count

    return gene_density_dict


def create_codon_interlaps(
    chrom: str,
    genome_seq: str,
    codons: list[str]
) -> tuple[dict[tuple[str, str], InterLap], dict[str, list]]:
    """
    create interlaps around each codon
    """
    reverse_codons = [str(Seq(codon).reverse_complement()) for codon in codons]

    interlap_dict = {}
    codon_dict = {}

    for pos in range(len(genome_seq) - CODON_LENGTH + 1):
        codon = genome_seq[pos:pos + CODON_LENGTH]
        if codon in codons:
            interval_start = pos - CODON_INTERVAL_OFFSET
            interval_stop = pos + CODON_INTERVAL_OFFSET
            if interval_start < 0 or interval_stop > len(genome_seq) - CODON_LENGTH + 1:
                continue
            key = f"{chrom}:{interval_start}-{interval_stop}:+"
            if (chrom, "+") in interlap_dict:
                interlap_dict[(chrom, "+")].add((interval_start, interval_stop, key))
            else:
                inter = InterLap()
                inter.add((interval_start, interval_stop, key))
                interlap_dict[(chrom, "+")] = inter

            codon_dict[key] = [codon, 0]

        elif codon in reverse_codons:
            interval_start = pos
            interval_stop = pos + CODON_INTERVAL_SIZE - 1
            if interval_start < 0 or interval_stop > len(genome_seq) - CODON_LENGTH + 1:
                continue
            key = f"{chrom}:{interval_start}-{interval_stop}:-"
            if (chrom, "-") in interlap_dict:
                interlap_dict[(chrom, "-")].add((interval_start, interval_stop, key))
            else:
                inter = InterLap()
                inter.add((interval_start, interval_stop, key))
                interlap_dict[(chrom, "-")] = inter

            codon_dict[key] = [str(Seq(codon).reverse_complement()), 0]

    return interlap_dict, codon_dict

def get_frame(position: int) -> int:
    """
    get the reading frame for the given position
    """
    return position % CODON_LENGTH


def get_genome_information(
    start: int,
    stop: int,
    strand: str,
    genome_seq: str
) -> tuple[str, str, str, str, str]:
    """
    TODO fix for combined method
    retrieve information from genome including nucleotide sequence, start_codon, stop_codon, amino acid sequence, 15nt window
    """
    if strand == "+":
        nt_seq = genome_seq[start:stop + 1]
        aa_seq = str(Seq(nt_seq).translate(table=11, to_stop=False))

        nt_window = genome_seq[start - NT_WINDOW_SIZE:start]
    else:
        nt_seq = str(Seq(genome_seq[start:stop + 1]).reverse_complement())
        aa_seq = str(Seq(nt_seq).translate(table=11, to_stop=False))

        nt_window = str(Seq(genome_seq[stop + 1:stop + 1 + NT_WINDOW_SIZE]).reverse_complement())

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

    for gene_name, (gene_chrom, gene_start, gene_stop, gene_strand, _) in gene_dict.items():
        if gene_chrom != chrom or gene_strand != strand:
            continue

        # Annotated - exact match (highest priority)
        if gene_start == start_position and gene_stop == stop_position:
            return "Annotated", gene_name

        # Near-Annotated (second highest priority)
        if strand == "+":
            if abs(start_position - gene_start) < NEAR_ANNOTATED_THRESHOLD and stop_position == gene_stop:
                label = "Near_Annotated"
                gene_name_assigned = gene_name
                continue
        else:  # strand == "-"
            if abs(stop_position - gene_stop) < NEAR_ANNOTATED_THRESHOLD and start_position == gene_start:
                label = "Near_Annotated"
                gene_name_assigned = gene_name
                continue

        if label != "Near_Annotated":
            # Internal_Inframe (third priority)
            if strand == "+":
                if start_position > gene_start and stop_position <= gene_stop and start_position % CODON_LENGTH == gene_start % CODON_LENGTH:
                    label = "Internal_Inframe"
                    gene_name_assigned = gene_name
                    continue
            else:  # strand == "-"
                if stop_position < gene_stop and start_position >= gene_start and stop_position % CODON_LENGTH == gene_stop % CODON_LENGTH:
                    label = "Internal_Inframe"
                    gene_name_assigned = gene_name
                    continue

            # N-terminal extension (fourth priority)
            if strand == "+":
                if stop_position == gene_stop and start_position < gene_start:
                    label = "N-terminal_extension"
                    gene_name_assigned = gene_name
            else:  # strand == "-"
                if start_position == gene_start and stop_position > gene_stop:
                    label = "N-terminal_extension"
                    gene_name_assigned = gene_name

        # Internal-OutofFrame (lowest priority among overlapping categories)
        if start_position >= gene_start and stop_position <= gene_stop and start_position % CODON_LENGTH != gene_start % CODON_LENGTH:
            # Only set if no higher priority label already assigned
            if label not in ["Near_Annotated", "Internal_Inframe", "N-terminal_extension"]:
                label = "Internal_OutofFrame"
                gene_name_assigned = gene_name

    return label, gene_name_assigned if gene_name_assigned else f"{chrom}:{start_position}-{stop_position}:{strand}"


def calculate_utr_distance(
    start_position: int,
    stop_position: int,
    gene_name: str,
    gene_dict: dict[str, tuple[str, int, int, str, int]]
) -> tuple[int | float, int | float]:
    """
    calculate the relative density fiveprime and threeprime distance
    """
    start_position += 1
    stop_position += 1

    if gene_name in gene_dict:
        fiveprime_dist = start_position - gene_dict[gene_name][1]
        threeprime_dist = gene_dict[gene_name][2] - start_position
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
    calculate the relative density
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
    headers: tuple[str, str, str]
) -> pd.DataFrame:
    """
    Generate the final dataframe to be written to file.
    This contains RPKM, TE, nucleotide and aminoacid sequences and more.
    """
    tis_header, tts_header, ribo_header = headers

    te_header = expr.get_te_header(wildcards)

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count",
              f"{tis_header}_peak_height", f"{tts_header}_peak_height", f"{ribo_header}_peak_height",
              "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq",
              f"{tis_header}_relative_density", f"{tts_header}_relative_density", f"{ribo_header}_relative_density",
              "5'-distance", "3'-distance"] + [f"{card}_rpkm" for card in wildcards] + \
              [f"{cond}_TE" for cond in te_header]

    name_list = [f"s{x}" for x in range(len(header))]
    n_tuple = collections.namedtuple('Pandas', name_list)

    result_rows = []
    for (chrom, strand) in detected_orfs_dict.keys():
        for (start, stop) in detected_orfs_dict[(chrom, strand)].keys():
            rpm_tis, rpm_tts, rpm_ribo = detected_orfs_dict[(chrom, strand)][(start, stop)]

            if np.isnan(rpm_tis) and np.isnan(rpm_tts):
                continue

            if gene_tis_dict != {}:
                gene_type, gene_name = get_gene_information(chrom, start, stop, strand, gene_tis_dict)
            else:
                gene_type, gene_name = get_gene_information(chrom, start, stop, strand, gene_tts_dict)

            nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(start, stop, strand, genome[chrom])

            if aa_seq.count("*") > 1:
                continue

            rpkm_list = []
            te_list = []
            if read_count_dict != {}:
                for idx, val in enumerate(read_count_dict[(chrom, start, stop, strand)]):
                    rpkm_list.append(expr.calculate_rpkm(accepted_read_list[idx][chrom], val, len(nt_seq)))

                te_list = expr.calculate_te(rpkm_list, wildcards)

            identifier = f"{chrom}:{start + 1}-{stop + 1}:{strand}"
            codon_count = int(len(nt_seq) / CODON_LENGTH)

            if gene_tis_dict != {}:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_name, gene_tis_dict)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_name, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_name, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_name, gene_type, gene_ribo_dict)
            else:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_name, gene_tts_dict)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_name, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_name, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_name, gene_type, gene_ribo_dict)

            result = [gene_type, identifier, chrom, start + 1, stop + 1, strand, gene_name, codon_count,
                      rpm_tis, rpm_tts, rpm_ribo, start_codon, stop_codon,
                      nt_window, nt_seq, aa_seq, relative_density_tis, relative_density_tts, relative_density_ribo,
                      fiveprime_dist, threeprime_dist] + rpkm_list + te_list

            result_rows.append(n_tuple(*result))

    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])
    df_results = df_results.sort_values(by=["Genome", "Start", "Stop", "Strand"])

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


VALID_MAPPING_TYPES = {"fiveprime", "threeprime", "global", "centered"}
def base_mapping(mapping: str) -> str:
    """
    Retrieve basename of mapping (only useful for HRIBO output)
    """
    for mapping_type in VALID_MAPPING_TYPES:
        if mapping_type in mapping:
            return mapping_type

    raise ValueError(
        f"Unknown mapping type: {mapping}. Valid types: {VALID_MAPPING_TYPES}"
    )