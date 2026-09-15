#!/usr/bin/env python

"""
Module for prediction of ORF boundaries based on detected codon positions and read coverage.
"""

from typing import Optional

from Bio.Seq import Seq
import numpy as np

from lib import misc
from lib.coordinates import parse_coordinate_id
import lib.messaging as msg



CODON_LENGTH = misc.CODON_LENGTH
CODON_INTERVAL_OFFSET = misc.CODON_INTERVAL_OFFSET

def screen_positions_for_tss(
    alignment_position_dict: dict[tuple[str, str], dict[int, int]],
    codon_interlap_dict: dict,
    codon_dict: dict[str, list],
    min_peak_height: float,
    peak_height_operator: str
) -> dict[str, list]:
    """
    screen over wig file and update the according codon entries
    """
    for (chrom, strand) in alignment_position_dict:
        for position, read_count in alignment_position_dict[(chrom, strand)].items():

            if read_count <= min_peak_height:
                continue

            if (chrom, strand) not in codon_interlap_dict:
                continue

            matching_codons = list(codon_interlap_dict[(chrom, strand)].find((position, position)))
            for match in matching_codons:
                if peak_height_operator == "sum":
                    codon_dict[match[2]][1] += read_count
                elif peak_height_operator == "max":
                    if read_count > codon_dict[match[2]][1]:
                        codon_dict[match[2]][1] = read_count
                else:
                    raise ValueError("Invalid method! Use either 'sum' or 'max'!")

    return codon_dict




def search_codon_forward(
    cur_position: int,
    genome_seq: str,
    match_codons: list[str]
) -> Optional[int]:
    """
    search for the next matching codon 3nt at a time.
    Return: Position of matching codon or None if not found
    """
    nt = genome_seq[cur_position:cur_position + 3]
    while nt not in match_codons:
        cur_position += 3
        if cur_position > len(genome_seq) - 2:
            break
        nt = genome_seq[cur_position:cur_position + 3]

    if nt not in match_codons:
        return None

    return cur_position

def search_codon_reverse(
    cur_position: int,
    genome_seq: str,
    match_codons: list[str]
) -> Optional[int]:
    """
    search for the next matching codon in reverse 3nt at a time.
    Return: Position of matching codon or None if not found
    """
    nt = genome_seq[cur_position:cur_position + 3]
    while nt not in match_codons:
        cur_position -= 3
        if cur_position < 0:
            break
        nt = genome_seq[cur_position:cur_position + 3]

    if nt not in match_codons:
        return None

    return cur_position

#    if method == "TIS" or method == "RIBO":
#        search_codons = start_codons
#        match_codons = stop_codons
#    else:
#        search_codons = stop_codons
#        match_codons = start_codons
def search_longest_forward(
    cur_position: int,
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str]
) -> Optional[int]:
    """
    search for the match codon in reverse that is following the last inframe search codon.
    Return: Position of matching codon or None if not found
    """
    furthest = None
    for position in range(cur_position + CODON_LENGTH, len(genome_seq) - CODON_LENGTH + 1, CODON_LENGTH):
        codon = genome_seq[position:position + CODON_LENGTH]
        if codon in search_codons:
            return furthest
        if codon in match_codons:
            furthest = position
    if furthest is None and genome_seq[cur_position:cur_position + CODON_LENGTH] in match_codons:
        return cur_position
    return furthest


def search_longest_reverse(
    cur_position: int,
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str]
) -> Optional[int]:
    """
    search for the match codon that is following the last inframe search codon.
    Example: Search for the next stop codon then take the furthest upstream start codon between the current and that stop codon.
    Return: Position of matching codon or None if not found
    """
    furthest = None
    for position in range(cur_position - CODON_LENGTH, -1, -CODON_LENGTH):
        codon = genome_seq[position:position + CODON_LENGTH]
        if codon in search_codons:
            break
        if codon in match_codons:
            furthest = position
    return furthest


def _circular_codon(genome_seq: str, position: int) -> str:
    """Return a forward-reference codon, wrapping across the origin."""
    reference_length = len(genome_seq)
    return "".join(genome_seq[(position + offset) % reference_length]
                   for offset in range(CODON_LENGTH))


def _search_circular_next(
    cur_position: int,
    genome_seq: str,
    match_codons: list[str],
    direction: int,
) -> Optional[int]:
    """Return the in-frame distance to the next codon within one revolution."""
    for distance in range(0, len(genome_seq) - CODON_LENGTH + 1, CODON_LENGTH):
        position = (cur_position + direction * distance) % len(genome_seq)
        if _circular_codon(genome_seq, position) in match_codons:
            return distance
    return None


def _search_circular_furthest(
    cur_position: int,
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str],
    direction: int,
) -> Optional[int]:
    """Find the furthest in-frame match before the next search codon."""
    furthest = None
    for distance in range(CODON_LENGTH, len(genome_seq) - CODON_LENGTH + 1, CODON_LENGTH):
        position = (cur_position + direction * distance) % len(genome_seq)
        codon = _circular_codon(genome_seq, position)
        if codon in search_codons:
            break
        if codon in match_codons:
            furthest = distance
    if furthest is None and _circular_codon(genome_seq, cur_position) in match_codons:
        return 0
    return furthest


def _detect_circular_orf(
    interval_start: int,
    strand: str,
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str],
    reverse_search_codons: list[str],
    reverse_match_codons: list[str],
    method: str,
    tts_start_selection: str,
) -> Optional[tuple[int, int]]:
    """Return canonical zero-based virtual ORF coordinates on a circle."""
    reference_length = len(genome_seq)
    codon_position = ((interval_start + CODON_INTERVAL_OFFSET) % reference_length
                      if strand == "+" else interval_start % reference_length)

    if method in {"TIS", "RIBO"}:
        if strand == "+":
            distance = _search_circular_next(
                codon_position, genome_seq, match_codons, 1,
            )
            if distance is None:
                return None
            return codon_position, codon_position + distance + CODON_LENGTH - 1
        distance = _search_circular_next(
            codon_position, genome_seq, reverse_match_codons, -1,
        )
        if distance is None:
            return None
        start = (codon_position - distance) % reference_length
        return start, start + distance + CODON_LENGTH - 1

    if method != "TTS":
        raise ValueError(msg.error(f"Error! Unknown prediction method: {method}"))

    if tts_start_selection not in {"furthest_inframe", "next_inframe"}:
        raise ValueError(msg.error(
            f"Error! Unknown TTS start selection method: {tts_start_selection} "
            "expected:(furthest_inframe, next_inframe)"
        ))
    if strand == "+":
        if tts_start_selection == "next_inframe":
            distance = _search_circular_next(
                codon_position, genome_seq, match_codons, -1,
            )
        else:
            distance = _search_circular_furthest(
                codon_position, genome_seq, search_codons, match_codons, -1,
            )
        if distance is None:
            return None
        start = (codon_position - distance) % reference_length
        return start, start + distance + CODON_LENGTH - 1

    if tts_start_selection == "next_inframe":
        distance = _search_circular_next(
            codon_position, genome_seq, reverse_match_codons, 1,
        )
    else:
        distance = _search_circular_furthest(
            codon_position, genome_seq, reverse_search_codons,
            reverse_match_codons, 1,
        )
    if distance is None:
        return None
    return codon_position, codon_position + distance + CODON_LENGTH - 1

def detect_potential_orfs(
    codon_dict: dict[str, list],
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str],
    method: str,
    detected_orfs_dict: dict,
    tts_start_selection: str,
    circular: bool = False,
) -> dict:
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_search_codons = [str(Seq(codon).reverse_complement()) for codon in search_codons]
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]

    if circular and len(genome_seq) < CODON_LENGTH:
        return detected_orfs_dict

    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, interval_start, _, strand = parse_coordinate_id(key)
        if circular:
            coordinates = _detect_circular_orf(
                interval_start, strand, genome_seq, search_codons, match_codons,
                reverse_search_codons, reverse_match_codons, method,
                tts_start_selection,
            )
            if coordinates is None:
                continue
            out_start, out_stop = coordinates
        elif method in ["TIS", "RIBO"]:
            if strand == "+":
                cur_start = interval_start + 2
                cur_position = cur_start

                cur_position = search_codon_forward(cur_position, genome_seq, match_codons)
                if cur_position is None:
                    continue

                cur_stop = cur_position + 2

            elif strand == "-":
                cur_start = interval_start + 2
                cur_position = cur_start - 2

                cur_position = search_codon_reverse(cur_position, genome_seq, reverse_match_codons)
                if cur_position is None:
                    continue

                cur_stop = cur_position

        elif method == "TTS":
            if strand == "+":
                cur_stop = interval_start + 4
                cur_position = cur_stop - 2

                if tts_start_selection == "next_inframe":
                    cur_position = search_codon_reverse(cur_position, genome_seq, match_codons)
                    if cur_position is None:
                        continue
                elif tts_start_selection == "furthest_inframe":

                    cur_position = search_longest_reverse(cur_position, genome_seq, search_codons, match_codons)
                    if cur_position is None:
                        continue
                else:
                    raise ValueError(msg.error(f"Error! Unknown TTS start selection method: {tts_start_selection} expected:(furthest_inframe, next_inframe)"))

                cur_start = cur_position

            elif strand == "-":
                cur_stop = interval_start
                cur_position = cur_stop

                if tts_start_selection == "next_inframe":
                    cur_position = search_codon_forward(cur_position, genome_seq, reverse_match_codons)
                    if cur_position is None:
                        continue
                elif tts_start_selection == "furthest_inframe":
                    cur_position = search_longest_forward(cur_position, genome_seq, reverse_search_codons, reverse_match_codons)
                    if cur_position is None:
                        continue
                else:
                    raise ValueError(msg.error(f"Error! Unknown TTS start selection method: {tts_start_selection} expected:(furthest_inframe, next_inframe)"))

                cur_start = cur_position + 2

        if not circular:
            if strand == "+":
                out_start, out_stop = cur_start, cur_stop
            else:
                out_start, out_stop = cur_stop, cur_start

        if (chrom, strand) in detected_orfs_dict:
            if (out_start, out_stop) in detected_orfs_dict[(chrom, strand)]:
                tmp_tis_val = detected_orfs_dict[(chrom, strand)][(out_start, out_stop)][0]
                tmp_tts_val = detected_orfs_dict[(chrom, strand)][(out_start, out_stop)][1]
                tmp_ribo_val = detected_orfs_dict[(chrom, strand)][(out_start, out_stop)][2]

                if method == "TIS":
                    tmp_tis_val = val[1]
                elif method == "TTS":
                    tmp_tts_val = val[1]
                else:
                    tmp_ribo_val = val[1]

                detected_orfs_dict[(chrom, strand)][(out_start, out_stop)] = (tmp_tis_val, tmp_tts_val, tmp_ribo_val)

            else:
                if method == "TIS":
                    detected_orfs_dict[(chrom, strand)][(out_start, out_stop)] = (val[1], np.nan, np.nan)
                elif method == "TTS":
                    detected_orfs_dict[(chrom, strand)][(out_start, out_stop)] = (np.nan, val[1], np.nan)
                else:
                    detected_orfs_dict[(chrom, strand)][(out_start, out_stop)] = (np.nan, np.nan, val[1])
        else:
            if method == "TIS":
                detected_orfs_dict[(chrom, strand)] = {(out_start, out_stop): (val[1], np.nan, np.nan)}
            elif method == "TTS":
                detected_orfs_dict[(chrom, strand)] = {(out_start, out_stop): (np.nan, val[1], np.nan)}
            else:
                detected_orfs_dict[(chrom, strand)] = {(out_start, out_stop): (np.nan, np.nan, val[1])}

    return detected_orfs_dict
