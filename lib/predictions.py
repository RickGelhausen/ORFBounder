#!/usr/bin/env python

"""
Module for prediction of ORF boundaries based on detected codon positions and read coverage.
"""

from Bio.Seq import Seq

import numpy as np

from lib import misc
import lib.messaging as msg

from typing import Optional

CODON_LENGTH = misc.CODON_LENGTH

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
                    msg.error("Invalid method! Use either 'sum' or 'max'!")

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
    loop_counter = 0
    original_position = cur_position
    nt = genome_seq[cur_position:cur_position + 3]
    while nt not in search_codons or loop_counter == 0:
        cur_position += 3
        loop_counter += 1
        if cur_position > len(genome_seq):
            return search_codon_forward(original_position, genome_seq, match_codons)

        nt = genome_seq[cur_position:cur_position + 3]

    while nt not in match_codons:
        cur_position -= 3
        if cur_position == original_position:
            return None

        nt = genome_seq[cur_position:cur_position + 3]

    return cur_position


def search_longest_reverse(
    cur_position: int,
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str]
) -> Optional[int]:
    """
    search for the match codon that is following the last inframe search codon.
    Return: Position of matching codon or None if not found
    """
    loop_counter = 0
    original_position = cur_position
    nt = genome_seq[cur_position:cur_position + 3]
    while nt not in search_codons or loop_counter == 0:
        loop_counter += 1
        cur_position -= 3
        if cur_position < 0:
            return search_codon_reverse(original_position, genome_seq, match_codons)

        nt = genome_seq[cur_position:cur_position + 3]

    while nt not in match_codons:
        cur_position += 3
        if cur_position == original_position:
            return None

        nt = genome_seq[cur_position:cur_position + 3]

    return cur_position

def detect_potential_orfs(
    codon_dict: dict[str, list],
    genome_seq: str,
    search_codons: list[str],
    match_codons: list[str],
    method: str,
    detected_orfs_dict: dict,
    tts_start_selection: str
) -> dict:
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_search_codons = [str(Seq(codon).reverse_complement()) for codon in search_codons]
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]

    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, mid, strand = key.split(":")
        interval_start, _ = mid.split("-")
        if method == "TIS" or method == "RIBO":
            if strand == "+":
                cur_start = int(interval_start) + 2
                cur_position = cur_start

                cur_position = search_codon_forward(cur_position, genome_seq, match_codons)
                if cur_position is None:
                    continue

                cur_stop = cur_position + 2

            elif strand == "-":
                cur_start = int(interval_start) + 2
                cur_position = cur_start - 2

                cur_position = search_codon_reverse(cur_position, genome_seq, reverse_match_codons)
                if cur_position is None:
                    continue

                cur_stop = cur_position

        elif method == "TTS":
            if strand == "+":
                cur_stop = int(interval_start) + 4
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
                cur_stop = int(interval_start)
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

def convert_codon_dict(
    codon_dict_tis: dict,
    codon_dict_tts: dict,
    offset_tis: int,
    offset_tts: int
) -> tuple[dict, dict]:
    """
    Reformat the codon_dicts.
    Convert the positions in the codon_dicts to the real positions using the offsets.
    """
    start_codon_dict = {}
    for key, val in codon_dict_tis.items():
        chrom, mid, strand = key[0].split(":")
        offset_tis = key[1]
        interval_start, _ = mid.split("-")
        if val[1] < 1:
            continue

        if strand == "+":
            cur_start = int(interval_start) - offset_tis + 2
        else:
            cur_start = int(interval_start) + offset_tis + 2

        if offset_tis in start_codon_dict:
            if (chrom, strand) in start_codon_dict[offset_tis]:
                start_codon_dict[offset_tis][(chrom, strand)].append((cur_start, val[1]))
            else:
                start_codon_dict[offset_tis][(chrom, strand)] = [(cur_start, val[1])]
        else:
            start_codon_dict[offset_tis] = {(chrom, strand): [(cur_start, val[1])]}

    stop_codon_dict = {}
    for key, val in codon_dict_tts.items():
        chrom, mid, strand = key[0].split(":")
        offset_tts = key[1]
        interval_start, _ = mid.split("-")
        if val[1] < 1:
            continue

        if strand == "+":
            cur_stop = int(interval_start) - offset_tts + 4
        else:
            cur_stop = int(interval_start) + offset_tts

        if offset_tts in stop_codon_dict:
            if (chrom, strand) in stop_codon_dict[offset_tts]:
                stop_codon_dict[offset_tts][(chrom, strand)].append((cur_stop, val[1]))
            else:
                stop_codon_dict[offset_tts][(chrom, strand)] = [(cur_stop, val[1])]
        else:
            stop_codon_dict[offset_tts] = {(chrom, strand): [(cur_stop, val[1])]}

    return start_codon_dict, stop_codon_dict

def combined_data_detection(
    codon_dict_tis: dict,
    codon_dict_tts: dict,
    offset_tis: int,
    offset_tts: int,
    max_orf_length: int
) -> dict:
    """
    Use the detected codons from TIS and TTS to find combined results.
    """
    start_codon_dict, stop_codon_dict = convert_codon_dict(codon_dict_tis, codon_dict_tts, offset_tis, offset_tts)

    keys = set()
    keys.update(start_codon_dict.keys())
    keys.update(stop_codon_dict.keys())

    predictions = {}
    for (chrom, strand) in sorted(list(keys)):
        if (chrom, strand) not in start_codon_dict or (chrom, strand) not in stop_codon_dict:
            continue

        if strand == "+":
            for stop, stop_rpm in sorted(stop_codon_dict[(chrom, strand)], key=lambda x: x[0]):
                for start, start_rpm in sorted(start_codon_dict[(chrom, strand)], key=lambda x: x[0]):
                    if start >= stop:
                        break
                    if misc.get_frame(start) != misc.get_frame(stop - 2) or abs(start - stop + 1) > max_orf_length:
                        continue
                    out_start, out_stop = start, stop
                    if (chrom, strand) in predictions:
                        predictions[(chrom, strand)][(out_start, out_stop)] = (start_rpm, stop_rpm)
                    else:
                        predictions[(chrom, strand)] = {(out_start, out_stop): (start_rpm, stop_rpm)}
        else:
            for start, start_rpm in sorted(start_codon_dict[(chrom, strand)], key=lambda x: x[0]):
                for stop, stop_rpm in sorted(stop_codon_dict[(chrom, strand)], key=lambda x: x[0]):
                    if stop >= start:
                        break
                    if misc.get_frame(start - 2) != misc.get_frame(stop) or abs(start - stop + 1) > max_orf_length:
                        continue
                    out_start, out_stop = stop, start
                    if (chrom, strand) in predictions:
                        predictions[(chrom, strand)][(out_start, out_stop)] = (start_rpm, stop_rpm)
                    else:
                        predictions[(chrom, strand)] = {(out_start, out_stop): (start_rpm, stop_rpm)}

    return predictions
