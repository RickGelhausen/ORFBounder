#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO

import collections
import pandas as pd

import lib.misc as misc

def screen_wig_for_tss(wig_file_data, codon_interlap, codon_dict, min_peak_height, peak_height_calculation):
    """
    screen over wig file and update the according codon entries
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))
        if read_count <= min_peak_height:
            continue

        for offset in codon_interlap.keys():
            matching_codons = list(codon_interlap[offset].find((position, position)))
            for match in matching_codons:
                if peak_height_calculation == "sum":
                    codon_dict[(match[2], offset)][1] += read_count
                elif peak_height_calculation == "max":
                    if read_count > codon_dict[match[2]][1]:
                        codon_dict[(match[2], offset)][1] = read_count
                else:
                    msg.error("Invalid method! Use either 'sum' or 'max'!")

    return codon_dict

def screen_area_for_tss(wig_file_data, codon_interlap, codon_dict):
    """
    screen over wig file and update the according codon entries
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))

        for offset in codon_interlap.keys():
            matching_codons = list(codon_interlap[offset].find((position, position)))
            for match in matching_codons:
                codon_dict[(match[2],offset)][1] += read_count

    return codon_dict

def search_codon_forward(cur_position, genome_seq, match_codons):
    """
    search for the next matching codon 3nt at a time.
    Return: Position of matching codon or -1 if not found
    """

    nt=genome_seq[cur_position:cur_position+3]
    while nt not in match_codons:
        cur_position+=3
        if cur_position > len(genome_seq)-2:
            break
        nt = genome_seq[cur_position:cur_position+3]

    if nt not in match_codons:
        return -1

    return cur_position

def search_codon_reverse(cur_position, genome_seq, match_codons):
    """
    search for the next matching codon in reverse 3nt at a time.
    Return: Position of matching codon or -1 if not found
    """
    nt=genome_seq[cur_position:cur_position+3]
    while nt not in match_codons:
        cur_position-=3
        if cur_position < 0:
            break
        nt = genome_seq[cur_position:cur_position+3]

    if nt not in match_codons:
        return -1

    return cur_position

def search_longest_reverse(cur_position, genome_seq, search_codons, match_codons):
    """
    search for the match codon that is following the last inframe search codon.
    Return: Position of matching codon or -1 if not found
    """
    loop_counter = 0
    original_position = cur_position
    nt=genome_seq[cur_position:cur_position+3]
    while nt not in search_codons or loop_counter == 0:
        loop_counter += 1
        cur_position -= 3
        if cur_position < 0:
            return search_codon_reverse(original_position, genome_seq, match_codons)

        nt = genome_seq[cur_position:cur_position+3]

    while nt not in match_codons:
        cur_position += 3
        if cur_position == original_position:
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    return cur_position

def search_longest_forward(cur_position, genome_seq, search_codons, match_codons):
    """
    search for the match codon in reverse that is following the last inframe search codon.
    Return: Position of matching codon or -1 if not found
    """
    loop_counter = 0
    original_position = cur_position
    nt=genome_seq[cur_position:cur_position+3]
    while nt not in search_codons or loop_counter == 0:
        cur_position += 3
        loop_counter+=1
        if cur_position > len(genome_seq):
            return search_codon_forward(original_position, genome_seq, match_codons)

        nt = genome_seq[cur_position:cur_position+3]

    while nt not in match_codons:
        cur_position -= 3
        if cur_position == original_position:
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    return cur_position

def detect_potential_ORFs(codon_dict, genome_seq, search_codons, match_codons, p_offset, method, detected_ORFs_dict, TTS_start_selection):
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_search_codons = [str(Seq(codon).reverse_complement()) for codon in search_codons]
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]
    rows_all = []

    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, mid, strand = key[0].split(":")
        offset = key[1]
        interval_start, interval_stop = mid.split("-")
        if method == "TIS":
            if strand == "+":
                cur_start = int(interval_start) - offset + 2
                cur_position = cur_start

                cur_position = search_codon_forward(cur_position, genome_seq, match_codons)
                if cur_position == -1:
                    continue

                cur_stop = cur_position + 2

            elif strand == "-":
                cur_start = int(interval_start) + offset + 2
                cur_position = cur_start - 2

                cur_position = search_codon_reverse(cur_position, genome_seq, reverse_match_codons)
                if cur_position == -1:
                    continue

                cur_stop = cur_position

        elif method == "TTS":
            if strand == "+":
                cur_stop = int(interval_start) - offset + 4
                cur_position = cur_stop - 2

                if TTS_start_selection == "next_inframe":
                    cur_position = search_codon_reverse(cur_position, genome_seq, match_codons)
                    if cur_position == -1:
                        continue
                elif TTS_start_selection == "furthest_inframe":
                    cur_position = search_longest_reverse(cur_position, genome_seq, search_codons, match_codons)
                    if cur_position == -1:
                        continue
                else:
                    msg.error("Error! Unknown TTS start selection method: %s expected:{furthest_inframe, next_inframe}" % TTS_start_selection)

                cur_start = cur_position

            elif strand == "-":
                cur_stop = int(interval_start) + offset
                cur_position = cur_stop

                if TTS_start_selection == "next_inframe":
                    cur_position = search_codon_forward(cur_position, genome_seq, reverse_match_codons)
                    if cur_position == -1:
                        continue
                elif TTS_start_selection == "furthest_inframe":
                    cur_position = search_longest_forward(cur_position, genome_seq, reverse_search_codons, reverse_match_codons)
                    if cur_position == -1:
                        continue
                else:
                    msg.error("Error! Unknown TTS start selection method: %s expected:{furthest_inframe, next_inframe}" % TTS_start_selection)

                cur_start = cur_position + 2

        if strand == "+":
            out_start, out_stop = cur_start, cur_stop
        else:
            out_start, out_stop = cur_stop, cur_start
        # if offset in detected_ORFs_dict:
        #     if (chrom, strand) in detected_ORFs_dict[offset]:
        #         if (out_start, out_stop) in detected_ORFs_dict[offset][(chrom, strand)]:
        #             if method == "TIS":
        #                 detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)] = (val[1], detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)][1])
        #             else:
        #                 detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)] = (detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)][0], val[1])
        #         else:
        #             if method == "TIS":
        #                 detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)] = (val[1], -1)
        #             else:
        #                 detected_ORFs_dict[offset][(chrom, strand)][(out_start, out_stop)] = (-1, val[1])
        #     else:
        #         if method == "TIS":
        #             detected_ORFs_dict[offset][(chrom, strand)] = {(out_start, out_stop) : (val[1], -1)}
        #         else:
        #             detected_ORFs_dict[offset][(chrom, strand)] = {(out_start, out_stop) : (-1, val[1])}
        # else:
        #     if method == "TIS":
        #         detected_ORFs_dict[offset] = {(chrom, strand) : {(out_start, out_stop) : (val[1], -1)}}
        #     else:
        #         detected_ORFs_dict[offset] = {(chrom, strand) : {(out_start, out_stop) : (-1, val[1])}}

        if (chrom, strand) in detected_ORFs_dict:
            if (out_start, out_stop) in detected_ORFs_dict[(chrom, strand)]:
                if method == "TIS":
                    tmp_val = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][0]
                    tmp_val.append(val[1])
                    tmp_offset = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][2]
                    tmp_offset.append(offset)
                    tmp_tts_val = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][1]
                    tmp_tts_offset = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][3]
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (tmp_val, tmp_tts_val, tmp_offset, tmp_tts_offset)
                else:
                    tmp_val = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][1]
                    tmp_val.append(val[1])
                    tmp_offset = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][3]
                    tmp_offset.append(offset)
                    tmp_tis_val = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][0]
                    tmp_tis_offset = detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][2]
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (tmp_tis_val, tmp_val, tmp_tis_offset, tmp_offset)
            else:
                if method == "TIS":
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = ([val[1]], [], [offset], [])
                else:
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = ([], [val[1]], [], [offset])
        else:
            if method == "TIS":
                detected_ORFs_dict[(chrom, strand)] = {(out_start, out_stop) : ([val[1]], [], [offset], [])}
            else:
                detected_ORFs_dict[(chrom, strand)] = {(out_start, out_stop) : ([], [val[1]], [], [offset])}

    return detected_ORFs_dict

def convert_codon_dict(codon_dict_TIS, codon_dict_TTS, offset_TIS, offset_TTS):
    """
    Reformat the codon_dicts.
    Convert the positions in the codon_dicts to the real positions using the offsets.
    """

    start_codon_dict = {}
    for key, val in codon_dict_TIS.items():
        chrom, mid, strand = key[0].split(":")
        offset_TIS = key[1]
        interval_start, interval_stop = mid.split("-")
        if val[1] < 1:
            continue

        if strand == "+":
            cur_start = int(interval_start) - offset_TIS + 2
        else:
            cur_start = int(interval_start) + offset_TIS + 2

        if offset_TIS in start_codon_dict:
            if (chrom, strand) in start_codon_dict[offset_TIS]:
                start_codon_dict[offset_TIS][(chrom, strand)].append((cur_start, val[1]))
            else:
                start_codon_dict[offset_TIS][(chrom, strand)] = [(cur_start, val[1])]
        else:
            start_codon_dict[offset_TIS] = {(chrom, strand) : [(cur_start, val[1])] }

    stop_codon_dict = {}
    for key, val in codon_dict_TTS.items():
        chrom, mid, strand = key[0].split(":")
        offset_TTS = key[1]
        interval_start, interval_stop = mid.split("-")
        if val[1] < 1:
            continue

        if strand == "+":
            cur_stop = int(interval_start) - offset_TTS + 4
        else:
            cur_stop = int(interval_start) + offset_TTS

        if offset_TTS in stop_codon_dict:
            if (chrom, strand) in stop_codon_dict[offset_TTS]:
                stop_codon_dict[offset_TTS][(chrom, strand)].append((cur_stop, val[1]))
            else:
                stop_codon_dict[offset_TTS][(chrom, strand)] = [(cur_stop, val[1])]
        else:
            stop_codon_dict[offset_TTS] = { (chrom, strand) : [(cur_stop, val[1])] }


    return start_codon_dict, stop_codon_dict

def combined_data_detection(codon_dict_TIS, codon_dict_TTS, offset_TIS, offset_TTS, max_ORF_length):
    """
    Use the detected codons from TIS and TTS to find combined results.
    """

    start_codon_dict, stop_codon_dict \
                = convert_codon_dict(codon_dict_TIS, codon_dict_TTS, offset_TIS, offset_TTS)

    keys = set()
    keys.update(start_codon_dict.keys())
    keys.update(stop_codon_dict.keys())

    predictions = {}
    for (chrom, strand) in sorted(list(keys)):
        if (chrom, strand) not in start_codon_dict or (chrom, strand) not in stop_codon_dict:
            continue

        if strand == "+":
            for stop, stop_rpm in sorted(stop_codon_dict[(chrom, strand)], key=lambda x : x[0]):
                for start, start_rpm in sorted(start_codon_dict[(chrom, strand)], key=lambda x : x[0]):
                    if start >= stop:
                        break
                    if misc.get_frame(start) != misc.get_frame(stop-2) or abs(start-stop+1) > max_ORF_length:
                        continue
                    out_start, out_stop = start, stop
                    if (chrom, strand) in predictions:
                        predictions[(chrom, strand)][(out_start, out_stop)] = (start_rpm, stop_rpm)
                    else:
                        predictions[(chrom, strand)] = {(out_start, out_stop) : (start_rpm, stop_rpm)}
        else:
            for start, start_rpm in sorted(start_codon_dict[(chrom, strand)], key=lambda x : x[0]):
                for stop, stop_rpm in sorted(stop_codon_dict[(chrom, strand)], key=lambda x : x[0]):
                    if stop >= start:
                        break
                    if misc.get_frame(start-2) != misc.get_frame(stop) or abs(start-stop+1) > max_ORF_length:
                        continue
                    out_start, out_stop = stop, start
                    if (chrom, strand) in predictions:
                        predictions[(chrom, strand)][(out_start, out_stop)] = (start_rpm, stop_rpm)
                    else:
                        predictions[(chrom, strand)] = {(out_start, out_stop) : (start_rpm, stop_rpm)}

    return predictions

# def combined_data_detection(tis_predictions, tts_predictions, max_ORF_length):
#     """
#     Use the detected TIS start position and TTS stop positions to find potentially quality ORFs.
#     """
#
#     combined_ORFs_dict = {}
#     for chrom, strand in tis_predictions.keys():
#         try:
#             for start, _, start_rpm, _ in sorted(tis_predictions[(chrom,strand)], key=lambda x: x[0]):
#                 for _, stop, _, stop_rpm in sorted(tts_predictions[(chrom,strand)], key=lambda x: x[1]):
#                     if misc.get_frame(start) == misc.get_frame(stop):
#                         if start < stop and stop - start + 1 <= max_ORF_length:
#                             if (chrom, strand) in combined_ORFs_dict:
#                                 combined_ORFs_dict[(chrom, strand)].append((start, stop, start_rpm, stop_rpm))
#                             else:
#                                 combined_ORFs_dict[(chrom, strand)] = []
#                         else:
#                             continue
#         except KeyError:
#             continue
#
#     return combined_ORFs_dict
