#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import collections
import pandas as pd

import lib.misc as misc

def screen_wig_for_tss(wig_file_data, codon_interlap, codon_dict, read_count_threshold):
    """
    screen over wig file and update the according codon entries
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))
        # if read_count > x here could be a readcount restriction
        if read_count <= read_count_threshold:# change here if interval changes
            continue
        matching_codons = list(codon_interlap.find((position, position)))
        for match in matching_codons:
            codon_dict[match[2]][1] += read_count
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

def search_longest_reverse(cur_stop, genome_seq, search_codons, match_codons):
    """
    search for the match codon that is following the last inframe search codon.
    Return: Position of matching codon or -1 if not found
    """
    loop_counter = 0

    cur_position = cur_stop - 2
    nt=genome_seq[cur_position:cur_position+3]
    while nt not in search_codons or loop_counter == 0:
        loop_counter+=1
        cur_position -= 3
        if cur_position < 0:
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    while nt not in match_codons:
        cur_position += 3
        if cur_position == cur_stop - 2:
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    return cur_position

def search_longest_forward(cur_stop, genome_seq, search_codons, match_codons):
    """
    search for the match codon in reverse that is following the last inframe search codon.
    Return: Position of matching codon or -1 if not found
    """
    loop_counter = 0

    cur_position = cur_stop
    nt=genome_seq[cur_position:cur_position+3]
    while nt not in search_codons or loop_counter == 0:
        cur_position += 3
        loop_counter+=1
        if cur_position > len(genome_seq):
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    while nt not in match_codons:
        cur_position -= 3
        if cur_position == cur_stop:
            return -1

        nt = genome_seq[cur_position:cur_position+3]

    return cur_position

def detect_potential_ORFs(codon_dict, genome_seq, search_codons, match_codons, p_offset, method, detected_ORFs_dict, longest_potential_ORF=True):
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_search_codons = [str(Seq(codon).reverse_complement()) for codon in search_codons]
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]
    rows_all = []

    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, mid, strand = key.split(":")
        interval_start, interval_stop = mid.split("-")
        if method == "TIS":
            if strand == "+":
                cur_start = int(interval_start) - p_offset + 2
                cur_position = cur_start

                cur_position = search_codon_forward(cur_position, genome_seq, match_codons)
                if cur_position == -1:
                    continue

                cur_stop = cur_position + 2

            elif strand == "-":
                cur_start = int(interval_start) + p_offset + 2
                cur_position = cur_start - 2

                cur_position = search_codon_reverse(cur_position, genome_seq, reverse_match_codons)
                if cur_position == -1:
                    continue

                cur_stop = cur_position

        elif method == "TTS":
            if strand == "+":
                cur_stop = int(interval_start) - p_offset + 4
                cur_position = cur_stop - 2

                if not longest_potential_ORF:
                    cur_position = search_codon_reverse(cur_position, genome_seq, match_codons)
                    if cur_position == -1:
                        continue
                else:
                    cur_position = search_longest_reverse(cur_position, genome_seq, search_codons, match_codons)
                    if cur_position == -1:
                        continue

                cur_start = cur_position

            elif strand == "-":
                cur_stop = int(interval_start) + p_offset
                cur_position = cur_stop

                if not longest_potential_ORF:
                    cur_position = search_codon_forward(cur_position, genome_seq, reverse_match_codons)
                    if cur_position == -1:
                        continue
                else:
                    cur_position = search_longest_forward(cur_position, genome_seq, reverse_search_codons, reverse_match_codons)
                    if cur_position == -1:
                        continue

                cur_start = cur_position + 2

        if strand == "+":
            out_start, out_stop = cur_start, cur_stop
        else:
            out_start, out_stop = cur_stop, cur_start

        if (chrom, strand) in detected_ORFs_dict:
            if (out_start, out_stop) in detected_ORFs_dict[(chrom, strand)]:
                if method == "TIS":
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (val[1], detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][1])
                else:
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)][0], val[1])
            else:
                if method == "TIS":
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (val[1], -1)
                else:
                    detected_ORFs_dict[(chrom, strand)][(out_start, out_stop)] = (-1, val[1])
        else:
            if method == "TIS":
                detected_ORFs_dict[(chrom, strand)] = {(out_start, out_stop) : (val[1], -1)}
            else:
                detected_ORFs_dict[(chrom, strand)] = {(out_start, out_stop) : (-1, val[1])}

    return detected_ORFs_dict


def combined_data_detection(tis_predictions, tts_predictions, max_ORF_length):
    """
    Use the detected TIS start position and TTS stop positions to find potentially quality ORFs.
    """

    combined_ORFs_dict = {}
    for chrom, strand in tis_predictions.keys():
        try:
            for start, _, start_rpm, _ in sorted(tis_predictions[(chrom,strand)], key=lambda x: x[0]):
                for _, stop, _, stop_rpm in sorted(tts_predictions[(chrom,strand)], key=lambda x: x[1]):
                    if misc.get_frame(start) == misc.get_frame(stop):
                        if start < stop and stop - start + 1 <= max_ORF_length:
                            if (chrom, strand) in combined_ORFs_dict:
                                combined_ORFs_dict[(chrom, strand)].append((start, stop, start_rpm, stop_rpm))
                            else:
                                combined_ORFs_dict[(chrom, strand)] = []
                        else:
                            continue
        except KeyError:
            continue

    return combined_ORFs_dict
