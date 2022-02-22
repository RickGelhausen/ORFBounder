#!/usr/bin/env python
import re
import collections
import numpy as np
import pandas as pd

from collections import deque

from Bio.Seq import Seq
from Bio import SeqIO

from interlap import InterLap
import lib.expression as expr
import lib.messaging as msg


def generate_annotation_dict(annotation_path):
    """
    create dictionary from annotation.
    key : (gene_id, locus_tag, name, gene_name)
    """

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
        read_list = [getattr(row, "_%s" %x) for x in range(9, len(row))]

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            msg.error("Error: invalid gff, wrongly formatted attribute fields.\n%s" % attribute_list)

        if feature.lower() == "cds":
            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            name = ""
            if "name" in attribute_list:
                name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                name = attribute_list[attribute_list.index("gene_name")+1]

            gene_id = ""
            if "gene_id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("gene_id")+1]
            elif "id" in attribute_list:
                gene_id = attribute_list[attribute_list.index("id")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            cds_dict[new_key] = (gene_id, locus_tag, name, read_list, old_locus_tag)
        elif feature.lower() in ["gene","pseudogene"]:
            gene_name = ""
            if "name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("name")+1]
            elif "gene_name" in attribute_list:
                gene_name = attribute_list[attribute_list.index("gene_name")+1]

            locus_tag = ""
            if "locus_tag" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("locus_tag")+1]
            elif "gene_id" in attribute_list:
                locus_tag = attribute_list[attribute_list.index("gene_id")+1]

            old_locus_tag = ""
            if "old_locus_tag" in attribute_list:
                old_locus_tag = attribute_list[attribute_list.index("old_locus_tag")+1]

            new_key = "%s:%s-%s:%s" % (chromosome, start, stop, strand)
            gene_dict[new_key] = (gene_name, locus_tag, old_locus_tag)

    for key in cds_dict.keys():
        gene_name = ""
        gene_id, locus_tag, name, read_list, old_locus_tag = cds_dict[key]

        if key in gene_dict:
            gene_name, gene_locus_tag, gene_old_locus_tag = gene_dict[key]

            if locus_tag == "":
                locus_tag = gene_locus_tag
            if old_locus_tag == "":
                old_locus_tag = gene_old_locus_tag

        annotation_dict[key] = (gene_id, locus_tag, name, read_list, gene_name, old_locus_tag)

    return annotation_dict

def annotation_interlap(annotation_file):
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

def calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict):
    """
    calculate density for each annotated gene
    """

    for (chrom, strand) in alignment_position_dict:
        for position, read_count in alignment_position_dict[(chrom, strand)].items():
            matching_genes = list(annotation_interlap_dict[(chrom, strand)].find((position, position)))
            for gene in matching_genes:
                gene_density_dict[gene[2]][4] += read_count

    return gene_density_dict

def create_codon_interlaps(chrom, genome_seq, codons):
    """
    create interlaps around each codon
    """
    reverse_codons = [str(Seq(codon).reverse_complement()) for codon in codons]

    interlap_dict = {}
    codon_dict = {}
    for pos in range(len(genome_seq)-2):
        codon = genome_seq[pos:pos+3]
        if codon in codons:
            interval_start = pos - 2
            interval_stop = pos + 2
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "+")
            if (chrom, "+") in interlap_dict:
                interlap_dict[(chrom, "+")].add((interval_start, interval_stop, key))
            else:
                inter = InterLap()
                inter.add((interval_start, interval_stop, key))
                interlap_dict[(chrom, "+")] = inter

            codon_dict[key] = [codon, 0]
        elif codon in reverse_codons:
            interval_start = pos
            interval_stop = pos + 4
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "-")
            if (chrom, "-") in interlap_dict:
                interlap_dict[(chrom, "-")].add((interval_start, interval_stop, key))
            else:
                inter = InterLap()
                inter.add((interval_start, interval_stop, key))
                interlap_dict[(chrom, "-")] = inter

            codon_dict[key] = [str(Seq(codon).reverse_complement()), 0]

    return interlap_dict, codon_dict

def get_frame(position):
    """
    get the reading from for the given position
    """
    return position % 3


def get_genome_information(start, stop, strand, genome_seq, method):
    """
    TODO fix for combined method
    retrieve infomations from genome including nucleotide sequence, start_codon, stop_codon, amino acid sequence, 15nt window
    """

    if strand == "+":
        nt_seq = genome_seq[start:stop+1]
        aa_seq = str(Seq(nt_seq).translate(table=11, to_stop=False))

        if method == "TIS":
            nt_window = genome_seq[start-15:start]
        else:
            nt_window = genome_seq[stop+1:stop+16]
    else:
        nt_seq = str(Seq(genome_seq[start:stop+1]).reverse_complement())
        aa_seq = str(Seq(nt_seq).translate(table=11, to_stop=False))

        if method == "TIS":
            nt_window = str(Seq(genome_seq[stop+1:stop+16]).reverse_complement())
        else:
            nt_window = str(Seq(genome_seq[start-15:start]).reverse_complement())

    start_codon, stop_codon = nt_seq[:3], nt_seq[-3:]

    return nt_seq, aa_seq, nt_window, start_codon, stop_codon


def get_gene_information(chrom, start_position, stop_position, strand, gene_dict):
    """
    determine the gene_name and gene_type
    # TODO clean up this function, check NC_002163.1:46424-49027:+, NC_002163.1:46557-46580:+
    """

    start_position += 1
    stop_position += 1
    for gene_name, (gene_chrom, gene_start, gene_stop, gene_strand, _) in gene_dict.items():
        if gene_start == start_position and gene_stop == stop_position:
            return "Annotated", gene_name

    if strand == "-":
        start_position, stop_position = stop_position, start_position

    for gene_name, (gene_chrom, gene_start, gene_stop, gene_strand, _) in gene_dict.items():
        if gene_chrom != chrom:
            continue

        if strand == "+":
            if abs(start_position-gene_start)<10 and gene_strand==strand:
                return "Near_Annotated", gene_name

            elif stop_position==gene_stop and gene_strand==strand:
                if start_position > gene_start:
                    return "Internal_Inframe", gene_name
                else:
                    return "N-terminal_extension", gene_name
            else:
                if start_position >= gene_start and start_position <= gene_stop and gene_strand==strand:
                    return "Internal_OutofFrame", gene_name
        else:
            gene_start, gene_stop = gene_stop, gene_start

            if abs(start_position-gene_start)<10 and gene_strand==strand:
                return "Near_Annotated", gene_name

            elif stop_position==gene_stop and gene_strand==strand:
                if start_position < gene_start:
                    return "Internal_Inframe", gene_name
                else:
                    return "N-terminal_extension", gene_name
            else:
                if start_position <= gene_start and start_position >= gene_stop and gene_strand==strand:
                    return "Internal_OutofFrame", gene_name

    return "Unannotated", "%s:%s-%s:%s" % (chrom, start_position, stop_position, strand)

def calculate_utr_distance(start_position, stop_position, gene_name, gene_dict, method):
    """
    calculate the relative density fiveprime and threeprime distance
    """
    start_position += 1
    stop_position += 1
    if gene_name in gene_dict:
        if method == "TIS":
            fiveprime_dist = start_position - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - start_position
        else:
            fiveprime_dist = stop_position - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - stop_position
    else:
        fiveprime_dist, threeprime_dist = np.nan, np.nan

    return fiveprime_dist, threeprime_dist

def calculate_relative_density(rpm, gene_name, gene_type, gene_dict):
    """
    calculate the relative density
    """

    if gene_dict == {} or gene_type == "N-terminal_extension" or gene_name not in gene_dict:
        return np.nan

    else:
        density = 0.0

        gene_rpm = gene_dict[gene_name][4]
        if gene_rpm != 0:
            density = rpm / gene_rpm
        else:
            density = 0.0

        return density

def generate_result_dataframe(detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict, genome, read_count_dict, \
                            accepted_read_list, wildcards, method, headers):
    """
    Generate the final dataframe to be written to file.
    This contains RPKM, TE, nucleotide and aminoacid sequences and more.
    """

    tis_header, tts_header, ribo_header = headers

    te_header = expr.get_te_header(wildcards)

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "Locus_tag", "Codon_count", \
              tis_header + "_peak_height", tts_header + "_peak_height", ribo_header + "_peak_height", \
              "Start_codon", "Stop_codon", "15nt_window", "Nucleotide_Seq", "Amino_Acid_Seq", \
              tis_header + "_relative_density", tts_header + "_relative_density", ribo_header + "_relative_density", \
              "5'-distance", "3'-distance"] + [card + "_rpkm" for card in wildcards] +\
              [cond + "_TE" for cond in te_header]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

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
            nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(start, stop, strand, genome[chrom], method)

            if aa_seq.count("*") > 1:
                continue
            rpkm_list = []
            te_list = []
            if read_count_dict != {}:
                for idx, val in enumerate(read_count_dict[(chrom, start, stop, strand)]):
                    rpkm_list.append(expr.calculate_rpkm(accepted_read_list[idx][chrom], val, len(nt_seq)))

                te_list = expr.calculate_te(rpkm_list, wildcards)

            identifier = "%s:%s-%s:%s" % (chrom, start+1, stop+1, strand)
            codon_count = int(len(nt_seq)/3)

            if gene_tis_dict != {}:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_name, gene_tis_dict, method)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_name, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_name, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_name, gene_type, gene_ribo_dict)
            else:
                fiveprime_dist, threeprime_dist = calculate_utr_distance(start, stop, gene_name, gene_tts_dict, method)
                relative_density_tis = calculate_relative_density(rpm_tis, gene_name, gene_type, gene_tis_dict)
                relative_density_tts = calculate_relative_density(rpm_tts, gene_name, gene_type, gene_tts_dict)
                relative_density_ribo = calculate_relative_density(rpm_ribo, gene_name, gene_type, gene_ribo_dict)

            result = [gene_type, identifier, chrom, start+1, stop+1, strand, gene_name, codon_count, \
                      rpm_tis, rpm_tts, rpm_ribo, start_codon, stop_codon, \
                      nt_window, nt_seq, aa_seq, relative_density_tis, relative_density_tts, relative_density_ribo, \
                      fiveprime_dist, threeprime_dist] + rpkm_list + te_list

            result_rows.append(nTuple(*result))

    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])

    df_results = df_results.sort_values(by=["Genome", "Start", "Stop", "Strand"])

    return df_results

def dictionary_depth(dic):
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

def base_mapping(mapping):
    """
    retrieve basename of mapping (only useful for HRIBO output)
    """

    if "fiveprime" in mapping:
        return "fiveprime"
    elif "threeprime" in mapping:
        return "threeprime"
    elif "global" in mapping:
        return "global"
    elif "centered" in mapping:
        return "centered"
