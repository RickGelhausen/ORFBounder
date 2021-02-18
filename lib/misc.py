#!/usr/bin/env python
import re
import pandas as pd

from Bio.Seq import Seq
from interlap import InterLap


def calculate_density(wig_file_data, annotation_interlap, gene_dict):
    """
    calculate density for each annotated gene
    """

    for line in wig_file_data:
        position, read_count = line.rstrip().split(" ")
        position = int(position)-1
        read_count = abs(float(read_count))

        matching_genes = list(annotation_interlap.find((position, position)))
        for gene in matching_genes:
            gene_dict[gene[2]][4] += read_count

    return gene_dict

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
        read_list = [getattr(row, "_%s" %x) for x in range(9,len(row))]

        attribute_list = [x.strip(" ") for x in re.split('[;=]', attributes) if x != ""]

        if len(attribute_list) % 2 == 0:
            for i in range(len(attribute_list)):
                if i % 2 == 0:
                    attribute_list[i] = attribute_list[i].lower()
        else:
            print(attribute_list)
            sys.exit("error, invalid gff, wrongly formatted attribute fields.")

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

def create_codon_interlaps(chrom, genome_seq, codons, p_offset):
    """
    create interlaps around each codon, incorporating the offset
    """
    reverse_codons = [str(Seq(codon).reverse_complement()) for codon in codons]

    fwd_codon_interlap = InterLap()
    rev_codon_interlap = InterLap()
    codon_dict = {}
    for pos in range(len(genome_seq)-2):
        codon = genome_seq[pos:pos+3]
        if codon in codons:
            interval_start = pos + p_offset - 2
            interval_stop = pos + p_offset + 2
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "+")
            fwd_codon_interlap.add((interval_start, interval_stop, key))
            codon_dict[key] = [codon, 0]
        elif codon in reverse_codons:
            interval_start = pos - p_offset
            interval_stop = pos - p_offset + 4
            if interval_start < 0 or interval_stop > len(genome_seq)-2:
                continue
            key = "%s:%s-%s:%s" % (chrom, interval_start, interval_stop, "-")
            rev_codon_interlap.add((interval_start, interval_stop, key))
            codon_dict[key] = [str(Seq(codon).reverse_complement()), 0]

    return fwd_codon_interlap, rev_codon_interlap, codon_dict

def annotation_interlap(annotation_file, method):
    """
    create an interlap object for the annotation
    """
    annotation_dict = generate_annotation_dict(annotation_file)

    annotation_fwd_interlap = InterLap()
    annotation_rev_interlap = InterLap()
    gene_dict = {}
    a_codon_pos = {}

    for key, val in annotation_dict.items():
        genome, mid, strand = key.split(":")
        start, stop = mid.split("-")
        start, stop = int(start), int(stop)

        if val[1] != "":
            locus_tag = val[1]
        else:
            locus_tag = val[4]

        if strand == "+":
            annotation_fwd_interlap.add((start, stop, locus_tag))
        else:
            annotation_rev_interlap.add((start, stop, locus_tag))

        if locus_tag not in gene_dict:
            gene_dict[locus_tag] = [genome, start, stop, strand, 0]

            if method == "TIS":
                if strand == "+":
                    a_codon_pos[(genome, int(start), strand)] = (start, stop, locus_tag)
                else:
                    a_codon_pos[(genome, int(stop)-2, strand)] = (start, stop, locus_tag)
            else:
                if strand == "+":
                    a_codon_pos[(genome, int(stop)-2, strand)] = (start, stop, locus_tag)
                else:
                    a_codon_pos[(genome, int(start), strand)] = (start, stop, locus_tag)

    return annotation_fwd_interlap, annotation_rev_interlap, gene_dict, a_codon_pos


def get_frame(position):
    """
    get the reading from for the given position
    """
    return position % 3


def get_genome_information(start, stop, strand, genome_seq, method):
    """
    retrieve infomations from genome including nucleotide sequence, start_codon, stop_codon, amino acid sequence, 15nt window
    """

    if strand == "+":
        nt_seq = genome_seq[start-1:stop]
        aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))

        if method == "TIS":
            nt_window = genome_seq[start-15:start]
        else:
            nt_window = genome_seq[stop+1:stop+16]


    else:
        nt_seq = str(Seq(genome_seq[start-1:stop]).reverse_complement())
        aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))

        if method == "TIS":
            nt_window = str(Seq(genome_seq[stop+1:stop+16]).reverse_complement())
        else:
            nt_window = str(Seq(genome_seq[start-15:start]).reverse_complement())

    start_codon, stop_codon = nt_seq[:3], nt_seq[-3:]

    return nt_seq, aa_seq, nt_window, start_codon, stop_codon


def get_gene_information(chrom, start_position, stop_position, strand, gene_dict):
    # val : genome, start, stop, strand, 0
    type = "Unannotated"
    for key, val in gene_dict.items():

        if val[0] != chrom:
            continue

        gene_start, gene_stop = val[1], val[2]

        if strand == "+":
            if abs(start_position-gene_start)<10 and val[3]==strand:
                type="Near_Annotated"
                break
            elif stop_position==gene_stop and val[3]==strand:
                if start_position > gene_start:
                    type="Internal_Inframe"
                    break
                else:
                    type="N-terminal_extension"
                    break
            else:
                if start_position >= gene_start and start_position <= gene_stop and val[3]==strand:
                    type="Internal_OutofFrame"
                    break
        else:
            if abs(start_position-gene_start)<10 and val[3]==strand:
                type="Near_Annotated"
                break
            elif stop_position==gene_stop and val[3]==strand:
                if start_position < gene_start:
                    type="Internal_Inframe"
                    break
                else:
                    type="N-terminal_extension"
                    break
            else:
                if start_position <= gene_start and start_position >= gene_stop and val[3]==strand:
                    type="Internal_OutofFrame"
                    break

    if type == "Unannotated":
        key = "%s:%s-%s:%s" % (chrom, start_position, stop_position, strand)

    return type, key


def compute_additional_information(start, stop, rpm_start, rpm_stop, gene_name, gene_type, gene_dict):
    """
    calculate the relative density fiveprime and threeprime distance
    """

    if gene_name in gene_dict:
        gene_rpm = gene_dict[gene_name][4]
        if gene_rpm != 0:
            relative_density_start = rpm_start / gene_rpm
            relative_density_stop = rpm_stop / gene_rpm
        else:
            relative_density_start, relative_density_stop = "NaN", "NaN"

        if method == "TIS":
            fiveprime_dist = start - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - start
        else:
            fiveprime_dist = stop - gene_dict[gene_name][1]
            threeprime_dist = gene_dict[gene_name][2] - stop
    else:
        fiveprime_dist, threeprime_dist, relative_density_start, relative_density_stop = "NaN", "NaN", "NaN", "NaN"

    if gene_type=="N-terminal_extension":
        relative_density_start, relative_density_stop = "NaN", "NaN"

    return fiveprime_dist, threeprime_dist, relative_density_start, relative_density_stop
    
