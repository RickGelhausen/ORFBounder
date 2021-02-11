#!/usr/bin/env python

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

def get_gene_information(chrom, start_position, stop_position, strand, gene_dict):
    # val : genome, start, stop, strand, 0
    type = "Unannotated"
    for key, val in gene_dict.items():

        if val[0] != chrom:
            continue

        if strand == "+":
            gene_start, gene_stop = val[1], val[2]

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
            gene_start, gene_stop = val[2], val[1]

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



def detect_potential_ORFs(codon_dict, gene_dict, genome_seq, match_codons, a_codon_pos, p_offset, method):
    """
    for each relavent codon site, find a matching orf region
    """
    reverse_match_codons = [str(Seq(codon).reverse_complement()) for codon in match_codons]

    rows_all = []

    header = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand", "locus_tag", "codon_count", "peak_height", "start_codon", "stop_codon", "15nt window", "nt_seq", "aa_seq", "relative_density", "5'-distance", "3'-distance"]
    name_list = ["s%s" % str(x) for x in range(len(header))]
    nTuple = collections.namedtuple('Pandas', name_list)

    detected_codons_list = []
    result_rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue

        chrom, mid, strand = key.split(":")
        interval_start, interval_stop = mid.split("-")

        if method == "TIS":
            if strand == "+":
                cur_start = int(interval_start) - p_offset + 2
                cur_position = cur_start

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = genome_seq[cur_start-1:cur_stop]
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in match_codons:
                        cur_position+=3multiple
                        if cur_position > len(genome_seq)-2:
                            breakCDS
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq += nt

                    if nt not in match_codons:
                        continue

                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_stop = cur_position + 2

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

            elif strand == "-":
                cur_start = int(interval_start) + p_offset + 2
                cur_position = cur_start - 2

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = str(Seq(genome_seq[cur_start-1:cur_stop]).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in reverse_match_codons:
                        cur_position-=3
                        if cur_position < 0:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq = nt + nt_seq

                    if nt not in reverse_match_codons:
                        continue

                    nt_seq = str(Seq(nt_seq).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_stop = cur_position

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

        elif method == "TTS":
            if strand == "+":
                cur_stop = int(interval_start) - p_offset + 4
                cur_position = cur_stop - 2

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"_internal_out.gff
                    nt_seq = genome_seq[cur_start-1:cur_stop]
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt
                    while nt not in match_codons:
                        cur_position-=3
                        if cur_position < 0:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq = nt + nt_seq

                    if nt not in match_codons:
                        continue

                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_start = cur_position

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)

            elif strand == "-":
                cur_stop = int(interval_start) + p_offset
                cur_position = cur_stop

                if (chrom, cur_position+1, strand) in a_codon_pos:
                    cur_start, cur_stop, gene_name = a_codon_pos[(chrom, cur_position+1, strand)]
                    gene_type = "Annotated"
                    nt_seq = str(Seq(genome_seq[cur_start-1:cur_stop]).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                else:
                    nt=genome_seq[cur_position:cur_position+3]
                    nt_seq=nt# change here if interval changes
                    while nt not in reverse_match_codons:
                        cur_position+=3
                        if cur_position > len(genome_seq)-2:
                            break
                        nt = genome_seq[cur_position:cur_position+3]
                        nt_seq += nt

                    if nt not in reverse_match_codons:
                        continue

                    nt_seq = str(Seq(nt_seq).reverse_complement())
                    aa_seq = str(Seq(nt_seq, generic_dna).translate(table=11, to_stop=False))
                    cur_start = cur_position + 2

                    gene_type, gene_name = get_gene_information(chrom, cur_start+1, cur_stop+1, strand, gene_dict)


        if aa_seq.count("*") > 1:
            continue

        if gene_type != "Annotated":
            if strand == "+":
                out_start, out_stop = cur_start+1, cur_stop+1
            else:
                out_start, out_stop = cur_stop+1, cur_start+1
        else:
            out_start, out_stop = cur_start, cur_stop

        rpm = val[1]
        if gene_name in gene_dict:
            gene_rpm = gene_dict[gene_name][4]
            if gene_rpm != 0:
                relative_density = rpm / gene_rpm
            else:
                relative_density = "NaN"

            if method == "TIS":
                fiveprime_dist = out_start - gene_dict[gene_name][1]
                threeprime_dist = gene_dict[gene_name][2] - out_start
            else:
                fiveprime_dist = out_stop - gene_dict[gene_name][1]
                threeprime_dist = gene_dict[gene_name][2] - out_stop
        else:
            fiveprime_dist, threeprime_dist, relative_density = "NaN", "NaN", "NaN"

        start_codon, stop_codon = nt_seq[:3], nt_seq[-3:]

        if gene_type=="N-terminal_extension":
            relative_density = "NaN"

        if method == "TIS":
            if strand == "+":
                nt_window = genome_seq[out_start-16:out_start-1]
            else:
                nt_window = str(Seq(genome_seq[out_stop:out_stop+15]).reverse_complement())
        else:
            if strand == "+":
                nt_window = genome_seq[out_stop:out_stop+15]
            else:
                nt_window = str(Seq(genome_seq[out_start-16:out_start-1]).reverse_complement())

        if method == "TIS":
            detected_codons_list.append(out_start)
        else:
            detected_codons_list.append(out_stop)

        unique_id="%s:%s-%s:%s" % (chrom, out_start, out_stop, strand)
        result = [gene_type, unique_id, chrom, out_start, out_stop, strand, gene_name, int(len(nt_seq)/3), rpm, start_codon, stop_codon, nt_window, nt_seq, aa_seq, relative_density, fiveprime_dist, threeprime_dist]
        result_rows.append(nTuple(*result))

    df_results = pd.DataFrame.from_records(result_rows, columns=[header[x] for x in range(len(header))])

    return df_results, detected_codons_list
