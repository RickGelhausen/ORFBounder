#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def generate_genome_dict(genome_file):
    """
    read a genome fasta file into a dictionary
    """
    genome_file = SeqIO.parse(genome_file, "fasta")
    genome_dict = dict()
    for entry in genome_file:
        genome_dict[str(entry.id)] = (str(entry.seq), str(entry.seq.complement()))

    return genome_dict

def handle_input(args):
    """
    Check if input is valid.
    """

    if args.fwd_wig_file_TIS != "" and args.rev_wig_file_TIS != "" and args.fwd_wig_file_TTS != "" and args.rev_wig_file_TTS != "":
        return "allSites"

    if args.fwd_wig_file_TIS != "" and args.rev_wig_file_TIS != "":
        return "TIS"

    if args.fwd_wig_file_TTS != "" and args.rev_wig_file_TTS != "":
        return "TTS"

    sys.exit("Please ensure to either provide 2 TIS files, 2 TTS files OR both")

def load_wig(wig_path):
    with open(wig_path, 'r') as wig_file:
        chromosome = ""
        wig_data_dict = {}
        for line in wig_file.readlines():
            line = line.rstrip()

            if line[0].isdigit() and line[0] != "0":
                if chromosome not in wig_data_dict.keys():
                    sys.exit("Incomplete header in wig file! Missing chrom= field!")

                wig_data_dict[chromosome].append(line)

            elif "chrom=" in line:
                tmp = re.split('[ =]', line)
                chromosome = tmp[tmp.index("chrom")+1]
                if chromosome not in wig_data_dict:
                    wig_data_dict[chromosome] = []

    return wig_data_dict

def write_gff_file(dataframe_out, output_path, output_basename, method):
    """
    write a dataframe to a gff file
    """
    with open(os.path.join(output_path, method, output_basename), "w") as f:
        f.write("##gff-version 3\n")
    with open(os.path.join(output_path, method, output_basename), "a") as f:
        dataframe_out.to_csv(f, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE)

def write_codon_interval_gff(output_path, output_basename, codon_dict, p_offset, method):
    """
    Create a gff3 file with all codon intervals.
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    rows = []
    for key, val in codon_dict.items():
        if val[1] <= 0:
            continue
        chrom, mid, strand = key.split(":")
        start, stop = mid.split("-")

        if method == "TIS":
            if strand == "+":
                cur_position = int(start) - p_offset + 2
            elif strand == "-":
                cur_position = int(start) + p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Start_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)
        else:# change here if interval changes
            if strand == "+":
                cur_position = int(start) - p_offset + 2
            elif strand == "-":
                cur_position = int(start) + p_offset + 2

            attribute = "ID=%s;Peak_height=%s;Name=%s;Stop_codon=%s;Original_position=%s" % (key, val[1], val[0], val[0], cur_position)

        rows.append(nTuple_gff(chrom, "ORFBounder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    write_gff_file(df, output_path, output_basename, method)

def write_results_to_output_files(df_results, output_path, output_basename, split_gff, method):
    """
    write output to file,
    if split_gff == True then write one additional output file for each gene_type
    """

    nTuple_gff = collections.namedtuple('Pandas', ["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    gff_all = []
    gff_annotated = []
    gff_unannotated = []
    gff_near_annotated = []
    gff_internal_inframe = []
    gff_n_terminal = []
    gff_internal_out = []

    for row in df_results.itertuples(index=False, name='Pandas'):
        gene_type = getattr(row, "Type")
        identifier = getattr(row, "Identifier")
        start_codon = getattr(row, "start_codon")
        stop_codon = getattr(row, "stop_codon")
        peak_height = float(getattr(row, "peak_height"))
        aa_length = len(getattr(row, "aa_seq"))
        locus_tag = getattr(row, "locus_tag")

        chrom, mid, strand = identifier.split(":")
        start, stop = mid.split("-")

        attribute = "ID=%s;Name=%s;Peak_height=%s;Start_codon=%s;Stop_codon=%s;AA_length=%s;Type=%s" % (identifier, locus_tag, peak_height, start_codon, stop_codon, aa_length, gene_type)
        cur_tuple = nTuple_gff(chrom, "ORFBounder", "CDS", int(start), int(stop), ".", strand, ".", attribute)

        gff_all.append(cur_tuple)
        if split_gff:
            if gene_type == "Annotated":
                gff_annotated.append(cur_tuple)
            elif gene_type == "Unannotated":
                gff_unannotated.append(cur_tuple)
            elif gene_type == "Near_Annotated":
                gff_near_annotated.append(cur_tuple)
            elif gene_type == "Internal_Inframe":
                gff_internal_inframe.append(cur_tuple)
            elif gene_type == "N-terminal_extension":
                gff_n_terminal.append(cur_tuple)
            elif gene_type == "Internal_OutofFrame":
                gff_internal_out.append(cur_tuple)

    print("Generating gff files...")
    df_all = pd.DataFrame.from_records(gff_all, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
    write_gff_file(df_all, df_ output_path, "results_gff/%s.gff" % output_basename, method)

    if split_gff:
        df_annotated = pd.DataFrame.from_records(gff_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        df_unannotated = pd.DataFrame.from_records(gff_unannotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        df_near_annotated = pd.DataFrame.from_records(gff_near_annotated, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        df_internal_inframe = pd.DataFrame.from_records(gff_internal_inframe, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        df_n_terminal = pd.DataFrame.from_records(gff_n_terminal, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])
        df_internal_out = pd.DataFrame.from_records(gff_internal_out, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

        write_gff_file(df_annotated, output_path, "results_gff/%s_annotated.gff" % output_basename, method)
        write_gff_file(df_unannotated, output_path, "results_gff/%s_unannotated.gff" % output_basename, method)
        write_gff_file(df_near_annotated, output_path, "results_gff/%s_near_annotated.gff" % output_basename, method)
        write_gff_file(df_internal_inframe, output_path, "results_gff/%s_internal_inframe.gff" % output_basename, method)
        write_gff_file(df_n_terminal, output_path, "results_gff/%s_n_terminal.gff" % output_basename, method)
        write_gff_file(df_internal_out, output_path, "results_gff/%s_internal_out.gff" % output_basename, method)
    print("Done.")

    print("Generating output_table...")
    out_csv = os.path.join(output_path, method, "result_tables", "%s.csv" % output_basename)
    if not os.path.isfile(out_csv):
        df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    else:
        df_results.to_csv(out_csv, sep="\t", index=False, quoting=csv.QUOTE_NONE, header=False, mode="a")
    print("Done.")
