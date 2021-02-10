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

        rows.append(nTuple_gff(chrom, "TTS_finder", "codon_interval", int(start)+1, int(stop)+1, ".", strand, ".", attribute))

    df = pd.DataFrame.from_records(rows, columns=["chromosome","source","type","start","stop","score","strand","phase","attribute"])

    write_gff_file(df, output_path, output_basename, method)

def write_results_to_output_files(df_results, output_path, output_basename, split_gff, method):
    """
    Write output to file, if split_gff
    """
