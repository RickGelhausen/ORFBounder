def init_read_count_dict(read_count_dict, result_dict):
    """
    Collect intervals needed for read_counting.
    adds all intervals that are not yet in the read_count_dict from the predictions_dict.
    """

    for (chrom, strand) in result_dict.keys():
        for (start, stop) in result_dict[(chrom, strand)].keys():
            if (chrom, start, stop, strand) not in read_count_dict:
                read_count_dict[(chrom, start, stop, strand)] = []

    return read_count_dict

def create_interlap_dict():
    """
    create a dictionary with interlap objects for the current bam file.
    """

def calculate_read_count():
    """
    run over all available bam files and add read_counts for each interval in the interval dict.
    """
