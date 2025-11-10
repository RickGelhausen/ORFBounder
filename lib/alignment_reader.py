"""
Class to read alignment files and create position dictionaries or interlap dictionaries
for read counting.
"""

from pathlib import Path

import re
import os
import sys
import pysam
import interlap


import lib.messaging as msg

class PositionReader:
    """
    Read sam/bam file and create a position dict for requested read lengths
    read_lengths : list (e.g. [30,32,40])
    mapping_mode : threeprime, fiveprime, centered, global
    """
    def __init__(self, alignment_file_path, read_length_dict, mapping_mode, offset_dict):
        self.alignment_file_path = alignment_file_path
        self.read_length_dict = read_length_dict
        self.mapping_mode = mapping_mode


        self.wildcard = re.split(r'_|\.', os.path.basename(alignment_file_path))[0]

        if not self.read_length_dict:
            self.read_lengths = None
        elif self.wildcard in self.read_length_dict:
            self.read_lengths = self.read_length_dict[self.wildcard]
        elif "default" in self.read_length_dict:
            self.read_lengths = self.read_length_dict["default"]
        else:
            msg.warning("Warning no default value given for read-lengths. Using all read lengths.")
            self.read_lengths = None

        if self.wildcard in offset_dict:
            self.offset_dict = offset_dict[self.wildcard]
        elif "default" in offset_dict:
            self.offset_dict = offset_dict["default"]
        else:
            msg.error(f"Error: Offsets and default value missing for wildcard: {self.wildcard}")

        self.reads_position_dict = {}
        self.no_accepted_reads_dict = {}

        msg.message(f"Reading read positions from alignment file: {os.path.basename(alignment_file_path)}")
        self._read_alignment_file()

    def _read_alignment_file(self):
        """
        read the alignment file using pysam
        """

        alignment_file = pysam.AlignmentFile(self.alignment_file_path)
        try:
            for read in alignment_file.fetch():
                chrom = read.reference_name

                if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                    continue

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                strand = "-" if read.is_reverse else "+"

                if self.read_lengths is not None:
                    if str(read_length) not in self.read_lengths:
                        continue

                if chrom in self.no_accepted_reads_dict:
                    self.no_accepted_reads_dict[chrom] += 1
                else:
                    self.no_accepted_reads_dict[chrom] = 1

                if read_length in self.offset_dict:
                    offset = self.offset_dict[read_length]
                elif "default" in self.offset_dict:
                    offset = self.offset_dict["default"]
                else:
                    msg.error(f"Error: Offsets and default value missing for (wildcard, readlength): {self.wildcard}, {read_length}")
                    sys.exit()

                clip_length = 11
                positions = []
                center_length = None
                if strand == "-":
                    if self.mapping_mode == "threeprime":
                        positions = [start + offset]
                    elif self.mapping_mode == "fiveprime":
                        positions = [stop + offset]
                    elif self.mapping_mode == "centered":
                        center_start = start + clip_length
                        center_stop = stop - clip_length
                        positions = [i + offset for i in range(center_start, center_stop + 1, 1)]
                        center_length = center_stop - center_start + 1
                    else:
                        positions = [i + offset for i in range(start, stop + 1, 1)]
                else:
                    if self.mapping_mode == "threeprime":
                        positions = [stop - offset]
                    elif self.mapping_mode == "fiveprime":
                        positions = [start - offset]
                    elif self.mapping_mode == "centered":
                        center_start = start + clip_length
                        center_stop = stop - clip_length
                        positions = [i + offset for i in range(center_start, center_stop + 1, 1)]
                        center_length = center_stop - center_start + 1
                    else:
                        positions = [i - offset for i in range(start, stop + 1, 1)]

                for pos in positions:
                    if pos < 0:
                        continue

                    if self.mapping_mode == "centered":
                        if (chrom, strand) in self.reads_position_dict:
                            if pos in self.reads_position_dict[(chrom, strand)]:
                                self.reads_position_dict[(chrom, strand)][pos] += (1 / center_length)
                            else:
                                self.reads_position_dict[(chrom, strand)][pos] = (1 / center_length)
                        else:
                            self.reads_position_dict[(chrom, strand)] = {pos : (1 / center_length)}

                    else:
                        if (chrom, strand) in self.reads_position_dict:
                            if pos in self.reads_position_dict[(chrom, strand)]:
                                self.reads_position_dict[(chrom, strand)][pos] += 1
                            else:
                                self.reads_position_dict[(chrom, strand)][pos] = 1
                        else:
                            self.reads_position_dict[(chrom, strand)] = {pos : 1}

        except ValueError:
            msg.error("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")
            sys.exit()

    def normalize_read_counts(self, normalization, min_read_count_dict):
        """
        Normalize all read counts for all positions
        """

        if normalization != "raw":
            for (chrom, strand) in self.reads_position_dict:
                for position in self.reads_position_dict[(chrom, strand)]:

                    if normalization == "min":

                        if chrom not in min_read_count_dict:
                            msg.warning(f"No minimum read count found for {chrom}. Skipping normalization.")

                        else:
                            min_read_count = min_read_count_dict[chrom]

                            if min_read_count > 0:
                                self.reads_position_dict[(chrom, strand)][position] *= (min_read_count / self.no_accepted_reads_dict[chrom])
                            else:
                                msg.error(f"Error: You chose min as normalization factor, but the provided minimal read count is not valid: {min_read_count}!")
                    elif normalization == "mil":
                        self.reads_position_dict[(chrom, strand)][position] *= (1000000 / self.no_accepted_reads_dict[chrom])
                    else:
                        msg.error(f"Error: Given normalization method is not supported: {normalization}. Supported methods: (raw, min, mil)")
        else:
            msg.message("Skipping normalization!")

    def to_wig(self, file_path):
        """
        Create two wiggle format files (+/-) for the used positions.
        """

        for (chrom, strand) in self.reads_position_dict:
            orientation = "reverse" if strand == "-" else "forward"
            cur_file = Path(file_path) / f"{Path(self.alignment_file_path).stem}_{chrom}_{orientation}.wig"
            Path(cur_file.parent).mkdir(parents=True, exist_ok=True)
            with open(cur_file, "w", encoding="utf-8") as f:
                f.write(f"track type=wiggle_0 name={Path(cur_file).name}\nvariableStep chrom={chrom} span=1\n")

                for position in sorted(self.reads_position_dict[(chrom,strand)].keys()):
                    value = self.reads_position_dict[(chrom,strand)][position]
                    f.write(f"{int(position)+1} {float(value)}\n")


    def output(self):
        """
        return the position dictionary and number of accepted reads dictionary"""
        return self.reads_position_dict, self.no_accepted_reads_dict

class IntervalReader():
    """
    Read sam/bam file and create an interlap dictionary for all reads or a selected amount of read lengths
    read_lengths : list (e.g. [30,32,40])
    rpkm_all_reads : if True, all mapped reads are used for calculating the RPKM
                     else, only the mapped reads of the given read-lengths are used
    """
    def __init__(self, alignment_file_path, read_length_dict, rpkm_all_reads):
        self.alignment_file_path = alignment_file_path
        self.read_length_dict = read_length_dict
        self.rpkm_all_reads = rpkm_all_reads

        self.wildcard = re.split('_|\.', os.path.basename(alignment_file_path))[0]
        if self.read_length_dict is None:
            self.read_lengths = None
        elif self.wildcard in self.read_length_dict:
            self.read_lengths = self.read_length_dict[self.wildcard]
        elif "default" in self.read_length_dict:
            self.read_lengths = self.read_length_dict["default"]
        else:
            msg.warning("Warning no default value given for read-lengths. Using all read lengths.")
            self.read_lengths = None

        if self.read_lengths is None:
            self.rpkm_all_reads = True

        self.no_accepted_reads_dict = {}
        self.reads_interlap_dict = {}

        self._read_alignment_file()

    def _read_alignment_file(self):
        """
        read the alignment file using pysam
        """

        alignment_file = pysam.AlignmentFile(self.alignment_file_path)
        tmp_dict = {}
        try:
            for read in alignment_file.fetch():
                chrom = read.reference_name

                if read.get_tag("NH") > 1 or read.mapping_quality < 0 or read.is_unmapped:
                    continue

                start = read.reference_start
                read_length = read.query_length # query read length
                stop = start + read_length - 1

                strand = "-" if read.is_reverse else "+"

                if not self.rpkm_all_reads:
                    if str(read_length) not in self.read_lengths:
                        continue

                if chrom in self.no_accepted_reads_dict:
                    self.no_accepted_reads_dict[chrom] += 1
                else:
                    self.no_accepted_reads_dict[chrom] = 1

                if (chrom, strand) in tmp_dict:
                    tmp_dict[(chrom, strand)].append((start, stop))
                else:
                    tmp_dict[(chrom, strand)] = [(start, stop)]

            for key, val in tmp_dict.items():
                inter = interlap.InterLap()
                inter.update(val)
                self.reads_interlap_dict[key] = inter

        except ValueError:
            sys.exit("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

    def output(self):
        """
        return the interlap dictionary and number of accepted reads dictionary
        """
        return self.reads_interlap_dict, self.no_accepted_reads_dict
