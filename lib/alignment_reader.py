import re
import os, sys
import pysam
import interlap

import lib.messaging as msg

class PositionReader:
    """
    Read sam/bam file and create a position dict for requested read lengths
    read_lengths : list (e.g. [30,32,40])
    mapping_mode : threeprime, fiveprime, centered, global
    """
    def __init__(self, alignment_file_path, read_lengths, mapping_mode, offset_dict):
        self.alignment_file_path = alignment_file_path
        self.read_lengths = read_lengths
        self.mapping_mode = mapping_mode

        self.wildcard = re.split('_|\.', os.path.basename(alignment_file_path))[0]
        if self.wildcard in offset_dict:
            self.offset_dict = offset_dict[self.wildcard]
        elif "default" in offset_dict:
            self.offset_dict = offset_dict["default"]
        else:
            msg.error("Error: Offsets and default value missing for wildcard: %s" % self.wildcard)

        self.reads_position_dict = {}
        self.no_accepted_reads_dict = {}

        msg.message("Reading read positions from alignment file: %s" % os.path.basename(alignment_file_path))
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

                start, stop = read.reference_start, read.reference_end-1
                read_length = stop - start + 1

                strand = "-" if read.is_reverse else "+"

                if self.read_lengths != -1:
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
                    msg.error("Error: Offsets and default value missing for (wildcard, readlength): %s, %s" % (self.wildcard, read_length))

                clip_length = 11
                positions = []
                if strand == "-":
                    if self.mapping_mode == "threeprime":
                        positions = [start + offset]
                    elif self.mapping_mode == "fiveprime":
                        positions = [stop + offset]
                    elif self.mapping_mode == "centered":
                        # mid = round((read_length)/2)
                        # positions = [start+mid-1 + offset, start+mid + offset, start+mid+1 + offset]
                        center_start = start + clip_length
                        center_stop = stop - clip_length
                        positions = [i + offset for i in range(center_start, center_stop+1, 1)]
                        center_length = center_stop - center_start + 1
                    else:
                        positions = [i + offset for i in range(start, stop+1, 1)]
                else:
                    if self.mapping_mode == "threeprime":
                        positions = [stop - offset]
                    elif self.mapping_mode == "fiveprime":
                        positions = [start - offset]
                    elif self.mapping_mode == "centered":
                        # mid = round((read_length)/2)
                        # positions = [start+mid-1 - offset, start+mid - offset, start+mid+1 - offset]
                        center_start = start + clip_length
                        center_stop = stop - clip_length
                        positions = [i + offset for i in range(center_start, center_stop+1, 1)]
                        center_length = center_stop - center_start + 1
                    else:
                        positions = [i - offset for i in range(start, stop+1, 1)]

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
            sys.exit("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")

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
                            if  min_read_count > 0:
                                self.reads_position_dict[(chrom, strand)][position] *= (min_read_count / self.no_accepted_reads_dict[(chrom, strand)])
                            else:
                                msg.error(f"Error: You chose min as normalization factor, but the provided minimal read count is not valid: {min_read_count}!")
                    elif normalization == "mil":
                        self.reads_position_dict[(chrom, strand)][position] *= (1000000 / self.no_accepted_reads_dict[(chrom, strand)])
                    else:
                        msg.error(f"Error: Given normalization method is not supported: {normalization}. Supported methods: (raw, min, mil)")
        else:
            msg.message("Skipping normalization!")

    def to_wig(self, file_path):
        """
        Create two wiggle format files (+/-) for the used positions.
        """

        for (chrom, strand) in self.reads_position_dict.keys():
            s = "reverse" if strand == "-" else "forward"
            cur_file = f"{os.path.splitext(file_path)[0]}_{chrom}_{s}.wig"
            with open(cur_file, "w") as f:
                f.write(f"track type=wiggle_0 name={cur_file}\nvariableStep chrom={chrom} span=1\n")

                for position in sorted(self.reads_position_dict[(chrom,strand)].keys()):
                    value = self.reads_position_dict[(chrom,strand)][position]
                    f.write(f"{int(position)+1} {float(value)}\n")


    def output(self):
        return self.reads_position_dict, self.no_accepted_reads_dict

class IntervalReader():
    """
    Read sam/bam file and create an interlap dictionary for all reads or a selected amount of read lengths
    read_lengths : list (e.g. [30,32,40])
    rpkm_all_reads : if True, all mapped reads are used for calculating the RPKM
                     else, only the mapped reads of the given read-lengths are used
    """
    def __init__(self, alignment_file_path, read_lengths, rpkm_all_reads):
        self.alignment_file_path = alignment_file_path
        self.read_lengths = read_lengths
        self.rpkm_all_reads = rpkm_all_reads

        if self.read_lengths == -1:
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

                start, stop = read.reference_start, read.reference_end-1
                read_length = stop - start + 1

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
        return self.reads_interlap_dict, self.no_accepted_reads_dict
