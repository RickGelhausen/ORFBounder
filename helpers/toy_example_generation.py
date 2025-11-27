#!/usr/bin/env python3
"""Generate a proper working toy example for ORFBounder."""

import os
import json
from pathlib import Path
import pysam


class ToyExampleGenerator:
    def __init__(self, output_dir="toy_example"):
        self.output_dir = Path(output_dir)
        self.data_dir = self.output_dir / "data"
        self.bam_dir = self.data_dir / "bam"
        self.config_dir = self.output_dir / "config"

        self.bam_dir.mkdir(parents=True, exist_ok=True)
        self.config_dir.mkdir(parents=True, exist_ok=True)

        # Use standard bacterial chromosome identifier
        self.chrom = "NC_000913.3"

        # Define 3 genes at realistic positions (not starting at 1)
        self.genes = [
            {"name": "gene001", "start": 100, "stop": 162},  # 63bp = 21 codons
            {"name": "gene002", "start": 300, "stop": 362},  # 63bp
            {"name": "gene003", "start": 500, "stop": 562},  # 63bp
        ]

    def generate_genome(self):
        """Generate genome with 3 genes at specific positions."""
        print("Generating genome...")

        # Total length 1000bp
        sequence = ["N"] * 1000

        # Gene 1 at position 100-162 (0-indexed: 99-161)
        gene1 = "ATGAAACCGGGTTTCAAGGCCATTTTGAACCCGGGTTTCAAGCCATTTTGAAACCGGGTAG"
        for i, base in enumerate(gene1):
            sequence[99 + i] = base

        # Gene 2 at position 300-362 (0-indexed: 299-361)
        gene2 = "ATGTTACCGGGTTTCAAGGCCATTTTCAACCCGGGTTTCAAGCCATTTACAAACCGGAAGTAA"
        for i, base in enumerate(gene2):
            sequence[299 + i] = base

        # Gene 3 at position 500-562 (0-indexed: 499-561)
        gene3 = "ATGGGACCGGGTTTCAAGGCCATTTTGAACCCGGGTTTCAAGCCATTTTGAAACCGGGTGA"
        for i, base in enumerate(gene3):
            sequence[499 + i] = base

        sequence = "".join(sequence)

        genome_file = self.data_dir / "genome.fa"
        with open(genome_file, 'w', encoding='utf-8') as f:
            f.write(f">{self.chrom}\n")
            # Write in 60bp lines
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

        print(f"  Created: {genome_file}")
        return sequence

    def generate_annotation(self):
        """Generate GFF3 annotation matching the genome."""
        print("Generating annotation...")

        gff_file = self.data_dir / "annotation.gff"
        with open(gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write(f"{self.chrom}\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID={self.chrom};Name={self.chrom}\n")

            for gene in self.genes:
                # Gene feature
                f.write(f"{self.chrom}\tRefSeq\tgene\t{gene['start']}\t{gene['stop']}\t.\t+\t.\t"
                       f"ID=gene-{gene['name']};Name={gene['name']};locus_tag={gene['name']}\n")
                # CDS feature
                f.write(f"{self.chrom}\tRefSeq\tCDS\t{gene['start']}\t{gene['stop']}\t.\t+\t0\t"
                       f"ID=cds-{gene['name']};Parent=gene-{gene['name']};Name={gene['name']};locus_tag={gene['name']}\n")

        print(f"  Created: {gff_file}")

    def generate_bam_files(self, genome_seq):
        """Generate BAM files with reads at start codons."""
        print("Generating BAM files...")

        # With offset 0, reads should start exactly at the start codon positions
        # Gene starts: 100, 300, 500
        # Use 30nt reads

        tis_reads = {
            "TIS-WT-1": [
                # Many reads at each start codon
                (100, genome_seq[99:129]),  # Gene 1 start
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (300, genome_seq[299:329]),  # Gene 2 start
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (500, genome_seq[499:529]),  # Gene 3 start
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
            ],
            "TIS-WT-2": [
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (100, genome_seq[99:129]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (300, genome_seq[299:329]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
                (500, genome_seq[499:529]),
            ],
        }

        # RIBO reads distributed across the genes
        ribo_reads = {
            "RIBO-WT-1": [
                (105, genome_seq[104:134]),
                (110, genome_seq[109:139]),
                (115, genome_seq[114:144]),
                (120, genome_seq[119:149]),
                (305, genome_seq[304:334]),
                (310, genome_seq[309:339]),
                (315, genome_seq[314:344]),
                (505, genome_seq[504:534]),
                (510, genome_seq[509:539]),
                (515, genome_seq[514:544]),
            ],
            "RIBO-WT-2": [
                (108, genome_seq[107:137]),
                (113, genome_seq[112:142]),
                (118, genome_seq[117:147]),
                (123, genome_seq[122:152]),
                (308, genome_seq[307:337]),
                (313, genome_seq[312:342]),
                (318, genome_seq[317:347]),
                (508, genome_seq[507:537]),
                (513, genome_seq[512:542]),
                (518, genome_seq[517:547]),
            ],
        }

        all_samples = {**tis_reads, **ribo_reads}

        for sample_name, reads in all_samples.items():
            sam_file = self.data_dir / f"{sample_name}.sam"

            with open(sam_file, 'w', encoding='utf-8') as f:
                f.write("@HD\tVN:1.0\tSO:coordinate\n")
                f.write(f"@SQ\tSN:{self.chrom}\tLN:1000\n")
                f.write("@PG\tID:toy\tPN:toy\tVN:1.0\n")

                for i, (pos, seq) in enumerate(reads, 1):
                    qual = "I" * len(seq)
                    cigar = f"{len(seq)}M"
                    f.write(f"read_{i:03d}\t0\t{self.chrom}\t{pos}\t255\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNH:i:1\n")

            # Convert to BAM
            bam_file = self.bam_dir / f"{sample_name}.bam"
            pysam.view("-bS", "-o", str(bam_file), str(sam_file), catch_stdout=False)
            sorted_bam = self.bam_dir / f"{sample_name}.sorted.bam"
            pysam.sort("-o", str(sorted_bam), str(bam_file), catch_stdout=False)
            os.remove(bam_file)
            os.rename(sorted_bam, bam_file)
            pysam.index(str(bam_file))
            os.remove(sam_file)

            print(f"  Created: {bam_file}")

    def generate_configs(self):
        """Generate all config files."""
        print("Generating config files...")

        # Offsets - using 0 for simplicity
        with open(self.config_dir / "offsets.json", 'w', encoding='utf-8') as f:
            json.dump({"default": {"default": 0}}, f, indent=2)

        # Read lengths
        with open(self.config_dir / "read_lengths.json", 'w', encoding='utf-8') as f:
            json.dump({"default": "30"}, f)

        # Mapped counts - matching the actual read counts
        with open(self.config_dir / "mapped_counts.tsv", 'w', encoding='utf-8') as f:
            f.write(f"TIS-WT-1\t{self.chrom}\t24\n")
            f.write(f"TIS-WT-2\t{self.chrom}\t20\n")
            f.write(f"RIBO-WT-1\t{self.chrom}\t10\n")
            f.write(f"RIBO-WT-2\t{self.chrom}\t10\n")

        # Main config
        with open(self.config_dir / "config.tsv", 'w', encoding='utf-8') as f:
            f.write("\t".join([
                "experiment_name", "annotation_file_path", "genome_file_path",
                "TIS_folder_path", "TTS_folder_path", "RIBO_folder_path",
                "read_length_json", "normalization_method", "mapped_counts_file_path",
                "mapping_method", "offset_file_path", "start_codons", "stop_codons",
                "alignment_folder_path", "min_peak_height", "peak_height_operator",
                "tts_start_selection", "max_ORF_length", "rpkm_read_usage", "gff_output_mode"
            ]) + "\n")

            f.write("\t".join([
                "toy_example", "toy_example/data/annotation.gff", "toy_example/data/genome.fa",
                "toy_example/data/bam", "toy_example/data/bam", "toy_example/data/bam",
                "toy_example/config/read_lengths.json", "raw", "toy_example/config/mapped_counts.tsv",
                "fiveprime", "toy_example/config/offsets.json", "ATG", "TAG,TAA,TGA",
                "toy_example/data/bam", "1", "sum", "", "", "all", "combined"
            ]) + "\n")

        print("  Created config files")

    def generate_readme(self):
        """Generate README."""
        with open(self.output_dir / "README.md", 'w', encoding='utf-8') as f:
            f.write("""# ORFBounder Toy Example

## Structure
- 3 genes at positions 100-162, 300-362, 500-562
- Chromosome: NC_000913.3 (1000bp)
- TIS reads concentrated at start codons (positions 100, 300, 500)
- RIBO reads distributed across gene bodies
- Offset: 0 (reads start exactly at start codon)

## Run
```bash
uv run call_orfbounder.py -c toy_example/config/config.tsv -r toy_example/results
```

## Expected Results
Should detect 3 ORFs corresponding to the 3 annotated genes.
""")

    def generate_all(self):
        """Generate everything."""
        print("=" * 60)
        print("ORFBounder Toy Example Generator")
        print("=" * 60)

        genome_seq = self.generate_genome()
        self.generate_annotation()
        self.generate_bam_files(genome_seq)
        self.generate_configs()
        self.generate_readme()

        print("=" * 60)
        print(f"✓ Generated: {self.output_dir}")
        print("=" * 60)
        print(f"\nuv run call_orfbounder.py -c {self.output_dir}/config/config.tsv -r {self.output_dir}/results\n")


if __name__ == "__main__":
    ToyExampleGenerator().generate_all()