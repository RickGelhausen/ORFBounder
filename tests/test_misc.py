"""
Unit tests for miscellaneous utility functions
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from collections import deque, OrderedDict
from unittest.mock import patch, Mock
from interlap import InterLap

from lib.misc import (
    CODON_LENGTH,
    NT_WINDOW_SIZE,
    NEAR_ANNOTATED_THRESHOLD,
    CODON_INTERVAL_OFFSET,
    CODON_INTERVAL_SIZE,
    VALID_MAPPING_TYPES,
    generate_annotation_dict,
    annotation_interlap,
    calculate_density,
    create_codon_interlaps,
    get_frame,
    get_genome_information,
    get_gene_information,
    calculate_utr_distance,
    calculate_relative_density,
    generate_result_dataframe,
    dictionary_depth,
    base_mapping
)


class TestConstants:
    """Tests for module constants"""

    def test_codon_length(self):
        """Test CODON_LENGTH is 3"""
        assert CODON_LENGTH == 3

    def test_nt_window_size(self):
        """Test NT_WINDOW_SIZE is 15"""
        assert NT_WINDOW_SIZE == 15

    def test_near_annotated_threshold(self):
        """Test NEAR_ANNOTATED_THRESHOLD is 10"""
        assert NEAR_ANNOTATED_THRESHOLD == 10

    def test_codon_interval_offset(self):
        """Test CODON_INTERVAL_OFFSET is 2"""
        assert CODON_INTERVAL_OFFSET == 2

    def test_codon_interval_size(self):
        """Test CODON_INTERVAL_SIZE is 5"""
        assert CODON_INTERVAL_SIZE == 5

    def test_valid_mapping_types(self):
        """Test VALID_MAPPING_TYPES contains expected values"""
        assert "fiveprime" in VALID_MAPPING_TYPES
        assert "threeprime" in VALID_MAPPING_TYPES
        assert "global" in VALID_MAPPING_TYPES
        assert "centered" in VALID_MAPPING_TYPES


class TestGenerateAnnotationDict:
    """Tests for generate_annotation_dict function"""

    @pytest.fixture
    def sample_gff_content(self):
        """Sample GFF3 content"""
        return """##gff-version 3
chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1;locus_tag=gene1;name=protein1;gene_id=g1
chr1\tORFBounder\tgene\t100\t201\t.\t+\t.\tgene_name=protein1;locus_tag=gene1
chr2\tORFBounder\tCDS\t300\t500\t.\t-\t.\tID=cds2;locus_tag=gene2
"""

    def test_generate_annotation_dict_basic(self, tmp_path, sample_gff_content):
        """Test basic annotation parsing"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(sample_gff_content)

        result = generate_annotation_dict(gff_file)

        assert "chr1:100-201:+" in result
        assert "chr2:300-500:-" in result

    def test_generate_annotation_dict_cds_values(self, tmp_path, sample_gff_content):
        """Test CDS values are parsed correctly"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(sample_gff_content)

        result = generate_annotation_dict(gff_file)

        gene_id, locus_tag, name, read_list, gene_name, old_locus_tag = result["chr1:100-201:+"]
        assert gene_id == "g1"
        assert locus_tag == "gene1"
        assert name == "protein1"
        assert gene_name == "protein1"
        assert old_locus_tag == ""
        assert read_list == []

    def test_generate_annotation_dict_nonexistent_file(self, tmp_path):
        """Test with non-existent file raises error"""
        gff_file = tmp_path / "nonexistent.gff"

        with pytest.raises(FileNotFoundError, match=f"Annotation file does not exist: {gff_file}"):
            generate_annotation_dict(gff_file)

    def test_generate_annotation_dict_invalid_attributes(self, tmp_path):
        """Test with invalid attribute format raises error"""
        gff_content = """chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1;invalid_format
"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(gff_content)

        with pytest.raises(ValueError, match="Error: invalid gff, wrongly formatted attribute fields."):
            generate_annotation_dict(gff_file)

    def test_generate_annotation_dict_gene_cds_merge(self, tmp_path):
        """Test that gene and CDS features are merged correctly"""
        gff_content = """chr1\tORFBounder\tgene\t100\t201\t.\t+\t.\tgene_name=myprotein;locus_tag=mygene;old_locus_tag=old123
chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1
"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(gff_content)

        result = generate_annotation_dict(gff_file)

        gene_id, locus_tag, name, read_list, gene_name, old_locus_tag = result["chr1:100-201:+"]
        assert locus_tag == "mygene"
        assert gene_name == "myprotein"
        assert old_locus_tag == "old123"

    def test_generate_annotation_dict_pseudogene(self, tmp_path):
        """Test pseudogene feature handling"""
        gff_content = """chr1\tORFBounder\tpseudogene\t100\t201\t.\t+\t.\tname=pseudo1;locus_tag=ps1
chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1
"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(gff_content)

        result = generate_annotation_dict(gff_file)

        assert "chr1:100-201:+" in result
        gene_id, locus_tag, name, read_list, gene_name, old_locus_tag = result["chr1:100-201:+"]
        assert gene_name == "pseudo1"


    def test_generate_annotation_dict_alternative_fields(self, tmp_path):
        """Test alternative field names (gene_name, gene_id)"""
        gff_content = """chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tid=cds1;gene_name=altname
"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(gff_content)

        result = generate_annotation_dict(gff_file)

        gene_id, locus_tag, name, read_list, gene_name, old_locus_tag = result["chr1:100-201:+"]
        assert gene_id == "cds1"
        assert name == "altname"


class TestAnnotationInterlap:
    """Tests for annotation_interlap function"""

    @pytest.fixture
    def sample_gff_content(self):
        """Sample GFF3 content"""
        return """chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1;locus_tag=gene1
chr1\tORFBounder\tCDS\t300\t401\t.\t-\t.\tID=cds2;locus_tag=gene2
chr2\tORFBounder\tCDS\t500\t601\t.\t+\t.\tID=cds3;locus_tag=gene3
"""

    def test_annotation_interlap_creates_interlaps(self, tmp_path, sample_gff_content):
        """Test interlap objects are created"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(sample_gff_content)

        interlap_dict, gene_dict = annotation_interlap(gff_file)

        assert ("chr1", "+") in interlap_dict
        assert ("chr1", "-") in interlap_dict
        assert ("chr2", "+") in interlap_dict

    def test_annotation_interlap_gene_density_dict(self, tmp_path, sample_gff_content):
        """Test gene density dict is created correctly"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(sample_gff_content)

        interlap_dict, gene_dict = annotation_interlap(gff_file)

        assert "gene1" in gene_dict
        assert "gene2" in gene_dict
        assert "gene3" in gene_dict

        # Check structure: [chromosome, start, stop, strand, density]
        assert gene_dict["gene1"] == ["chr1", 100, 201, "+", 0]
        assert gene_dict["gene2"] == ["chr1", 300, 401, "-", 0]

    def test_annotation_interlap_find_overlaps(self, tmp_path, sample_gff_content):
        """Test that interlap can find overlapping intervals"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(sample_gff_content)

        interlap_dict, gene_dict = annotation_interlap(gff_file)

        # Find overlaps with position in gene1
        overlaps = list(interlap_dict[("chr1", "+")].find((150, 150)))

        assert len(overlaps) > 0
        assert overlaps[0][2] == "gene1"

    def test_annotation_interlap_uses_gene_name_fallback(self, tmp_path):
        """Test fallback to gene_name when locus_tag is empty"""
        gff_content = """chr1\tORFBounder\tgene\t100\t201\t.\t+\t.\tgene_name=fallback_name
chr1\tORFBounder\tCDS\t100\t201\t.\t+\t.\tID=cds1
"""
        gff_file = tmp_path / "annotation.gff"
        gff_file.write_text(gff_content)

        interlap_dict, gene_dict = annotation_interlap(gff_file)

        assert "fallback_name" in gene_dict


class TestCalculateDensity:
    """Tests for calculate_density function"""

    def test_calculate_density_basic(self):
        """Test basic density calculation"""
        alignment_position_dict = {
            ("chr1", "+"): {150: 10, 160: 20}
        }

        interlap = InterLap()
        interlap.add((100, 201, "gene1"))
        annotation_interlap_dict = {("chr1", "+"): interlap}

        gene_density_dict = {
            "gene1": ["chr1", 100, 201, "+", 0]
        }

        result = calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        # Should add 10 + 20 = 30 to gene1's density
        assert result["gene1"][4] == 30

    def test_calculate_density_no_overlap(self):
        """Test when positions don't overlap with genes"""
        alignment_position_dict = {
            ("chr1", "+"): {50: 10, 60: 20}  # Outside gene range
        }

        interlap = InterLap()
        interlap.add((100, 201, "gene1"))
        annotation_interlap_dict = {("chr1", "+"): interlap}

        gene_density_dict = {
            "gene1": ["chr1", 100, 201, "+", 0]
        }

        result = calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        # Density should remain 0
        assert result["gene1"][4] == 0

    def test_calculate_density_multiple_genes(self):
        """Test with multiple overlapping genes"""
        alignment_position_dict = {
            ("chr1", "+"): {150: 10}
        }

        interlap = InterLap()
        interlap.add((100, 200, "gene1"))
        interlap.add((140, 250, "gene2"))
        annotation_interlap_dict = {("chr1", "+"): interlap}

        gene_density_dict = {
            "gene1": ["chr1", 100, 200, "+", 0],
            "gene2": ["chr1", 140, 250, "+", 0]
        }

        result = calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        # Position 150 overlaps both genes
        assert result["gene1"][4] == 10
        assert result["gene2"][4] == 10

    def test_calculate_density_missing_chromosome(self):
        """Test when chromosome is not in annotation"""
        alignment_position_dict = {
            ("chr2", "+"): {150: 10}
        }

        interlap = InterLap()
        interlap.add((100, 200, "gene1"))
        annotation_interlap_dict = {("chr1", "+"): interlap}

        gene_density_dict = {
            "gene1": ["chr1", 100, 200, "+", 0]
        }

        result = calculate_density(alignment_position_dict, annotation_interlap_dict, gene_density_dict)

        # No change since chr2 not in annotation
        assert result["gene1"][4] == 0


class TestCreateCodonInterlapsBasic:
    """Basic functionality tests"""

    def test_empty_sequence(self):
        """Test with empty genome sequence"""
        chrom = "chr1"
        genome_seq = ""
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(interlap_dict) == 0
        assert len(codon_dict) == 0

    def test_sequence_shorter_than_codon(self):
        """Test with sequence shorter than codon length"""
        chrom = "chr1"
        genome_seq = "AT"  # Only 2 nucleotides
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(interlap_dict) == 0
        assert len(codon_dict) == 0

    def test_exact_codon_length_sequence(self):
        """Test with sequence exactly one codon long"""
        chrom = "chr1"
        genome_seq = "ATG"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should be filtered due to boundary constraints
        # interval_start = 0 - 2 = -2 (negative, filtered)
        assert len(codon_dict) == 0

    def test_no_matching_codons(self):
        """Test when no codons match"""
        chrom = "chr1"
        genome_seq = "AAAAAAAAAAAA"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(interlap_dict) == 0
        assert len(codon_dict) == 0


class TestCreateCodonInterlapsForwardStrand:
    """Tests for forward strand codon detection"""

    def test_single_atg_forward(self):
        """Test single ATG on forward strand"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"  # ATG at position 4
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert ("chr1", "+") in interlap_dict
        assert len(codon_dict) == 1

        # Check interval boundaries
        # interval_start = pos - CODON_INTERVAL_OFFSET = 4 - 2 = 2
        # interval_stop = pos + CODON_INTERVAL_OFFSET = 4 + 2 = 6
        key = "chr1:2-6:+"
        assert key in codon_dict
        assert codon_dict[key] == ["ATG", 0]

    def test_multiple_atg_forward(self):
        """Test multiple ATGs on forward strand"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNATGNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find 2 ATGs
        assert len(codon_dict) == 2

        # First ATG at position 4
        assert "chr1:2-6:+" in codon_dict
        # Second ATG at position 10
        assert "chr1:8-12:+" in codon_dict

    def test_consecutive_atg_forward(self):
        """Test consecutive ATGs"""
        chrom = "chr1"
        genome_seq = "NNNATGATGNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # ATG at position 3 and 6
        assert len(codon_dict) == 2
        assert "chr1:1-5:+" in codon_dict
        assert "chr1:4-8:+" in codon_dict

    def test_overlapping_start_codons_forward(self):
        """Test overlapping start codons (ATGTG contains both ATG and GTG)"""
        chrom = "chr1"
        genome_seq = "NNNATGTGNN"  # ATG at 3, GTG at 5
        codons = ["ATG", "GTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find both
        assert len(codon_dict) == 2
        assert "chr1:1-5:+" in codon_dict  # ATG at pos 3
        assert "chr1:3-7:+" in codon_dict  # GTG at pos 5

    def test_forward_at_sequence_start(self):
        """Test ATG at very start of sequence"""
        chrom = "chr1"
        genome_seq = "ATGNNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 0 - 2 = -2 (negative, should be filtered)
        assert len(codon_dict) == 0

    def test_forward_at_sequence_start_one(self):
        """Test ATG at very start of sequence"""
        chrom = "chr1"
        genome_seq = "NATGNNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 1 - 2 = -1 (negative, should be filtered)
        assert len(codon_dict) == 0

    def test_forward_near_sequence_start(self):
        """Test ATG near start but not filtered"""
        chrom = "chr1"
        genome_seq = "NNATGNNNNN"  # ATG at position 2
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 2 - 2 = 0 (valid)
        # interval_stop = 2 + 2 = 4 (valid)
        assert len(codon_dict) == 1
        assert "chr1:0-4:+" in codon_dict

    def test_forward_at_sequence_end(self):
        """Test ATG at end of sequence"""
        chrom = "chr1"
        genome_seq = "NNNNATG"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # ATG at position 4
        # interval_stop = 4 + 2 = 6
        # len(genome_seq) - CODON_LENGTH + 1 = 7 - 3 + 1 = 5
        # 6 > 5, should be filtered
        assert len(codon_dict) == 0

    def test_forward_near_sequence_end(self):
        """Test ATG near end but not filtered"""
        chrom = "chr1"
        genome_seq = "NNATGNN"  # ATG at position 2
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_stop = 2 + 2 = 4
        # len(genome_seq) - CODON_LENGTH + 1 = 7 - 3 + 1 = 5
        # 4 <= 5, valid
        assert len(codon_dict) == 1


class TestCreateCodonInterlapsReverseStrand:
    """Tests for reverse strand codon detection"""

    def test_single_cat_reverse(self):
        """Test single CAT (reverse complement of ATG) on reverse strand"""
        chrom = "chr1"
        genome_seq = "NNNNCATNNNNN"  # CAT at position 4
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert ("chr1", "-") in interlap_dict
        assert len(codon_dict) == 1

        # For reverse strand:
        # interval_start = pos = 4
        # interval_stop = pos + CODON_INTERVAL_SIZE - 1 = 4 + 5 - 1 = 8
        key = "chr1:4-8:-"
        assert key in codon_dict
        assert codon_dict[key] == ["ATG", 0]  # Stored as original codon

    def test_multiple_cat_reverse(self):
        """Test multiple CATs on reverse strand"""
        chrom = "chr1"
        genome_seq = "NNNNCATNNNCATNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find 2 reverse complement matches
        assert len(codon_dict) == 2
        assert "chr1:4-8:-" in codon_dict
        assert "chr1:10-14:-" in codon_dict

    def test_reverse_gtg_tac(self):
        """Test GTG reverse complement (CAC)"""
        chrom = "chr1"
        genome_seq = "NNNNCACNNNNN"  # CAC is reverse complement of GTG
        codons = ["GTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(codon_dict) == 1
        key = "chr1:4-8:-"
        assert key in codon_dict
        assert codon_dict[key] == ["GTG", 0]

    def test_reverse_ttg_caa(self):
        """Test TTG reverse complement (CAA)"""
        chrom = "chr1"
        genome_seq = "NNNNCAANNNNN"  # CAA is reverse complement of TTG
        codons = ["TTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(codon_dict) == 1
        assert codon_dict["chr1:4-8:-"] == ["TTG", 0]

    def test_reverse_at_sequence_start(self):
        """Test reverse codon at sequence start"""
        chrom = "chr1"
        genome_seq = "CATNNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 0
        # interval_stop = 0 + 5 - 1 = 4
        # Valid, should be found
        assert len(codon_dict) == 1
        assert "chr1:0-4:-" in codon_dict

    def test_reverse_at_sequence_end(self):
        """Test reverse codon at sequence end"""
        chrom = "chr1"
        genome_seq = "NNNNCAT"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # CAT at position 4
        # interval_stop = 4 + 5 - 1 = 8
        # len(genome_seq) - CODON_LENGTH + 1 = 7 - 3 + 1 = 5
        # 8 > 5, should be filtered
        assert len(codon_dict) == 0

    def test_reverse_near_sequence_end(self):
        """Test reverse codon near end but valid"""
        chrom = "chr1"
        genome_seq = "CATNN"  # CAT at position 0
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_stop = 0 + 5 - 1 = 4
        # len(genome_seq)  = 5
        # 4 <= 5, should be filtered
        assert len(codon_dict) == 1

        genome_seq = "CATN"  # CAT at position 0
        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)
        # interval_stop = 0 + 5 - 1 = 4
        # len(genome_seq)  = 4
        # 4 <= 4, should be filtered
        assert len(codon_dict) == 0

class TestCreateCodonInterlapsBothStrands:
    """Tests for scenarios with both forward and reverse strand codons"""

    def test_both_strands_separate(self):
        """Test codons on both strands in different locations"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNCATNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert ("chr1", "+") in interlap_dict
        assert ("chr1", "-") in interlap_dict
        assert len(codon_dict) == 2

    def test_both_strands_overlapping(self):
        """Test when forward and reverse codons overlap"""
        chrom = "chr1"
        genome_seq = "NNATGCATNN"  # ATG at 2, CAT at 5
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Both should be found
        assert len(codon_dict) == 2
        assert "chr1:0-4:+" in codon_dict  # ATG forward
        assert "chr1:5-9:-" in codon_dict  # CAT reverse


class TestCreateCodonInterlapsMultipleCodonTypes:
    """Tests with multiple codon types"""

    def test_all_three_start_codons(self):
        """Test with ATG, GTG, and TTG"""
        chrom = "chr1"
        genome_seq = "NNATGNNGTGNNTTGNN"
        codons = ["ATG", "GTG", "TTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert len(codon_dict) == 3

        # Verify each codon is stored with correct type
        codon_types = set(v[0] for v in codon_dict.values())
        assert "ATG" in codon_types
        assert "GTG" in codon_types
        assert "TTG" in codon_types

    def test_mixed_strands_multiple_codons(self):
        """Test multiple codon types on both strands"""
        chrom = "chr1"
        genome_seq = "NNATGNCACNGTGNCAANTTGNN"
        # ATG at 2, CAC (GTG rev) at 6, GTG at 10, CAA (TTG rev) at 14, TTG at 18
        codons = ["ATG", "GTG", "TTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find all 5 codons
        assert len(codon_dict) == 5

        # Check both strands are represented
        assert ("chr1", "+") in interlap_dict
        assert ("chr1", "-") in interlap_dict


class TestCreateCodonInterlapsInterlapStructure:
    """Tests for interlap data structure correctness"""

    def test_interlap_contains_correct_intervals(self):
        """Test that interlap contains correct interval boundaries"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Get the interlap for positive strand
        interlap = interlap_dict[("chr1", "+")]

        # Query the interval
        overlaps = list(interlap.find((3, 3)))

        assert len(overlaps) == 1
        start, stop, key = overlaps[0]
        assert start == 2
        assert stop == 6
        assert key == "chr1:2-6:+"

    def test_interlap_can_find_overlaps(self):
        """Test that interlap can find overlapping intervals"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        interlap = interlap_dict[("chr1", "+")]

        # Query should find the interval
        assert len(list(interlap.find((4, 4)))) == 1  # Center of ATG
        assert len(list(interlap.find((2, 2)))) == 1  # Start of interval
        assert len(list(interlap.find((6, 6)))) == 1  # End of interval
        assert len(list(interlap.find((0, 1)))) == 0  # Before interval
        assert len(list(interlap.find((7, 8)))) == 0  # After interval

    def test_interlap_key_matches_codon_dict(self):
        """Test that interlap key matches codon_dict key"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        interlap = interlap_dict[("chr1", "+")]
        overlaps = list(interlap.find((4, 4)))

        key = overlaps[0][2]
        assert key in codon_dict

    def test_multiple_codons_in_same_interlap(self):
        """Test that multiple codons on same strand go into same interlap"""
        chrom = "chr1"
        genome_seq = "NNATGNNNATGNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should only have one interlap object for positive strand
        assert len(interlap_dict) == 1
        assert ("chr1", "+") in interlap_dict

        # But should have 2 entries in codon_dict
        assert len(codon_dict) == 2


class TestCreateCodonInterlapsCodonDict:
    """Tests for codon_dict structure and contents"""

    def test_codon_dict_key_format(self):
        """Test codon_dict key format is correct"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        for key in codon_dict:
            # Key should be "chrom:start-stop:strand"
            parts = key.split(":")
            assert len(parts) == 3
            assert parts[0] == chrom
            assert "-" in parts[1]
            assert parts[2] in ["+", "-"]

    def test_codon_dict_value_structure(self):
        """Test codon_dict value is [codon, count]"""
        chrom = "chr1"
        genome_seq = "NNNNATGNNNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        for value in codon_dict.values():
            assert isinstance(value, list)
            assert len(value) == 2
            assert isinstance(value[0], str)  # Codon sequence
            assert value[0] in codons  # Must be one of the input codons
            assert value[1] == 0  # Initial count

    def test_codon_dict_stores_original_codon(self):
        """Test that reverse complement codons are stored as original"""
        chrom = "chr1"
        genome_seq = "NNNNCATNNNNN"  # CAT = reverse complement of ATG
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should store as ATG, not CAT
        for value in codon_dict.values():
            assert value[0] == "ATG"
            assert value[0] != "CAT"


class TestCreateCodonInterlapsBoundaryConditions:
    """Tests for boundary conditions and edge cases"""

    def test_boundary_offset_equals_position(self):
        """Test when CODON_INTERVAL_OFFSET equals position"""
        chrom = "chr1"
        # ATG at position 2 (CODON_INTERVAL_OFFSET = 2)
        genome_seq = "NNATGNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 2 - 2 = 0 (valid)
        assert len(codon_dict) == 1
        assert "chr1:0-4:+" in codon_dict

    def test_boundary_offset_minus_one(self):
        """Test just before boundary"""
        chrom = "chr1"
        # ATG at position 1
        genome_seq = "NATGNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # interval_start = 1 - 2 = -1 (negative, filtered)
        assert len(codon_dict) == 0

    def test_very_long_sequence(self):
        """Test with long sequence to ensure no performance issues"""
        chrom = "chr1"
        # Create a long sequence with ATGs scattered throughout
        genome_seq = "N" * 100 + "ATG" + "N" * 100 + "ATG" + "N" * 100
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find both ATGs
        assert len(codon_dict) == 2

    def test_many_consecutive_codons(self):
        """Test with many consecutive codons"""
        chrom = "chr1"
        genome_seq = "NNATGATGATGATGNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find 4 ATGs
        assert len(codon_dict) == 4


class TestCreateCodonInterlapsRealWorldScenarios:
    """Tests simulating real-world usage"""

    def test_ecoli_start_codons(self):
        """Test with E. coli typical start codons"""
        chrom = "NC_000913"
        genome_seq = "NNATGNNGTGNNCAANNTTGNN"
        codons = ["ATG", "GTG", "TTG"]  # E. coli uses all three

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find all 4 forward codons
        assert len([k for k in codon_dict.keys() if ":+" in k]) == 3

    def test_overlapping_intervals(self):
        """Test when interval windows overlap"""
        chrom = "chr1"
        genome_seq = "NNATGATGNN"  # ATG at 2 and 5, intervals will overlap
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Both should be found even though intervals overlap
        assert len(codon_dict) == 2

    def test_realistic_gene_start(self):
        """Test realistic gene start region"""
        chrom = "chr1"
        # Typical sequence around start codon with Shine-Dalgarno upstream
        genome_seq = "AGGAGGTGATGAAAAAANNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Should find the ATG
        assert len(codon_dict) >= 1

        # Check ATG is found
        found_atg = any("ATG" == v[0] for v in codon_dict.values())
        assert found_atg


class TestCreateCodonInterlapsInvariantsAndProperties:
    """Tests for invariants and properties that should always hold"""

    def test_codon_dict_and_interlap_consistency(self):
        """Test that every codon_dict entry has corresponding interlap entry"""
        chrom = "chr1"
        genome_seq = "NNATGNNGTGNCATNCAANN"
        codons = ["ATG", "GTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Every key in codon_dict should be findable in interlap_dict
        for key in codon_dict.keys():
            chrom_key, interval, strand = key.split(":")
            start, stop = map(int, interval.split("-"))

            # Should be able to find it in the interlap
            overlaps = list(interlap_dict[(chrom_key, strand)].find((start, stop)))
            assert len(overlaps) >= 1

            # One of the overlaps should have this exact key
            keys_found = [o[2] for o in overlaps]
            assert key in keys_found

    def test_no_negative_positions(self):
        """Test that no intervals have negative positions"""
        chrom = "chr1"
        genome_seq = "ATGATGCATCATNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        # Check all keys have non-negative positions
        for key in codon_dict.keys():
            _, interval, _ = key.split(":")
            start, stop = map(int, interval.split("-"))
            assert start >= 0
            assert stop >= 0

    def test_intervals_respect_size_constraints(self):
        """Test that interval sizes are as expected"""
        chrom = "chr1"
        genome_seq = "NNATGNNGTGNCATNCAANN"
        codons = ["ATG", "GTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        for key in codon_dict.keys():
            _, interval, strand = key.split(":")
            start, stop = map(int, interval.split("-"))

            if strand == "+":
                # Forward: interval_stop - interval_start should equal 2*CODON_INTERVAL_OFFSET
                assert stop - start == 2 * CODON_INTERVAL_OFFSET
            else:
                # Reverse: should equal CODON_INTERVAL_SIZE - 1
                assert stop - start == CODON_INTERVAL_SIZE - 1

    def test_codon_count_always_zero(self):
        """Test that initial count is always 0"""
        chrom = "chr1"
        genome_seq = "NNATGNNGTGNCATNN"
        codons = ["ATG", "GTG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        for value in codon_dict.values():
            assert value[1] == 0

    def test_return_types(self):
        """Test that return types are correct"""
        chrom = "chr1"
        genome_seq = "NNATGNNN"
        codons = ["ATG"]

        interlap_dict, codon_dict = create_codon_interlaps(chrom, genome_seq, codons)

        assert isinstance(interlap_dict, dict)
        assert isinstance(codon_dict, dict)

        # Check interlap_dict keys are tuples
        for key in interlap_dict.keys():
            assert isinstance(key, tuple)
            assert len(key) == 2
            assert isinstance(key[0], str)  # chromosome
            assert key[1] in ["+", "-"]  # strand

        # Check interlap values are InterLap objects
        for value in interlap_dict.values():
            assert isinstance(value, InterLap)

class TestGetFrame:
    """Tests for get_frame function"""

    def test_get_frame_zero(self):
        """Test frame 0"""
        assert get_frame(0) == 0
        assert get_frame(3) == 0
        assert get_frame(6) == 0

    def test_get_frame_one(self):
        """Test frame 1"""
        assert get_frame(1) == 1
        assert get_frame(4) == 1
        assert get_frame(7) == 1

    def test_get_frame_two(self):
        """Test frame 2"""
        assert get_frame(2) == 2
        assert get_frame(5) == 2
        assert get_frame(8) == 2

    def test_get_frame_large_position(self):
        """Test with large position"""
        assert get_frame(100) == 1
        assert get_frame(999) == 0
        assert get_frame(1000) == 1


class TestGetGenomeInformation:
    """Tests for get_genome_information function"""

    def test_get_genome_information_plus_strand(self):
        """Test extraction on plus strand"""
        genome_seq = "A" * 20 + "ATGCCCCCCTAA" + "G" * 20
        start = 20
        stop = 31
        strand = "+"

        nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(
            start, stop, strand, genome_seq
        )

        assert nt_seq == "ATGCCCCCCTAA"
        assert start_codon == "ATG"
        assert stop_codon == "TAA"
        assert len(nt_window) == 15
        assert aa_seq == "MPP*"

    def test_get_genome_information_minus_strand(self):
        """Test extraction on minus strand"""
        genome_seq = "A" * 20 + "TTAGGGGGGCAT" + "G" * 20
        start = 20
        stop = 31
        strand = "-"

        nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(
            start, stop, strand, genome_seq
        )

        # Should be reverse complemented
        assert start_codon == "ATG"
        assert stop_codon == "TAA"
        assert len(nt_window) == 15
        assert nt_seq == "ATGCCCCCCTAA"
        assert aa_seq == "MPP*"

    def test_get_genome_information_translation(self):
        """Test amino acid translation"""
        genome_seq = "ATGAAATAA"  # ATG (M), AAA (K), TAA (stop)
        start = 0
        stop = 8
        strand = "+"

        nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(
            start, stop, strand, genome_seq
        )

        assert nt_seq == "ATGAAATAA"
        assert aa_seq == "MK*"

    def test_get_genome_information_window_size(self):
        """Test that 15nt window is correct size"""
        genome_seq = "A" * 20 + "ATGCCC" + "G" * 20
        start = 20
        stop = 25
        strand = "+"

        nt_seq, aa_seq, nt_window, start_codon, stop_codon = get_genome_information(
            start, stop, strand, genome_seq
        )

        assert len(nt_window) == 15


class TestGetGeneInformation:
    """Tests for get_gene_information function"""

    @pytest.fixture
    def sample_gene_dict(self):
        """Sample gene dictionary"""
        return {
            "gene1": ("chr1", 101, 200, "+", 0),
            "gene2": ("chr1", 301, 500, "-", 0),
            "gene3": ("chr2", 100, 300, "+", 0)
        }

    def test_get_gene_information_annotated(self, sample_gene_dict):
        """Test exact match returns Annotated"""
        gene_type, gene_name = get_gene_information(
            "chr1", 100, 199, "+", sample_gene_dict
        )

        assert gene_type == "Annotated"
        assert gene_name == "gene1"

    def test_get_gene_information_unannotated(self, sample_gene_dict):
        """Test no match returns Unannotated"""
        gene_type, gene_name = get_gene_information(
            "chr1", 600, 700, "+", sample_gene_dict
        )

        assert gene_type == "Unannotated"
        assert gene_name == "chr1:601-701:+"  # 1-based

    def test_get_gene_information_near_annotated_plus(self, sample_gene_dict):
        """Test near-annotated on plus strand"""
        # Start within 10bp, same stop
        gene_type, gene_name = get_gene_information(
            "chr1", 105, 199, "+", sample_gene_dict
        )

        assert gene_type == "Near_Annotated"
        assert gene_name == "gene1"

    def test_get_gene_information_near_annotated_minus(self, sample_gene_dict):
        """Test near-annotated on minus strand"""
        # Stop within 10bp, same start
        gene_type, gene_name = get_gene_information(
            "chr1", 300, 495, "-", sample_gene_dict
        )

        assert gene_type == "Near_Annotated"
        assert gene_name == "gene2"

    def test_get_gene_information_internal_inframe_plus(self, sample_gene_dict):
        """Test internal in-frame on plus strand"""
        # Start > gene start, stop <= gene stop, same frame
        gene_type, gene_name = get_gene_information(
            "chr1", 103, 190, "+", sample_gene_dict
        )

        # 104 % 3 == 101 % 3 (both = 2), so in-frame
        assert gene_type == "Internal_Inframe"
        assert gene_name == "gene1"

    def test_get_gene_information_n_terminal_extension_plus(self, sample_gene_dict):
        """Test N-terminal extension on plus strand"""
        # Start < gene start, same stop
        gene_type, gene_name = get_gene_information(
            "chr1", 80, 199, "+", sample_gene_dict
        )

        assert gene_type == "N-terminal_extension"
        assert gene_name == "gene1"

    def test_get_gene_information_n_terminal_extension_minus(self, sample_gene_dict):
        """Test N-terminal extension on minus strand"""
        # Stop > gene stop, same start
        gene_type, gene_name = get_gene_information(
            "chr1", 300, 520, "-", sample_gene_dict
        )

        assert gene_type == "N-terminal_extension"
        assert gene_name == "gene2"

    def test_get_gene_information_internal_outofframe(self, sample_gene_dict):
        """Test internal out-of-frame"""
        # Within gene bounds but different frame
        gene_type, gene_name = get_gene_information(
            "chr1", 102, 190, "+", sample_gene_dict
        )

        # 103 % 3 != 101 % 3, so out of frame
        assert gene_type == "Internal_OutofFrame"
        assert gene_name == "gene1"

    def test_get_gene_information_wrong_chromosome(self, sample_gene_dict):
        """Test with wrong chromosome"""
        gene_type, gene_name = get_gene_information(
            "chr3", 100, 199, "+", sample_gene_dict
        )

        assert gene_type == "Unannotated"

    def test_get_gene_information_wrong_strand(self, sample_gene_dict):
        """Test with wrong strand"""
        gene_type, gene_name = get_gene_information(
            "chr1", 100, 199, "-", sample_gene_dict
        )

        assert gene_type == "Unannotated"


class TestCalculateUtrDistance:
    """Tests for calculate_utr_distance function"""

    @pytest.fixture
    def sample_gene_dict(self):
        """Sample gene dictionary"""
        return {
            "gene1": ("chr1", 100, 300, "+", 0)
        }

    def test_calculate_utr_distance_basic(self, sample_gene_dict):
        """Test basic UTR distance calculation"""
        fiveprime_dist, threeprime_dist = calculate_utr_distance(
            99, 299, "gene1", sample_gene_dict  # 0-based input
        )

        # 100 (1-based) - 100 (gene start) = 0
        # 300 (gene stop) - 100 (1-based) = 200
        assert fiveprime_dist == 0
        assert threeprime_dist == 200

    def test_calculate_utr_distance_upstream(self, sample_gene_dict):
        """Test when ORF starts upstream of gene"""
        fiveprime_dist, threeprime_dist = calculate_utr_distance(
            79, 299, "gene1", sample_gene_dict
        )

        # 80 - 100 = -20
        assert fiveprime_dist == -20

    def test_calculate_utr_distance_downstream(self, sample_gene_dict):
        """Test when ORF starts downstream"""
        fiveprime_dist, threeprime_dist = calculate_utr_distance(
            149, 299, "gene1", sample_gene_dict
        )

        # 150 - 100 = 50
        assert fiveprime_dist == 50

    def test_calculate_utr_distance_gene_not_in_dict(self):
        """Test with gene not in dictionary"""
        gene_dict = {}

        fiveprime_dist, threeprime_dist = calculate_utr_distance(
            99, 299, "unknown_gene", gene_dict
        )

        assert np.isnan(fiveprime_dist)
        assert np.isnan(threeprime_dist)


class TestCalculateRelativeDensity:
    """Tests for calculate_relative_density function"""

    @pytest.fixture
    def sample_gene_dict(self):
        """Sample gene dictionary with densities"""
        return {
            "gene1": ("chr1", 100, 300, "+", 100.0),  # Gene RPM = 100
            "gene2": ("chr1", 400, 600, "+", 0.0)     # Gene RPM = 0
        }

    def test_calculate_relative_density_basic(self, sample_gene_dict):
        """Test basic relative density calculation"""
        density = calculate_relative_density(
            50.0, "gene1", "Annotated", sample_gene_dict
        )

        # 50 / 100 = 0.5
        assert density == 0.5

    def test_calculate_relative_density_gene_rpm_zero(self, sample_gene_dict):
        """Test when gene RPM is zero"""
        density = calculate_relative_density(
            50.0, "gene2", "Annotated", sample_gene_dict
        )

        assert density == 0.0

    def test_calculate_relative_density_n_terminal_returns_nan(self, sample_gene_dict):
        """Test N-terminal extension returns NaN"""
        density = calculate_relative_density(
            50.0, "gene1", "N-terminal_extension", sample_gene_dict
        )

        assert np.isnan(density)

    def test_calculate_relative_density_empty_dict(self):
        """Test with empty gene dict returns NaN"""
        density = calculate_relative_density(
            50.0, "gene1", "Annotated", {}
        )

        assert np.isnan(density)

    def test_calculate_relative_density_gene_not_in_dict(self, sample_gene_dict):
        """Test gene not in dictionary returns NaN"""
        density = calculate_relative_density(
            50.0, "unknown_gene", "Annotated", sample_gene_dict
        )

        assert np.isnan(density)


class TestGenerateResultDataframe:
    """Tests for generate_result_dataframe function"""

    @pytest.fixture
    def sample_data(self):
        """Sample data for dataframe generation"""
        detected_orfs_dict = {
            ("chr1", "+"): {
                (99, 200): (100.0, 50.0, 75.0)  # rpm_tis, rpm_tts, rpm_ribo
            }
        }

        gene_dict = {
            "gene1": ("chr1", 100, 201, "+", 100.0)
        }

        genome = {
            "chr1": "A" * 100 + "ATGAAATAA" + "G" * 100
        }

        read_count_dict = {
            ("chr1", 99, 200, "+"): [10, 20]
        }

        accepted_read_list = [
            {"chr1": 1000},
            {"chr1": 2000}
        ]

        wildcards = ["TIS-control-rep1", "RNATIS-control-rep1"]

        headers = ("TIS", "TTS", "RIBO")

        return (detected_orfs_dict, gene_dict, gene_dict, gene_dict,
                genome, read_count_dict, accepted_read_list, wildcards, headers)

    def test_generate_result_dataframe_basic(self, sample_data):
        """Test basic dataframe generation"""
        (detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
         genome, read_count_dict, accepted_read_list, wildcards, headers) = sample_data

        df = generate_result_dataframe(
            detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
            genome, read_count_dict, accepted_read_list, wildcards, headers
        )

        assert len(df) == 1
        assert "Type" in df.columns
        assert "Identifier" in df.columns

    def test_generate_result_dataframe_columns(self, sample_data):
        """Test expected columns are present"""
        (detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
         genome, read_count_dict, accepted_read_list, wildcards, headers) = sample_data

        df = generate_result_dataframe(
            detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
            genome, read_count_dict, accepted_read_list, wildcards, headers
        )

        expected_columns = ["Type", "Identifier", "Genome", "Start", "Stop", "Strand",
                          "Locus_tag", "Codon_count", "Start_codon", "Stop_codon",
                          "Nucleotide_Seq", "Amino_Acid_Seq"]

        for col in expected_columns:
            assert col in df.columns

    def test_generate_result_dataframe_filters_nan_orfs(self):
        """Test that ORFs with both NaN TIS and TTS are filtered"""
        detected_orfs_dict = {
            ("chr1", "+"): {
                (99, 200): (np.nan, np.nan, 75.0)  # Both TIS and TTS are NaN
            }
        }

        genome = {"chr1": "A" * 100 + "ATGAAATAA" + "G" * 100}
        gene_dict = {}

        df = generate_result_dataframe(
            detected_orfs_dict, gene_dict, gene_dict, gene_dict,
            genome, {}, [], [], ("TIS", "TTS", "RIBO")
        )

        assert len(df) == 0

    def test_generate_result_dataframe_filters_multiple_stops(self):
        """Test that ORFs with multiple stop codons are filtered"""
        detected_orfs_dict = {
            ("chr1", "+"): {
                (100, 111): (100.0, 50.0, 75.0)
            }
        }

        # Sequence with multiple stop codons
        genome = {"chr1": "A" * 100 + "ATGTAATAATAA" + "G" * 100}
        gene_dict = {}

        df = generate_result_dataframe(
            detected_orfs_dict, gene_dict, gene_dict, gene_dict,
            genome, {}, [], [], ("TIS", "TTS", "RIBO")
        )

        # Should be filtered due to multiple stops
        assert len(df) == 0

    def test_generate_result_dataframe_sorted(self, sample_data):
        """Test that dataframe is sorted correctly"""
        (detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
         genome, read_count_dict, accepted_read_list, wildcards, headers) = sample_data

        # Add another ORF at earlier position
        detected_orfs_dict[("chr1", "+")][(49, 150)] = (100.0, 50.0, 75.0)

        read_count_dict[("chr1", 49, 150, "+")] = [10, 20]
        df = generate_result_dataframe(
            detected_orfs_dict, gene_tis_dict, gene_tts_dict, gene_ribo_dict,
            genome, read_count_dict, accepted_read_list, wildcards, headers
        )

        # Should be sorted by genome, start, stop, strand
        assert df.iloc[0]["Start"] < df.iloc[1]["Start"]


class TestDictionaryDepth:
    """Tests for dictionary_depth function"""

    def test_dictionary_depth_flat(self):
        """Test depth of flat dictionary"""
        d = {"a": 1, "b": 2, "c": 3}

        assert dictionary_depth(d) == 1

    def test_dictionary_depth_nested_one_level(self):
        """Test depth with one level nesting"""
        d = {"a": {"b": 1}}

        assert dictionary_depth(d) == 2

    def test_dictionary_depth_nested_two_levels(self):
        """Test depth with two levels nesting"""
        d = {"a": {"b": {"c": 1}}}

        assert dictionary_depth(d) == 3

    def test_dictionary_depth_mixed(self):
        """Test depth with mixed nesting"""
        d = {"a": 1, "b": {"c": {"d": 2}}, "e": 3}

        # Maximum depth is 3
        assert dictionary_depth(d) == 3

    def test_dictionary_depth_empty(self):
        """Test depth of empty dictionary"""
        d = {}

        assert dictionary_depth(d) == 0


class TestBaseMapping:
    """Tests for base_mapping function"""

    def test_base_mapping_fiveprime(self):
        """Test fiveprime mapping extraction"""
        assert base_mapping("fiveprime") == "fiveprime"
        assert base_mapping("fiveprime_offset") == "fiveprime"
        assert base_mapping("some_fiveprime_suffix") == "fiveprime"

    def test_base_mapping_threeprime(self):
        """Test threeprime mapping extraction"""
        assert base_mapping("threeprime") == "threeprime"
        assert base_mapping("threeprime_offset") == "threeprime"

    def test_base_mapping_global(self):
        """Test global mapping extraction"""
        assert base_mapping("global") == "global"
        assert base_mapping("global_coverage") == "global"

    def test_base_mapping_centered(self):
        """Test centered mapping extraction"""
        assert base_mapping("centered") == "centered"
        assert base_mapping("centered_mode") == "centered"

    def test_base_mapping_invalid(self):
        """Test invalid mapping raises error"""
        with pytest.raises(ValueError) as exc_info:
            base_mapping("invalid_mapping")

        assert "Unknown mapping type" in str(exc_info.value)
        assert "invalid_mapping" in str(exc_info.value)

    def test_base_mapping_priority(self):
        """Test that first matching type is returned"""
        # If string contains multiple types, first match wins
        result = base_mapping("fiveprime_threeprime")
        assert result in VALID_MAPPING_TYPES