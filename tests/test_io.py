"""
Unit tests for utility functions in lib.io module
"""

import re
from unittest.mock import patch
import pytest
import json
import csv
from pathlib import Path
import pandas as pd
from io import StringIO

from lib.io import (
    generate_genome_dict, parse_read_lengths, parse_total_reads,
    parse_alignment_input, check_alignment_path_input, parse_offset_json,
    write_gff_file, write_codon_interval_gff, excel_writer,
    write_results_to_gff, write_results_to_table
)


class TestGenerateGenomeDict:
    """Tests for generate_genome_dict function"""

    def test_generate_genome_dict_single_sequence(self, tmp_path):
        """Test parsing single sequence FASTA file"""
        fasta_content = """>chr1
ATGCGATCGATCG
ATCGATCGATCGA
"""
        fasta_file = tmp_path / "genome.fasta"
        fasta_file.write_text(fasta_content)

        result = generate_genome_dict(fasta_file)

        assert "chr1" in result
        assert result["chr1"] == "ATGCGATCGATCGATCGATCGATCGA"
        assert len(result) == 1

    def test_generate_genome_dict_multiple_sequences(self, tmp_path):
        """Test parsing multiple sequences"""
        fasta_content = """>chr1
ATGCGATCG
>chr2
GCTAGCTAG
>chr3
TTTTTTTT
"""
        fasta_file = tmp_path / "genome.fasta"
        fasta_file.write_text(fasta_content)

        result = generate_genome_dict(fasta_file)

        assert len(result) == 3
        assert result["chr1"] == "ATGCGATCG"
        assert result["chr2"] == "GCTAGCTAG"
        assert result["chr3"] == "TTTTTTTT"

    def test_generate_genome_dict_empty_file(self, tmp_path):
        """Test with empty FASTA file returns empty dict"""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        result = generate_genome_dict(fasta_file)

        assert isinstance(result, dict)
        assert len(result) == 0

    def test_generate_genome_dict_multiline_sequences(self, tmp_path):
        """Test sequences split across multiple lines"""
        fasta_content = """>chr1
ATGC
GATC
GATC
>chr2
GCTA
"""
        fasta_file = tmp_path / "genome.fasta"
        fasta_file.write_text(fasta_content)

        result = generate_genome_dict(fasta_file)

        assert result["chr1"] == "ATGCGATCGATC"
        assert result["chr2"] == "GCTA"


class TestParseReadLengths:
    """Tests for parse_read_lengths function"""

    def test_parse_read_lengths_single_int(self, tmp_path):
        """Test parsing single integer value"""
        json_data = {"sample1": 30, "default": 28}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        assert result["sample1"] == ["30"]
        assert result["default"] == ["28"]

    def test_parse_read_lengths_comma_separated(self, tmp_path):
        """Test parsing comma-separated values"""
        json_data = {"sample1": "28,29,30", "default": "25,26"}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        assert "28" in result["sample1"]
        assert "29" in result["sample1"]
        assert "30" in result["sample1"]
        assert len(result["sample1"]) == 3
        assert "25" in result["default"]
        assert "26" in result["default"]

    def test_parse_read_lengths_range(self, tmp_path):
        """Test parsing range notation"""
        json_data = {"sample1": "28-32"}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        expected = ["28", "29", "30", "31", "32"]
        assert len(result["sample1"]) == 5
        for val in expected:
            assert val in result["sample1"]

    def test_parse_read_lengths_reverse_range(self, tmp_path):
        """Test range with reversed numbers (should handle correctly)"""
        json_data = {"sample1": "32-28"}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        # Should swap and expand correctly
        expected = ["28", "29", "30", "31", "32"]
        assert len(result["sample1"]) == 5
        for val in expected:
            assert val in result["sample1"]

    def test_parse_read_lengths_mixed_format(self, tmp_path):
        """Test mixed comma and range format"""
        json_data = {"sample1": "25,28-30,35"}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        expected = ["25", "28", "29", "30", "35"]
        assert len(result["sample1"]) == 5
        for val in expected:
            assert val in result["sample1"]

    def test_parse_read_lengths_empty_file(self, tmp_path):
        """Test with empty file returns None"""
        json_file = tmp_path / "empty.json"

        result = parse_read_lengths(json_file)

        assert result is None

    def test_parse_read_lengths_nonexistent_file(self):
        """Test with non-existent file raises error"""

        result = parse_read_lengths(Path("nonexistent.json"))
        assert result is None

    def test_parse_read_lengths_empty_file_warning(self, tmp_path):
        """Test that empty/missing file triggers warning"""
        json_file = tmp_path / "empty.json"

        with patch('lib.io.msg.warning') as mock_warning:
            result = parse_read_lengths(json_file)

            # Verify warning was called with expected message
            mock_warning.assert_called_once_with(
                "Warning: Empty read-lengths parameter given, using all available read lengths."
            )

            assert result is None

    def test_parse_read_lengths_sorted_output(self, tmp_path):
        """Test that output is sorted"""
        json_data = {"sample1": "35,28,30,29"}
        json_file = tmp_path / "read_lengths.json"
        json_file.write_text(json.dumps(json_data))

        result = parse_read_lengths(json_file)

        # Should be sorted
        assert result["sample1"] == ["28", "29", "30", "35"]


class TestParseTotalReads:
    """Tests for parse_total_reads function"""

    def test_parse_total_reads_min_normalization(self, tmp_path):
        """Test parsing with min normalization method"""
        content = "sample1\tchr1\t1000\nsample2\tchr1\t800\nsample1\tchr2\t1500\nsample2\tchr2\t1200\n"
        counts_file = tmp_path / "counts.txt"
        counts_file.write_text(content)

        result = parse_total_reads(counts_file, "min")

        assert result["chr1"] == 800
        assert result["chr2"] == 1200

    def test_parse_total_reads_min_with_single_sample(self, tmp_path):
        """Test with single sample per chromosome"""
        content = "sample1\tchr1\t1000\nsample1\tchr2\t1500\n"
        counts_file = tmp_path / "counts.txt"
        counts_file.write_text(content)

        result = parse_total_reads(counts_file, "min")

        assert result["chr1"] == 1000
        assert result["chr2"] == 1500

    def test_parse_total_reads_non_min_method_returns_none(self, tmp_path):
        """Test that non-min methods return None"""
        counts_file = tmp_path / "counts.txt"
        counts_file.touch()

        result = parse_total_reads(counts_file, "raw")
        assert result is None

        result = parse_total_reads(counts_file, "mil")
        assert result is None

    def test_parse_total_reads_missing_file_min_method(self):
        """Test min method with missing file raises error"""
        counts_file = Path("nonexistent.txt")

        with pytest.raises(FileNotFoundError):
            parse_total_reads(counts_file, "min")

    def test_parse_total_reads_empty_lines_ignored(self, tmp_path):
        """Test that empty lines are properly ignored"""
        content = "sample1\tchr1\t1000\n\n\nsample2\tchr1\t800\n\n"
        counts_file = tmp_path / "counts.txt"
        counts_file.write_text(content)

        result = parse_total_reads(counts_file, "min")

        assert result["chr1"] == 800


class TestParseAlignmentInput:
    """Tests for parse_alignment_input function"""

    def test_parse_alignment_input_both_files(self, tmp_path):
        """Test with both TIS and TTS files"""
        tis_file = tmp_path / "tis.bam"
        tts_file = tmp_path / "tts.bam"
        tis_file.touch()
        tts_file.touch()

        result = parse_alignment_input(tis_file, tts_file)

        assert result == "combined_methods"

    def test_parse_alignment_input_only_tis(self, tmp_path):
        """Test with only TIS file"""
        tis_file = tmp_path / "tis.bam"
        tis_file.touch()

        result = parse_alignment_input(tis_file, None)

        assert result == "TIS"

    def test_parse_alignment_input_only_tts(self, tmp_path):
        """Test with only TTS file"""
        tts_file = tmp_path / "tts.bam"
        tts_file.touch()

        result = parse_alignment_input(None, tts_file)

        assert result == "TTS"

    def test_parse_alignment_input_no_files(self):
        """Test with no files raises error"""
        with pytest.raises(FileNotFoundError, match="Error: Please ensure to either provide a TIS file, a TTS file or both!"):
            parse_alignment_input(None, None)

    def test_parse_alignment_input_tis_not_exists(self, tmp_path):
        """Test with non-existent TIS file raises error"""
        tis_file = tmp_path / "nonexistent.bam"

        with pytest.raises(FileNotFoundError, match="Error: Non-empty alignment file path given for TIS does not exist:"):
            parse_alignment_input(tis_file, None)

    def test_parse_alignment_input_tts_not_exists(self, tmp_path):
        """Test with non-existent TTS file raises error"""
        tts_file = tmp_path / "nonexistent.bam"

        with pytest.raises(FileNotFoundError, match="Error: Non-empty alignment file path given for TTS does not exist:"):
            parse_alignment_input(None, tts_file)

    def test_parse_alignment_input_both_tis_missing(self, tmp_path):
        """Test with TIS missing but TTS exists"""
        tis_file = tmp_path / "nonexistent_tis.bam"
        tts_file = tmp_path / "tts.bam"
        tts_file.touch()

        with pytest.raises(FileNotFoundError, match="Error: Non-empty alignment file path given for TIS does not exist:"):
            parse_alignment_input(tis_file, tts_file)


class TestCheckAlignmentPathInput:
    """Tests for check_alignment_path_input function"""

    def test_check_alignment_path_input_none_returns_none(self):
        """Test with None alignment path returns None"""
        result = check_alignment_path_input(None, None, None)
        assert result is None

    def test_check_alignment_path_input_not_directory(self, tmp_path):
        """Test with non-directory path returns None"""
        file_path = tmp_path / "file.txt"
        file_path.touch()

        result = check_alignment_path_input(file_path, None, None)
        assert result is None

    def test_check_alignment_path_input_with_tis(self, tmp_path):
        """Test finding TIS-related files"""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()

        # Create TIS and RNATIS files
        tis_file = alignment_dir / "TIS-control-rep1.bam"
        rnatis_file = alignment_dir / "RNATIS-control-rep1.bam"
        tis_file.touch()
        rnatis_file.touch()

        tis_input = tmp_path / "TIS-control-rep1.bam"
        tis_input.touch()

        result = check_alignment_path_input(alignment_dir, tis_input, None)

        assert result is not None
        assert len(result) == 2
        assert tis_file in result
        assert rnatis_file in result

    def test_check_alignment_path_input_with_tts(self, tmp_path):
        """Test finding TTS-related files"""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()

        tts_file = alignment_dir / "TTS-control-rep1.bam"
        rnatts_file = alignment_dir / "RNATTS-control-rep1.bam"
        tts_file.touch()
        rnatts_file.touch()

        tts_input = tmp_path / "TTS-control-rep1.bam"
        tts_input.touch()

        result = check_alignment_path_input(alignment_dir, None, tts_input)

        assert result is not None
        assert len(result) == 2
        assert tts_file in result
        assert rnatts_file in result

    def test_check_alignment_path_input_includes_ribo_rna(self, tmp_path):
        """Test that RIBO and RNA files are included"""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()

        tis_file = alignment_dir / "TIS-control-rep1.bam"
        ribo_file = alignment_dir / "RIBO-control-rep1.bam"
        rna_file = alignment_dir / "RNA-control-rep1.bam"

        tis_file.touch()
        ribo_file.touch()
        rna_file.touch()

        tis_input = tmp_path / "TIS-control-rep1.bam"
        tis_input.touch()

        result = check_alignment_path_input(alignment_dir, tis_input, None)

        assert result is not None
        assert len(result) == 3
        assert tis_file in result
        assert ribo_file in result
        assert rna_file in result

    def test_check_alignment_path_input_no_matches(self, tmp_path):
        """Test with no matching files returns None"""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()

        other_file = alignment_dir / "unrelated.bam"
        other_file.touch()

        tis_input = tmp_path / "TIS-control-rep1.bam"
        tis_input.touch()

        result = check_alignment_path_input(alignment_dir, tis_input, None)

        assert result is None

    def test_check_alignment_path_input_filters_sam_bam_only(self, tmp_path):
        """Test that only .sam and .bam files are considered"""
        alignment_dir = tmp_path / "alignments"
        alignment_dir.mkdir()

        tis_file = alignment_dir / "TIS-control-rep1.bam"
        txt_file = alignment_dir / "TIS-control-rep1.txt"
        tis_file.touch()
        txt_file.touch()

        tis_input = tmp_path / "TIS-control-rep1.bam"
        tis_input.touch()

        result = check_alignment_path_input(alignment_dir, tis_input, None)

        assert tis_file in result
        assert txt_file not in result


class TestParseOffsetJson:
    """Tests for parse_offset_json function"""

    def test_parse_offset_json_valid(self, tmp_path):
        """Test parsing valid offset JSON"""
        offset_data = {
            "sample1": {
                "28": 13,
                "29": 12,
                "30": 12,
                "default": 12
            },
            "default": {
                "28": 12,
                "default": 12
            }
        }

        json_file = tmp_path / "offsets.json"
        json_file.write_text(json.dumps(offset_data))

        result = parse_offset_json(json_file)

        assert result == offset_data
        assert result["sample1"]["28"] == 13
        assert result["sample1"]["29"] == 12
        assert result["default"]["default"] == 12

    def test_parse_offset_json_nonexistent(self, tmp_path):
        """Test with non-existent file raises error"""
        json_file = tmp_path / "nonexistent.json"

        with pytest.raises(FileNotFoundError, match="Error: Offset JSON file does not exist!"):
            parse_offset_json(json_file)

    def test_parse_offset_json_nested_structure(self, tmp_path):
        """Test parsing nested structure is preserved"""
        offset_data = {
            "sample1": {"28": 10, "29": 11, "30": 12},
            "sample2": {"28": 15, "default": 14}
        }

        json_file = tmp_path / "offsets.json"
        json_file.write_text(json.dumps(offset_data))

        result = parse_offset_json(json_file)

        assert "sample1" in result
        assert "sample2" in result
        assert result["sample1"]["28"] == 10
        assert result["sample2"]["default"] == 14

    def test_parse_offset_json_empty_file(self, tmp_path):
        """Test with empty file raises error"""
        json_file = tmp_path / "empty.json"
        json_file.write_text("")

        with pytest.raises(ValueError, match="Error: Offset JSON file is empty or invalid!"):
            parse_offset_json(json_file)

class TestWriteGffFile:
    """Tests for write_gff_file function"""

    @pytest.fixture
    def sample_dataframe(self):
        """Create sample DataFrame for GFF output"""
        data = {
            "chromosome": ["chr1", "chr2"],
            "source": ["ORFBounder", "ORFBounder"],
            "type": ["CDS", "CDS"],
            "start": [100, 200],
            "stop": [200, 300],
            "score": [".", "."],
            "strand": ["+", "-"],
            "phase": [".", "."],
            "attribute": ["ID=1;Name=gene1", "ID=2;Name=gene2"]
        }
        return pd.DataFrame(data)

    def test_write_gff_file_creates_file(self, tmp_path, sample_dataframe):
        """Test that GFF file is created"""
        output_path = tmp_path / "output"
        output_file = Path("test.gff")

        write_gff_file(sample_dataframe, output_path, output_file)

        full_path = output_path / output_file
        assert full_path.exists()

    def test_write_gff_file_has_header(self, tmp_path, sample_dataframe):
        """Test that GFF file has proper header"""
        output_path = tmp_path / "output"
        output_file = Path("test.gff")

        write_gff_file(sample_dataframe, output_path, output_file)

        full_path = output_path / output_file
        content = full_path.read_text()

        assert content.startswith("##gff-version 3\n")

    def test_write_gff_file_content(self, tmp_path, sample_dataframe):
        """Test GFF file content is correct"""
        output_path = tmp_path / "output"
        output_file = Path("test.gff")

        write_gff_file(sample_dataframe, output_path, output_file)

        full_path = output_path / output_file
        content = full_path.read_text()

        assert "chr1" in content
        assert "chr2" in content
        assert "ORFBounder" in content
        assert "ID=1;Name=gene1" in content
        assert "ID=2;Name=gene2" in content

    def test_write_gff_file_creates_nested_directories(self, tmp_path, sample_dataframe):
        """Test that nested directories are created"""
        output_path = tmp_path / "output"
        output_file = Path("nested") / "dir" / "test.gff"

        write_gff_file(sample_dataframe, output_path, output_file)

        full_path = output_path / output_file
        assert full_path.exists()
        assert full_path.parent.exists()


class TestWriteCodonIntervalGff:
    """Tests for write_codon_interval_gff function"""

    def test_write_codon_interval_gff_basic(self, tmp_path):
        """Test basic codon interval GFF writing"""
        codon_dict = {
            "chr1:100-103:+": ["ATG", 50],
            "chr1:200-203:-": ["ATG", 30]
        }

        output_path = tmp_path / "output"
        output_basename = "codons.gff"

        write_codon_interval_gff(output_path, output_basename, codon_dict)

        gff_file = output_path / output_basename
        assert gff_file.exists()

    def test_write_codon_interval_gff_empty(self, tmp_path):
        """Test basic codon interval GFF writing"""
        codon_dict = {}

        output_path = tmp_path / "output"
        output_basename = "codons.gff"
        with pytest.raises(ValueError, match="Error: Codon dictionary is empty, cannot create codon interval gff file!"):
            write_codon_interval_gff(output_path, output_basename, codon_dict)


    def test_write_codon_interval_gff_filters_zero_peaks(self, tmp_path):
        """Test that entries with peak_height <= 0 are filtered"""
        codon_dict = {
            "chr1:100-103:+": ["ATG", 50],
            "chr1:200-203:-": ["ATG", 0],
            "chr2:150-153:+": ["ATG", -5]
        }

        output_path = tmp_path / "output"
        output_basename = "codons.gff"

        write_codon_interval_gff(output_path, output_basename, codon_dict)

        gff_file = output_path / output_basename
        content = gff_file.read_text()

        # Should only contain the entry with peak height 50
        assert "chr1" in content
        assert "Peak_height=50" in content
        # The zero and negative peak heights should not be in the file
        lines = content.strip().split('\n')
        data_lines = [l for l in lines if not l.startswith('#')]
        assert len(data_lines) == 1

    def test_write_codon_interval_gff_format(self, tmp_path):
        """Test GFF format is correct"""
        codon_dict = {
            "chr1:100-103:+": ["ATG", 50]
        }

        output_path = tmp_path / "output"
        output_basename = "codons.gff"

        write_codon_interval_gff(output_path, output_basename, codon_dict)

        gff_file = output_path / output_basename
        content = gff_file.read_text()

        assert "##gff-version 3" in content
        assert "chr1" in content
        assert "ORFBounder" in content
        assert "codon_interval" in content
        assert "Start_codon=ATG" in content
        assert "Peak_height=50" in content


class TestExcelWriter:
    """Tests for excel_writer function"""

    @pytest.fixture
    def sample_dataframes(self):
        """Sample DataFrames for Excel output"""
        df1 = pd.DataFrame({
            "Nucleotide_Seq": ["ATGCCC", "ATGAAA"],
            "Amino_Acid_Seq": ["MP", "MK"],
            "Start_codon": ["ATG", "ATG"],
            "Count": [100, 200]
        })

        df2 = pd.DataFrame({
            "Gene": ["gene1", "gene2"],
            "Expression": [50.5, 75.3]
        })

        return {"Sheet1": df1, "Sheet2": df2}

    def test_excel_writer_creates_file(self, tmp_path, sample_dataframes):
        """Test that Excel file is created"""
        output_file = tmp_path / "output.xlsx"

        excel_writer(output_file, sample_dataframes)

        assert output_file.exists()

    def test_excel_writer_multiple_sheets(self, tmp_path, sample_dataframes):
        """Test that multiple sheets are created"""
        output_file = tmp_path / "output.xlsx"

        excel_writer(output_file, sample_dataframes)

        # Read back the Excel file
        xl_file = pd.ExcelFile(output_file)

        assert len(xl_file.sheet_names) == 2
        assert "Sheet1" in xl_file.sheet_names
        assert "Sheet2" in xl_file.sheet_names

    def test_excel_writer_preserves_data(self, tmp_path, sample_dataframes):
        """Test that data is preserved correctly"""
        output_file = tmp_path / "output.xlsx"

        excel_writer(output_file, sample_dataframes)

        # Read back and verify
        df1_read = pd.read_excel(output_file, sheet_name="Sheet1")
        df2_read = pd.read_excel(output_file, sheet_name="Sheet2")

        assert len(df1_read) == 2
        assert len(df2_read) == 2
        assert list(df1_read.columns) == ["Nucleotide_Seq", "Amino_Acid_Seq", "Start_codon", "Count"]
        assert list(df2_read.columns) == ["Gene", "Expression"]


class TestWriteResultsToGff:
    """Tests for write_results_to_gff function"""

    @pytest.fixture
    def sample_results_df(self):
        """Sample results DataFrame"""
        data = {
            "Gene_type": ["Annotated", "Unannotated", "Near_Annotated"],
            "Identifier": ["ID1", "ID2", "ID3"],
            "Chromosome": ["chr1", "chr1", "chr2"],
            "Start": [100, 200, 300],
            "Stop": [200, 300, 400],
            "Strand": ["+", "-", "+"],
            "Locus_tag": ["gene1", "gene2", "gene3"],
            "Codon_count": [33, 50, 25],
            "RPM_TIS": [100.5, 200.3, 150.0],
            "RPM_TTS": [95.2, 195.5, 145.0],
            "RPM_RIBO": [120.0, 220.0, 160.0],
            "Start_codon": ["ATG", "ATG", "GTG"],
            "Stop_codon": ["TAA", "TAG", "TGA"]
        }
        return pd.DataFrame(data)

    def test_write_results_to_gff_creates_main_file(self, tmp_path, sample_results_df):
        """Test that main GFF file is created"""
        output_path = tmp_path / "output" / "gff_per_condition"
        output_basename = "results"

        write_results_to_gff(sample_results_df, output_path, output_basename, split_gff=False)

        main_gff = output_path / f"{output_basename}.gff"
        assert main_gff.exists()

    def test_write_results_to_gff_contains_all_entries(self, tmp_path, sample_results_df):
        """Test that main file contains all entries"""
        output_path = tmp_path / "output" / "gff_per_condition"
        output_basename = "results"

        write_results_to_gff(sample_results_df, output_path, output_basename, split_gff=False)

        main_gff = output_path / f"{output_basename}.gff"
        content = main_gff.read_text()

        lines = [l for l in content.split('\n') if l and not l.startswith('#')]
        assert len(lines) == 3

    def test_write_results_to_gff_split_creates_separate_files(self, tmp_path, sample_results_df):
        """Test that split_gff creates separate files per type"""
        output_path = tmp_path / "output" / "gff_per_condition"
        output_basename = "results"

        write_results_to_gff(sample_results_df, output_path, output_basename, split_gff=True)

        annotated_gff = output_path / f"{output_basename}_annotated.gff"
        unannotated_gff = output_path / f"{output_basename}_unannotated.gff"
        near_annotated_gff = output_path / f"{output_basename}_near_annotated.gff"

        assert annotated_gff.exists()
        assert unannotated_gff.exists()
        assert near_annotated_gff.exists()

    def test_write_results_to_gff_no_split_only_main_file(self, tmp_path, sample_results_df):
        """Test that split_gff=False only creates main file"""
        output_path = tmp_path / "output" / "gff_per_condition"
        output_basename = "results"

        write_results_to_gff(sample_results_df, output_path, output_basename, split_gff=False)

        annotated_gff = output_path / f"{output_basename}_annotated.gff"

        # Split files should not exist
        assert not annotated_gff.exists()


class TestWriteResultsToTable:
    """Tests for write_results_to_table function"""

    @pytest.fixture
    def sample_results_df(self):
        """Sample results DataFrame"""
        data = {
            "Gene_type": ["Annotated", "Unannotated"],
            "Identifier": ["ID1", "ID2"],
            "Locus_tag": ["gene1", "gene2"],
            "Expression": [100.5, 200.3]
        }
        return pd.DataFrame(data)

    def test_write_results_to_table_creates_csv(self, tmp_path, sample_results_df):
        """Test CSV file creation"""
        output_path = tmp_path / "output"
        output_basename = "results"

        write_results_to_table(sample_results_df, output_path, output_basename)

        csv_file = output_path / f"{output_basename}.csv"
        assert csv_file.exists()

    def test_write_results_to_table_creates_xlsx(self, tmp_path, sample_results_df):
        """Test XLSX file creation"""
        output_path = tmp_path / "output"
        output_basename = "results"

        write_results_to_table(sample_results_df, output_path, output_basename)

        xlsx_file = output_path / f"{output_basename}.xlsx"
        assert xlsx_file.exists()

    def test_write_results_to_table_csv_content(self, tmp_path, sample_results_df):
        """Test CSV content is correct"""
        output_path = tmp_path / "output"
        output_basename = "results"

        write_results_to_table(sample_results_df, output_path, output_basename)

        csv_file = output_path / f"{output_basename}.csv"
        content = csv_file.read_text()

        assert "Gene_type" in content
        assert "Identifier" in content
        assert "Annotated" in content
        assert "Unannotated" in content

    def test_write_results_to_table_xlsx_has_cds_sheet(self, tmp_path, sample_results_df):
        """Test XLSX has CDS sheet"""
        output_path = tmp_path / "output"
        output_basename = "results"

        write_results_to_table(sample_results_df, output_path, output_basename)

        xlsx_file = output_path / f"{output_basename}.xlsx"
        xl_file = pd.ExcelFile(xlsx_file)

        assert "CDS" in xl_file.sheet_names