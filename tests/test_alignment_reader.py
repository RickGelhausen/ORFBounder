"""
Unit tests for PositionReader and IntervalReader classes
"""

from pathlib import Path

import re
from unittest.mock import Mock, patch
import pytest

from lib.alignment_reader import PositionReader, IntervalReader
from lib import messaging as msg

class TestPositionReader:
    """Tests for PositionReader class"""

    @pytest.fixture
    def mock_alignment_file(self, tmp_path):
        """Create a mock alignment file path"""
        bam_file = tmp_path / "TIS-WT-1.bam"
        bam_file.touch()
        return bam_file

    @pytest.fixture
    def basic_read_length_dict(self):
        """Basic read length dictionary"""
        return {"TIS-WT-1": ["28", "29", "30"], "default": ["28", "29", "30"]}

    @pytest.fixture
    def basic_offset_dict(self):
        """Basic offset dictionary"""
        return {"TIS-WT-1": {28: 12, 29: 12, 30: 14, "default": 12}, "default": {28: 12, 29: 12, 30: 12, "default": 12}}

    @pytest.fixture
    def mock_read(self):
        """Create a mock read object"""
        read = Mock()
        read.reference_name = "chr1"
        read.get_tag.return_value = 1  # NH tag
        read.mapping_quality = 20
        read.is_unmapped = False
        read.reference_start = 100
        read.query_length = 30
        read.is_reverse = False
        return read

    def test_init_with_wildcard_specific_read_lengths(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test initialization with wildcard-specific read lengths"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            assert reader.wildcard == "TIS-WT-1"
            assert reader.read_lengths == ["28", "29", "30"]
            assert reader.mapping_mode == "fiveprime"

    def test_init_with_default_read_lengths(self, tmp_path, basic_offset_dict):
        """Test initialization falls back to default read lengths"""
        bam_file = tmp_path / "unknown_sample.bam"
        bam_file.touch()

        read_length_dict = {"default": ["25", "26", "27"]}

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = PositionReader(
                bam_file,
                read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            assert reader.read_lengths == ["25", "26", "27"]

    def test_init_without_read_lengths(self, mock_alignment_file, basic_offset_dict):
        """Test initialization with no read length filtering"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = PositionReader(
                mock_alignment_file,
                {},
                "fiveprime",
                basic_offset_dict
            )

            assert reader.read_lengths is None

    def test_init_missing_offset_raises_error(self, mock_alignment_file, basic_read_length_dict):
        """Test initialization with missing offset raises ValueError"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            with pytest.raises(ValueError, match=msg.error("Error: Offsets and default value missing for wildcard: test")):
                PositionReader(
                    mock_alignment_file,
                    basic_read_length_dict,
                    "fiveprime",
                    {}
                )

    def test_read_alignment_file_filters_multi_mapped(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test that multi-mapped reads (NH > 1) are filtered out"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 2  # NH tag > 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            assert len(reader.reads_position_dict) == 0

    def test_read_alignment_file_filters_unmapped(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test that unmapped reads are filtered out"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            assert len(reader.reads_position_dict) == 0

    def test_read_alignment_file_filters_low_quality(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test that low quality mapped reads are filtered out"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = -1
        mock_read.is_unmapped = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            assert len(reader.reads_position_dict) == 0

    def test_fiveprime_mapping_plus_strand(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test fiveprime mapping mode on plus strand"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            # For plus strand, fiveprime: start - offset = 100 - 14 = 86
            assert ("chr1", "+") in reader.reads_position_dict
            assert 86 in reader.reads_position_dict[("chr1", "+")]
            assert reader.reads_position_dict[("chr1", "+")][86] == 1

    def test_fiveprime_mapping_minus_strand(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test fiveprime mapping mode on minus strand"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            # For minus strand, fiveprime: stop + offset = 129 + 14 = 143
            assert ("chr1", "-") in reader.reads_position_dict
            assert 143 in reader.reads_position_dict[("chr1", "-")]
            assert reader.reads_position_dict[("chr1", "-")][143] == 1

    def test_threeprime_mapping_plus_strand(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test threeprime mapping mode on plus strand"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "threeprime",
                basic_offset_dict
            )

            # For plus strand, threeprime: stop - offset = 129 - 14 = 115
            assert ("chr1", "+") in reader.reads_position_dict
            assert 115 in reader.reads_position_dict[("chr1", "+")]
            assert reader.reads_position_dict[("chr1", "+")][115] == 1

    def test_threeprime_mapping_minus_strand(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test threeprime mapping mode on minus strand"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "threeprime",
                basic_offset_dict
            )

            # For minus strand, threeprime: start + offset = 100 + 14 = 114
            assert ("chr1", "-") in reader.reads_position_dict
            assert 114 in reader.reads_position_dict[("chr1", "-")]
            assert reader.reads_position_dict[("chr1", "-")][114] == 1

    def test_centered_mapping_mode(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test centered mapping mode"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "centered",
                basic_offset_dict
            )

            # Centered mode clips 11 from each end
            # center_start = 100 + 11 = 111
            # center_stop = 129 - 11 = 118
            # center_length = 118 - 111 + 1 = 8
            # Each position gets 1/8 count
            assert ("chr1", "+") in reader.reads_position_dict
            # Check that positions are weighted
            for pos in reader.reads_position_dict[("chr1", "+")].values():
                assert pos == pytest.approx(1/8)

    def test_global_mapping_mode(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test global mapping mode"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "global",
                basic_offset_dict
            )

            # Global mode should cover all positions from start to stop
            assert ("chr1", "+") in reader.reads_position_dict
            # Should have 30 positions (with offset applied)
            assert len(reader.reads_position_dict[("chr1", "+")]) == 30

    def test_negative_positions_filtered(self, mock_alignment_file, basic_read_length_dict):
        """Test that negative positions are filtered out"""
        offset_dict = {"test": {30: 150}, "default": {30: 150}}

        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                offset_dict
            )

            # start - offset = 100 - 150 = -50 (should be filtered)
            # The position dict should be empty or not contain negative positions
            if ("chr1", "+") in reader.reads_position_dict:
                for pos in reader.reads_position_dict[("chr1", "+")].keys():
                    assert pos >= 0

    def test_multiple_reads_same_position(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test that multiple reads at the same position are counted"""
        mock_read1 = Mock()
        mock_read1.reference_name = "chr1"
        mock_read1.get_tag.return_value = 1
        mock_read1.mapping_quality = 20
        mock_read1.is_unmapped = False
        mock_read1.reference_start = 100
        mock_read1.query_length = 30
        mock_read1.is_reverse = False

        mock_read2 = Mock()
        mock_read2.reference_name = "chr1"
        mock_read2.get_tag.return_value = 1
        mock_read2.mapping_quality = 20
        mock_read2.is_unmapped = False
        mock_read2.reference_start = 100
        mock_read2.query_length = 30
        mock_read2.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read1, mock_read2]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            # Both reads should contribute to the same position
            assert reader.reads_position_dict[("chr1", "+")][86] == 2

    def test_normalize_read_counts_raw(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read):
        """Test normalization with raw mode (no normalization)"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            original_count = reader.reads_position_dict[("chr1", "+")][86]
            reader.normalize_read_counts("raw", {})

            # Should remain unchanged
            assert reader.reads_position_dict[("chr1", "+")][86] == original_count

    def test_normalize_read_counts_mil(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read):
        """Test normalization with per million reads"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            reader.normalize_read_counts("mil", {})

            # Should be normalized to per million
            # 1 read / 1 total read * 1000000 = 1000000
            assert reader.reads_position_dict[("chr1", "+")][86] == pytest.approx(1000000.0)

    def test_normalize_read_counts_min(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read):
        """Test normalization with minimum read count"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            min_read_count_dict = {"chr1": 500}
            reader.normalize_read_counts("min", min_read_count_dict)

            # Should be normalized: 1 * (500 / 1) = 500
            assert reader.reads_position_dict[("chr1", "+")][86] == pytest.approx(500.0)

    def test_normalize_read_counts_invalid_method(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read):
        """Test normalization with invalid method raises error"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )
            with pytest.raises(ValueError, match=(re.escape("Error: Given normalization method is not supported: invalid. Supported methods: (raw, min, mil)"))):
                reader.normalize_read_counts("invalid", {})

    def test_to_wig_output(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read, tmp_path):
        """Test wiggle file output"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            output_dir = tmp_path / "wig_output"
            reader.to_wig(output_dir)

            # Check that file was created
            expected_file = output_dir / "TIS-WT-1_chr1_forward.wig"
            assert expected_file.exists()

            # Check file contents
            content = expected_file.read_text()
            assert "track type=wiggle_0" in content
            assert "variableStep chrom=chr1" in content

    def test_output_method(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict, mock_read):
        """Test output method returns correct data structures"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = PositionReader(
                mock_alignment_file,
                basic_read_length_dict,
                "fiveprime",
                basic_offset_dict
            )

            position_dict, read_count_dict = reader.output()

            assert isinstance(position_dict, dict)
            assert isinstance(read_count_dict, dict)
            assert ("chr1", "+") in position_dict
            assert "chr1" in read_count_dict
            assert read_count_dict["chr1"] == 1

    def test_value_error_handling(self, mock_alignment_file, basic_read_length_dict, basic_offset_dict):
        """Test that ValueError is caught and exits with proper message"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.side_effect = ValueError("Index not found")

            with pytest.raises(ValueError, match=(re.escape("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index."))):
                    reader = PositionReader(
                        mock_alignment_file,
                        basic_read_length_dict,
                        "fiveprime",
                        basic_offset_dict
                    )


class TestIntervalReader:
    """Tests for IntervalReader class"""

    @pytest.fixture
    def mock_alignment_file(self, tmp_path):
        """Create a mock alignment file path"""
        bam_file = tmp_path / "TIS-WT-1.bam"
        bam_file.touch()
        return bam_file

    @pytest.fixture
    def basic_read_length_dict(self):
        """Basic read length dictionary"""
        return {"TIS-WT-1": ["28", "29", "30"], "default": ["28", "29", "30"]}

    @pytest.fixture
    def mock_read(self):
        """Create a mock read object"""
        read = Mock()
        read.reference_name = "chr1"
        read.get_tag.return_value = 1
        read.mapping_quality = 20
        read.is_unmapped = False
        read.reference_start = 100
        read.query_length = 30
        read.is_reverse = False
        return read

    def test_init_with_specific_read_lengths(self, mock_alignment_file, basic_read_length_dict):
        """Test initialization with specific read lengths"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=False
            )

            assert reader.wildcard == "TIS-WT-1"
            assert reader.read_lengths == ["28", "29", "30"]
            assert reader.rpkm_all_reads is False

    def test_init_with_rpkm_all_reads(self, mock_alignment_file, basic_read_length_dict):
        """Test initialization with rpkm_all_reads flag"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert reader.rpkm_all_reads is True

    def test_init_without_read_lengths_forces_rpkm_all(self, mock_alignment_file):
        """Test that None read_lengths forces rpkm_all_reads to True"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = IntervalReader(
                mock_alignment_file,
                None,
                rpkm_all_reads=False
            )

            assert reader.read_lengths is None
            assert reader.rpkm_all_reads is True

    def test_read_alignment_file_creates_interlap(self, mock_alignment_file, basic_read_length_dict, mock_read):
        """Test that interlap objects are created correctly"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert ("chr1", "+") in reader.reads_interlap_dict
            # Interlap object should be created
            assert hasattr(reader.reads_interlap_dict[("chr1", "+")], 'find')

    def test_filters_multi_mapped_reads(self, mock_alignment_file, basic_read_length_dict):
        """Test that multi-mapped reads are filtered"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 2  # NH > 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert len(reader.reads_interlap_dict) == 0

    def test_filters_unmapped_reads(self, mock_alignment_file, basic_read_length_dict):
        """Test that unmapped reads are filtered"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert len(reader.reads_interlap_dict) == 0

    def test_filters_low_quality_reads(self, mock_alignment_file, basic_read_length_dict):
        """Test that low quality reads are filtered"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = -1
        mock_read.is_unmapped = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert len(reader.reads_interlap_dict) == 0

    def test_filters_by_read_length_when_not_rpkm_all(self, mock_alignment_file, basic_read_length_dict):
        """Test that reads are filtered by length when rpkm_all_reads is False"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 50  # Not in the allowed list
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=False
            )

            # Read should be filtered out
            assert len(reader.reads_interlap_dict) == 0

    def test_accepts_all_lengths_when_rpkm_all(self, mock_alignment_file, basic_read_length_dict):
        """Test that all read lengths are accepted when rpkm_all_reads is True"""
        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 50  # Not in the allowed list
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            # Read should be accepted
            assert len(reader.reads_interlap_dict) > 0

    def test_counts_reads_correctly(self, mock_alignment_file, basic_read_length_dict, mock_read):
        """Test that accepted reads are counted correctly"""
        mock_read2 = Mock()
        mock_read2.reference_name = "chr1"
        mock_read2.get_tag.return_value = 1
        mock_read2.mapping_quality = 20
        mock_read2.is_unmapped = False
        mock_read2.reference_start = 200
        mock_read2.query_length = 30
        mock_read2.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read, mock_read2]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert reader.no_accepted_reads_dict["chr1"] == 2

    def test_handles_multiple_chromosomes(self, mock_alignment_file, basic_read_length_dict):
        """Test handling of reads from multiple chromosomes"""
        mock_read1 = Mock()
        mock_read1.reference_name = "chr1"
        mock_read1.get_tag.return_value = 1
        mock_read1.mapping_quality = 20
        mock_read1.is_unmapped = False
        mock_read1.reference_start = 100
        mock_read1.query_length = 30
        mock_read1.is_reverse = False

        mock_read2 = Mock()
        mock_read2.reference_name = "chr2"
        mock_read2.get_tag.return_value = 1
        mock_read2.mapping_quality = 20
        mock_read2.is_unmapped = False
        mock_read2.reference_start = 100
        mock_read2.query_length = 30
        mock_read2.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read1, mock_read2]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert ("chr1", "+") in reader.reads_interlap_dict
            assert ("chr2", "+") in reader.reads_interlap_dict
            assert reader.no_accepted_reads_dict["chr1"] == 1
            assert reader.no_accepted_reads_dict["chr2"] == 1

    def test_handles_both_strands(self, mock_alignment_file, basic_read_length_dict):
        """Test handling of reads from both strands"""
        mock_read_plus = Mock()
        mock_read_plus.reference_name = "chr1"
        mock_read_plus.get_tag.return_value = 1
        mock_read_plus.mapping_quality = 20
        mock_read_plus.is_unmapped = False
        mock_read_plus.reference_start = 100
        mock_read_plus.query_length = 30
        mock_read_plus.is_reverse = False

        mock_read_minus = Mock()
        mock_read_minus.reference_name = "chr1"
        mock_read_minus.get_tag.return_value = 1
        mock_read_minus.mapping_quality = 20
        mock_read_minus.is_unmapped = False
        mock_read_minus.reference_start = 100
        mock_read_minus.query_length = 30
        mock_read_minus.is_reverse = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read_plus, mock_read_minus]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            assert ("chr1", "+") in reader.reads_interlap_dict
            assert ("chr1", "-") in reader.reads_interlap_dict

    def test_output_method(self, mock_alignment_file, basic_read_length_dict, mock_read):
        """Test output method returns correct data structures"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = [mock_read]

            reader = IntervalReader(
                mock_alignment_file,
                basic_read_length_dict,
                rpkm_all_reads=True
            )

            interlap_dict, read_count_dict = reader.output()

            assert isinstance(interlap_dict, dict)
            assert isinstance(read_count_dict, dict)
            assert ("chr1", "+") in interlap_dict
            assert "chr1" in read_count_dict

    def test_value_error_handling(self, mock_alignment_file, basic_read_length_dict):
        """Test that ValueError is caught and exits with proper message"""
        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.side_effect = ValueError("Index not found")

            with pytest.raises(ValueError, match=(re.escape("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index."))):
                reader = IntervalReader(
                    mock_alignment_file,
                    basic_read_length_dict,
                    rpkm_all_reads=True
                )

    def test_default_read_length_fallback(self, tmp_path, basic_read_length_dict):
        """Test fallback to default read lengths for unknown wildcard"""
        bam_file = tmp_path / "unknown_sample.bam"
        bam_file.touch()

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_pysam.return_value.fetch.return_value = []

            reader = IntervalReader(
                bam_file,
                basic_read_length_dict,
                rpkm_all_reads=False
            )

            assert reader.read_lengths == ["28", "29", "30"]