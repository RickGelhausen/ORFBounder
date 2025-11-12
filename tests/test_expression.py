"""
Unit tests for expression calculation functions
"""

from collections import OrderedDict
import re
from unittest.mock import Mock, patch

import pytest
from interlap import InterLap
import numpy as np

from lib.expression import (
    OrderedCounter,
    RNAMAP,
    header_to_dictionary,
    get_te_header,
    calculate_rpkm,
    te_value,
    get_avg,
    te_value_to_dictionary,
    calculate_te,
    init_read_count_dict,
    create_interlap_dict,
    count_reads,
    retrieve_read_counts
)


class TestOrderedCounter:
    """Tests for OrderedCounter class"""

    def test_ordered_counter_basic(self):
        """Test that OrderedCounter maintains order"""
        counter = OrderedCounter()
        counter['a'] = 1
        counter['b'] = 2
        counter['c'] = 3

        keys = list(counter.keys())
        assert keys == ['a', 'b', 'c']

    def test_ordered_counter_counter_functionality(self):
        """Test that Counter functionality works"""
        counter = OrderedCounter(['a', 'b', 'a', 'c', 'b', 'a'])

        assert counter['a'] == 3
        assert counter['b'] == 2
        assert counter['c'] == 1

    def test_ordered_counter_most_common(self):
        """Test most_common method"""
        counter = OrderedCounter(['a', 'b', 'a', 'c', 'b', 'a'])

        most_common = counter.most_common(2)
        assert most_common[0] == ('a', 3)
        assert most_common[1] == ('b', 2)


class TestHeaderToDictionary:
    """Tests for header_to_dictionary function"""

    def test_header_to_dictionary_adds_new_entry(self):
        """Test adding a new header entry"""
        cur_dict = OrderedDict()
        wildcards = ["RNATIS-control-rep1", "TIS-control-rep1"]

        header_to_dictionary("TIS", "control", "rep1", wildcards, cur_dict)

        assert ("TIS", "control") in cur_dict
        assert cur_dict[("TIS", "control")] == ["rep1"]

    def test_header_to_dictionary_appends_to_existing(self):
        """Test appending to existing entry"""
        cur_dict = OrderedDict()
        cur_dict[("TIS", "control")] = ["rep1"]
        wildcards = ["RNATIS-control-rep2", "TIS-control-rep2"]

        header_to_dictionary("TIS", "control", "rep2", wildcards, cur_dict)

        assert cur_dict[("TIS", "control")] == ["rep1", "rep2"]

    def test_header_to_dictionary_no_matching_wildcard(self):
        """Test when RNA wildcard is not present"""
        cur_dict = OrderedDict()
        wildcards = ["TIS-control-rep1"]  # No RNATIS

        header_to_dictionary("TIS", "control", "rep1", wildcards, cur_dict)

        # Should not add to dictionary if RNA wildcard not present
        assert ("TIS", "control") not in cur_dict

    def test_header_to_dictionary_multiple_conditions(self):
        """Test with multiple conditions"""
        cur_dict = OrderedDict()
        wildcards = ["RNATIS-control-rep1", "RNATIS-treated-rep1"]

        header_to_dictionary("TIS", "control", "rep1", wildcards, cur_dict)
        header_to_dictionary("TIS", "treated", "rep1", wildcards, cur_dict)

        assert ("TIS", "control") in cur_dict
        assert ("TIS", "treated") in cur_dict


class TestGetTeHeader:
    """Tests for get_te_header function"""

    def test_get_te_header_single_method(self):
        """Test with single method and replicate"""
        wildcards = ["TIS-control-rep1", "RNATIS-control-rep1"]

        result = get_te_header(wildcards)

        assert result == ["TIS-control-rep1"]

    def test_get_te_header_multiple_replicates(self):
        """Test with multiple replicates"""
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TIS-control-rep2", "RNATIS-control-rep2"
        ]

        result = get_te_header(wildcards)

        assert "TIS-control-rep1" in result
        assert "TIS-control-rep2" in result
        assert len(result) == 2

    def test_get_te_header_multiple_methods(self):
        """Test with multiple methods"""
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TTS-control-rep1", "RNATTS-control-rep1",
            "RIBO-control-rep1", "RNA-control-rep1"
        ]

        result = get_te_header(wildcards)

        assert "TIS-control-rep1" in result
        assert "TTS-control-rep1" in result
        assert "RIBO-control-rep1" in result

    def test_get_te_header_filters_rna(self):
        """Test that RNA wildcards are filtered out"""
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "RNA-control-rep1"
        ]

        result = get_te_header(wildcards)

        # Only TIS should be in result, not RNA
        assert "TIS-control-rep1" in result
        assert "RNA-control-rep1" not in result

    def test_get_te_header_preserves_order(self):
        """Test that order is preserved"""
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TTS-control-rep1", "RNATTS-control-rep1"
        ]

        result = get_te_header(wildcards)

        assert result[0] == "TIS-control-rep1"
        assert result[1] == "TTS-control-rep1"

    def test_get_te_header_no_dash_filtered(self):
        """Test that wildcards without dashes are filtered"""
        wildcards = ["invalid", "TIS-control-rep1", "RNATIS-control-rep1"]

        result = get_te_header(wildcards)

        assert result == ["TIS-control-rep1"]


class TestCalculateRpkm:
    """Tests for calculate_rpkm function"""

    def test_calculate_rpkm_basic(self):
        """Test basic RPKM calculation"""
        result = calculate_rpkm(total_mapped=1000000, read_count=100, read_length=1000)

        expected = (100 * 1000000000) / (1000000 * 1000)
        assert result == expected

    def test_calculate_rpkm_zero_read_length(self):
        """Test with zero read length returns 0"""
        with patch('lib.expression.msg.warning') as mock_warning:
            result = calculate_rpkm(total_mapped=1000000, read_count=100, read_length=0)

            assert result == 0.0
            mock_warning.assert_called_once()
            assert "read_length: 0" in mock_warning.call_args[0][0]

    def test_calculate_rpkm_zero_total_mapped(self):
        """Test with zero total mapped returns 0"""
        with patch('lib.expression.msg.warning') as mock_warning:
            result = calculate_rpkm(total_mapped=0, read_count=100, read_length=1000)

            assert result == 0.0
            mock_warning.assert_called_once()
            assert "total_mapped: 0" in mock_warning.call_args[0][0]

    def test_calculate_rpkm_rounds_to_two_decimals(self):
        """Test that result is rounded to 2 decimal places"""
        result = calculate_rpkm(total_mapped=333333, read_count=123, read_length=456)

        # Should be rounded to 2 decimal places
        assert isinstance(result, float)
        # Check it has at most 2 decimal places
        print(result)
        assert len(str(result).split('.', maxsplit=1)[-1]) <= 2

    def test_calculate_rpkm_small_values(self):
        """Test with very small values"""
        result = calculate_rpkm(total_mapped=1000000000, read_count=1, read_length=1)

        assert result == 1.0

    def test_calculate_rpkm_large_values(self):
        """Test with large values"""
        result = calculate_rpkm(total_mapped=100, read_count=1000, read_length=100)

        expected = (1000 * 1000000000) / (100 * 100)
        assert result == expected

    def test_calculate_rpkm_test(self):
        """Test that result is rounded to 2 decimal places"""
        result = calculate_rpkm(total_mapped=333333, read_count=123, read_length=456)


        assert result == 809.21 # 809.211335527125

class TestTeValue:
    """Tests for te_value function"""

    def test_te_value_basic(self):
        """Test basic TE calculation"""
        result = te_value(ribo_count=100.0, rna_count=50.0)

        assert result == 2.0

    def test_te_value_both_zero_returns_nan(self):
        """Test that both zero returns NaN"""
        result = te_value(ribo_count=0.0, rna_count=0.0)

        assert np.isnan(result)

    def test_te_value_rna_zero_returns_nan(self):
        """Test that RNA zero returns NaN (avoid division by zero)"""
        result = te_value(ribo_count=100.0, rna_count=0.0)

        assert np.isnan(result)

    def test_te_value_ribo_zero_valid(self):
        """Test that RIBO zero with RNA non-zero returns 0"""
        result = te_value(ribo_count=0.0, rna_count=50.0)

        assert result == 0.0

    def test_te_value_fractional(self):
        """Test with fractional values"""
        result = te_value(ribo_count=33.5, rna_count=67.0)

        assert result == 0.5


class TestGetAvg:
    """Tests for get_avg function"""

    def test_get_avg_basic(self):
        """Test basic average calculation"""
        t_eff = [1.0, 2.0, 3.0]

        result = get_avg(t_eff)

        assert len(result) == 4
        assert result[-1] == 2.0  # Average of 1, 2, 3

    def test_get_avg_with_nan_values(self):
        """Test average calculation with NaN values"""
        t_eff = [1.0, np.nan, 3.0]

        result = get_avg(t_eff)

        # Should average only valid values (1.0 and 3.0)
        assert len(result) == 4
        assert result[-1] == 2.0

    def test_get_avg_all_nan(self):
        """Test with all NaN values returns NaN"""
        t_eff = [np.nan, np.nan, np.nan]

        result = get_avg(t_eff)

        assert len(result) == 4
        assert np.isnan(result[-1])

    def test_get_avg_single_value(self):
        """Test with single value"""
        t_eff = [5.0]

        result = get_avg(t_eff)

        assert result[-1] == 5.0

    def test_get_avg_empty_list(self):
        """Test with empty list"""
        t_eff = []

        result = get_avg(t_eff)

        assert len(result) == 1
        assert np.isnan(result[0])


class TestTeValueToDictionary:
    """Tests for te_value_to_dictionary function"""

    def test_te_value_to_dictionary_adds_new_entry(self):
        """Test adding new TE value"""
        read_dict = OrderedDict()
        read_dict[("TIS", "control", "rep1")] = 100.0
        read_dict[("RNATIS", "control", "rep1")] = 50.0

        cur_dict = OrderedDict()

        te_value_to_dictionary("TIS", "control", "rep1", read_dict, cur_dict)

        assert ("TIS", "control") in cur_dict
        assert cur_dict[("TIS", "control")] == [2.0]

    def test_te_value_to_dictionary_appends_to_existing(self):
        """Test appending to existing TE values"""
        read_dict = OrderedDict()
        read_dict[("TIS", "control", "rep1")] = 100.0
        read_dict[("RNATIS", "control", "rep1")] = 50.0
        read_dict[("TIS", "control", "rep2")] = 200.0
        read_dict[("RNATIS", "control", "rep2")] = 100.0

        cur_dict = OrderedDict()
        cur_dict[("TIS", "control")] = [2.0]

        te_value_to_dictionary("TIS", "control", "rep2", read_dict, cur_dict)

        assert len(cur_dict[("TIS", "control")]) == 2
        assert cur_dict[("TIS", "control")][1] == 2.0

    def test_te_value_to_dictionary_missing_rna(self):
        """Test when RNA value is missing"""
        read_dict = OrderedDict()
        read_dict[("TIS", "control", "rep1")] = 100.0
        # Missing RNATIS

        cur_dict = OrderedDict()

        te_value_to_dictionary("TIS", "control", "rep1", read_dict, cur_dict)

        # Should not add to dictionary if RNA is missing
        assert ("TIS", "control") not in cur_dict


class TestCalculateTe:
    """Tests for calculate_te function"""

    def test_calculate_te_single_replicate(self):
        """Test TE calculation with single replicate"""
        read_list = [100.0, 50.0]  # TIS, RNATIS
        wildcards = ["TIS-control-rep1", "RNATIS-control-rep1"]

        result = calculate_te(read_list, wildcards)

        assert len(result) == 1
        assert result[0] == 2.0

    def test_calculate_te_multiple_replicates(self):
        """Test TE calculation with multiple replicates"""
        read_list = [100.0, 50.0, 200.0, 100.0]  # TIS rep1, RNATIS rep1, TIS rep2, RNATIS rep2
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TIS-control-rep2", "RNATIS-control-rep2"
        ]

        result = calculate_te(read_list, wildcards)

        # Should have 3 values: rep1 TE, rep2 TE, and average
        assert len(result) == 3
        assert result[0] == 2.0
        assert result[1] == 2.0
        assert result[2] == 2.0  # Average

    def test_calculate_te_multiple_methods(self):
        """Test with multiple methods"""
        read_list = [
            100.0, 50.0,  # TIS, RNATIS
            150.0, 75.0   # TTS, RNATTS
        ]
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TTS-control-rep1", "RNATTS-control-rep1"
        ]

        result = calculate_te(read_list, wildcards)

        assert len(result) == 2
        assert result[0] == 2.0  # TIS TE
        assert result[1] == 2.0  # TTS TE

    def test_calculate_te_filters_rna_only(self):
        """Test that RNA-only wildcards are filtered"""
        read_list = [100.0, 50.0, 75.0]
        wildcards = ["TIS-control-rep1", "RNATIS-control-rep1", "RNA-control-rep1"]

        result = calculate_te(read_list, wildcards)

        # Should only calculate TE for TIS, not RNA
        assert len(result) == 1
        assert result[0] == 2.0

    def test_calculate_te_with_nan_values(self):
        """Test TE calculation with NaN values"""
        read_list = [100.0, 0.0, 200.0, 100.0]  # Second has 0 RNA -> NaN
        wildcards = [
            "TIS-control-rep1", "RNATIS-control-rep1",
            "TIS-control-rep2", "RNATIS-control-rep2"
        ]

        result = calculate_te(read_list, wildcards)

        assert len(result) == 3
        assert np.isnan(result[0])  # rep1 has NaN
        assert result[1] == 2.0     # rep2 is valid
        assert result[2] == 2.0     # Average of valid only

    def test_calculate_te_duplicate_keys_warns(self):
        """Test that duplicate keys trigger warning"""
        read_list = [100.0, 50.0, 100.0]
        wildcards = ["TIS-control-rep1", "RNATIS-control-rep1", "TIS-control-rep1"]

        with patch('lib.expression.msg.warning') as mock_warning:
            result = calculate_te(read_list, wildcards)

            mock_warning.assert_called_once()
            assert "multiple equal keys" in mock_warning.call_args[0][0]


class TestInitReadCountDict:
    """Tests for init_read_count_dict function"""

    def test_init_read_count_dict_adds_new_intervals(self):
        """Test that new intervals are added"""
        read_count_dict = {}
        result_dict = {
            ("chr1", "+"): {
                (100, 200): {"some": "data"},
                (300, 400): {"some": "data"}
            }
        }

        result = init_read_count_dict(read_count_dict, result_dict)

        assert ("chr1", 100, 200, "+") in result
        assert ("chr1", 300, 400, "+") in result
        assert result[("chr1", 100, 200, "+")] == []

    def test_init_read_count_dict_skips_existing(self):
        """Test that existing intervals are not overwritten"""
        read_count_dict = {
            ("chr1", 100, 200, "+"): [10, 20]
        }
        result_dict = {
            ("chr1", "+"): {
                (100, 200): {"some": "data"}
            }
        }

        result = init_read_count_dict(read_count_dict, result_dict)

        # Existing data should be preserved
        assert result[("chr1", 100, 200, "+")] == [10, 20]

    def test_init_read_count_dict_multiple_chromosomes(self):
        """Test with multiple chromosomes"""
        read_count_dict = {}
        result_dict = {
            ("chr1", "+"): {(100, 200): {}},
            ("chr2", "-"): {(300, 400): {}}
        }

        result = init_read_count_dict(read_count_dict, result_dict)

        assert ("chr1", 100, 200, "+") in result
        assert ("chr2", 300, 400, "-") in result


class TestCreateInterlapDict:
    """Tests for create_interlap_dict function"""

    def test_create_interlap_dict_basic(self, tmp_path):
        """Test basic interlap dict creation"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 1
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_samfile = Mock()
            mock_samfile.__enter__ = Mock(return_value=mock_samfile)
            mock_samfile.__exit__ = Mock(return_value=None)
            mock_samfile.fetch.return_value = [mock_read]
            mock_pysam.return_value = mock_samfile

            interlap_dict, total_mapped = create_interlap_dict(bam_file)

            assert ("chr1", "+") in interlap_dict
            assert "chr1" in total_mapped
            assert total_mapped["chr1"] == 1

    def test_create_interlap_dict_filters_multi_mapped(self, tmp_path):
        """Test that multi-mapped reads are filtered"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        mock_read = Mock()
        mock_read.reference_name = "chr1"
        mock_read.get_tag.return_value = 2
        mock_read.mapping_quality = 20
        mock_read.is_unmapped = False
        mock_read.reference_start = 100
        mock_read.query_length = 30
        mock_read.is_reverse = False

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_samfile = Mock()
            mock_samfile.__enter__ = Mock(return_value=mock_samfile)
            mock_samfile.__exit__ = Mock(return_value=None)
            mock_samfile.fetch.return_value = [mock_read]
            mock_pysam.return_value = mock_samfile

            interlap_dict, total_mapped = create_interlap_dict(bam_file)

            assert len(total_mapped) == 0

    def test_create_interlap_dict_both_strands(self, tmp_path):
        """Test with reads on both strands"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        mock_read_plus = Mock()
        mock_read_plus.reference_name = "chr1"
        mock_read_plus.get_tag.return_value = 1
        mock_read_plus.mapping_quality = 20
        mock_read_plus.is_unmapped = False
        mock_read_plus.reference_start = 100
        mock_read_plus.reference_end = 130
        mock_read_plus.query_length = 30
        mock_read_plus.is_reverse = False

        mock_read_minus = Mock()
        mock_read_minus.reference_name = "chr1"
        mock_read_minus.get_tag.return_value = 1
        mock_read_minus.mapping_quality = 20
        mock_read_minus.is_unmapped = False
        mock_read_minus.reference_start = 200
        mock_read_minus.reference_end = 230
        mock_read_minus.query_length = 28
        mock_read_minus.is_reverse = True

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_samfile = Mock()
            mock_samfile.__enter__ = Mock(return_value=mock_samfile)
            mock_samfile.__exit__ = Mock(return_value=None)
            mock_samfile.fetch.return_value = [mock_read_plus, mock_read_minus]
            mock_pysam.return_value = mock_samfile

            interlap_dict, total_mapped = create_interlap_dict(bam_file)

            assert ("chr1", "+") in interlap_dict
            assert ("chr1", "-") in interlap_dict
            assert total_mapped["chr1"] == 2

    def test_create_interlap_dict_value_error_exits(self, tmp_path):
        """Test that ValueError is caught and exits"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        with patch('pysam.AlignmentFile') as mock_pysam:
            mock_samfile = Mock()
            mock_samfile.__enter__ = Mock(return_value=mock_samfile)
            mock_samfile.__exit__ = Mock(return_value=None)
            mock_samfile.fetch.side_effect = ValueError("Index not found")
            mock_pysam.return_value = mock_samfile

            with pytest.raises(ValueError, match=re.escape("Error: Ensure that all bam files used for readcounting have an appropriate index file (.bam.bai). You can create them using samtools index.")):
                create_interlap_dict(bam_file)


class TestCountReads:
    """Tests for count_reads function"""

    def test_count_reads_basic(self):
        """Test basic read counting"""


        interlap = InterLap()
        interlap.update([(100, 150), (120, 170), (200, 250)])

        read_interlap_dict = {("chr1", "+"): interlap}

        # Query region that overlaps with first two intervals
        count = count_reads("chr1", 110, 160, "+", read_interlap_dict)

        assert count == 2

    def test_count_reads_no_overlap(self):
        """Test region with no overlapping reads"""

        interlap = InterLap()
        interlap.update([(100, 150), (200, 250)])

        read_interlap_dict = {("chr1", "+"): interlap}

        # Query region with no overlap
        count = count_reads("chr1", 160, 190, "+", read_interlap_dict)

        assert count == 0

    def test_count_reads_complete_overlap(self):
        """Test region completely overlapping multiple reads"""

        interlap = InterLap()
        interlap.update([(100, 150), (120, 170), (140, 190)])

        read_interlap_dict = {("chr1", "+"): interlap}

        # Query large region overlapping all
        count = count_reads("chr1", 50, 200, "+", read_interlap_dict)

        assert count == 3

    def test_count_reads_different_strand(self):
        """Test counting reads on different strand"""

        interlap = InterLap()
        interlap.update([(100, 150), (200, 250)])

        read_interlap_dict = {("chr1", "-"): interlap}

        interlap_plus = InterLap()
        interlap_plus.update([(300, 350)])
        read_interlap_dict[("chr1", "+")] = interlap_plus

        # Query on plus strand should find no reads
        count = count_reads("chr1", 100, 220, "+", read_interlap_dict)

        assert count == 0

class TestRetrieveReadCounts:
    """Tests for retrieve_read_counts function"""

    def test_retrieve_read_counts_basic(self, tmp_path):
        """Test basic read count retrieval"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        read_count_dict = {
            ("chr1", 100, 200, "+"): []
        }

        read_lengths = None
        all_reads_rpkm = True

        with patch('lib.expression.IntervalReader') as mock_reader:
            interlap = InterLap()
            interlap.update([(100, 200)])

            mock_instance = Mock()
            mock_instance.output.return_value = (
                {("chr1", "+"): interlap},
                {"chr1": 100}
            )
            mock_reader.return_value = mock_instance

            result_dict, accepted_reads = retrieve_read_counts(
                read_count_dict,
                [bam_file],
                read_lengths,
                all_reads_rpkm
            )

            assert len(result_dict[("chr1", 100, 200, "+")]) == 1
            assert len(accepted_reads) == 1
            assert accepted_reads[0]["chr1"] == 100

    def test_retrieve_read_counts_multiple_bam_files(self, tmp_path):
        """Test with multiple BAM files"""
        bam_file1 = tmp_path / "test1.bam"
        bam_file2 = tmp_path / "test2.bam"
        bam_file1.touch()
        bam_file2.touch()

        read_count_dict = {
            ("chr1", 100, 200, "+"): []
        }

        with patch('lib.expression.IntervalReader') as mock_reader:
            from interlap import InterLap
            interlap = InterLap()
            interlap.update([(100, 200)])

            mock_instance = Mock()
            mock_instance.output.return_value = (
                {("chr1", "+"): interlap},
                {"chr1": 50}
            )
            mock_reader.return_value = mock_instance

            result_dict, accepted_reads = retrieve_read_counts(
                read_count_dict,
                [bam_file1, bam_file2],
                None,
                True
            )

            # Should have counts from both files
            assert len(result_dict[("chr1", 100, 200, "+")]) == 2
            assert len(accepted_reads) == 2

    def test_retrieve_read_counts_multiple_intervals(self, tmp_path):
        """Test with multiple intervals"""
        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        read_count_dict = {
            ("chr1", 100, 200, "+"): [],
            ("chr1", 300, 400, "+"): []
        }

        with patch('lib.expression.IntervalReader') as mock_reader:
            from interlap import InterLap
            interlap = InterLap()
            interlap.update([(100, 200), (300, 400)])

            mock_instance = Mock()
            mock_instance.output.return_value = (
                {("chr1", "+"): interlap},
                {"chr1": 100}
            )
            mock_reader.return_value = mock_instance

            result_dict, accepted_reads = retrieve_read_counts(
                read_count_dict,
                [bam_file],
                None,
                True
            )

            # Both intervals should have counts
            assert len(result_dict[("chr1", 100, 200, "+")]) == 1
            assert len(result_dict[("chr1", 300, 400, "+")]) == 1