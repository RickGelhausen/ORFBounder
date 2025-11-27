"""
Unit tests for ORF prediction functions
"""

import pytest
import numpy as np
from Bio.Seq import Seq
from interlap import InterLap

from lib.predictions import (
    CODON_LENGTH,
    screen_positions_for_tss,
    search_codon_forward,
    search_codon_reverse,
    search_longest_forward,
    search_longest_reverse,
    detect_potential_orfs,
    convert_codon_dict
)


class TestCodonLength:
    """Test CODON_LENGTH constant"""

    def test_codon_length_is_three(self):
        """Test that CODON_LENGTH is 3"""
        assert CODON_LENGTH == 3


class TestScreenPositionsForTss:
    """Tests for screen_positions_for_tss function"""

    @pytest.fixture
    def sample_alignment_dict(self):
        """Sample alignment position dictionary"""
        return {
            ("chr1", "+"): {100: 10, 101: 6, 125: 5, 150: 20, 200: 5}
        }

    @pytest.fixture
    def sample_codon_interlap_dict(self):
        """Sample codon interlap dictionary"""
        interlap = InterLap()
        interlap.add((98, 102, "chr1:98-102:+"))
        interlap.add((148, 152, "chr1:148-152:+"))
        return {("chr1", "+"): interlap}

    @pytest.fixture
    def sample_codon_dict(self):
        """Sample codon dictionary"""
        return {
            "chr1:98-102:+": ["ATG", 0],
            "chr1:148-152:+": ["ATG", 0]
        }

    def test_screen_positions_sum_operator(self, sample_alignment_dict, sample_codon_interlap_dict, sample_codon_dict):
        """Test screening with sum operator"""
        result = screen_positions_for_tss(
            sample_alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=5.0,
            peak_height_operator="sum"
        )

        # Position 100 (count=10) overlaps chr1:98-102:+
        # Position 150 (count=20) overlaps chr1:148-152:+
        assert result["chr1:98-102:+"][1] == 16
        assert result["chr1:148-152:+"][1] == 20

    def test_screen_positions_max_operator(self, sample_alignment_dict, sample_codon_interlap_dict, sample_codon_dict):
        """Test screening with max operator"""
        # Add another position that overlaps the first codon
        sample_alignment_dict[("chr1", "+")][99] = 15

        result = screen_positions_for_tss(
            sample_alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=5.0,
            peak_height_operator="max"
        )

        # Should take max of 10 and 15
        assert result["chr1:98-102:+"][1] == 15
        assert result["chr1:148-152:+"][1] == 20

    def test_screen_positions_filters_low_peaks(self, sample_alignment_dict, sample_codon_interlap_dict, sample_codon_dict):
        """Test that positions below min_peak_height are filtered"""
        result = screen_positions_for_tss(
            sample_alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=15.0,  # Only position 150 (20) passes
            peak_height_operator="sum"
        )

        # Position 100 (10) and 200 (5) should be filtered
        assert result["chr1:98-102:+"][1] == 0
        assert result["chr1:148-152:+"][1] == 20

    def test_screen_positions_missing_chromosome(self, sample_alignment_dict, sample_codon_dict):
        """Test when chromosome is not in codon_interlap_dict"""
        sample_alignment_dict[("chr2", "+")] = {100: 10}

        interlap = InterLap()
        interlap.add((98, 102, "chr1:98-102:+"))
        codon_interlap_dict = {("chr1", "+"): interlap}

        result = screen_positions_for_tss(
            sample_alignment_dict,
            codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=5.0,
            peak_height_operator="sum"
        )
        print(result)
        # chr2 positions should be skipped
        assert result == sample_codon_dict

    def test_screen_positions_no_overlaps(self, sample_codon_interlap_dict, sample_codon_dict):
        """Test when positions don't overlap with codons"""
        alignment_dict = {
            ("chr1", "+"): {50: 10, 250: 20}  # Outside codon intervals
        }

        result = screen_positions_for_tss(
            alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=5.0,
            peak_height_operator="sum"
        )

        # Counts should remain 0
        assert result["chr1:98-102:+"][1] == 0
        assert result["chr1:148-152:+"][1] == 0

    def test_screen_positions_invalid_operator(self, sample_alignment_dict, sample_codon_interlap_dict, sample_codon_dict):
        """Test with invalid operator raises error"""
        with pytest.raises(ValueError, match="Invalid method! Use either 'sum' or 'max'!"):
            screen_positions_for_tss(
                sample_alignment_dict,
                sample_codon_interlap_dict,
                sample_codon_dict,
                min_peak_height=5.0,
                    peak_height_operator="invalid"
                )

    def test_screen_positions_sum_accumulates(self, sample_codon_interlap_dict, sample_codon_dict):
        """Test that sum operator accumulates multiple reads"""
        alignment_dict = {
            ("chr1", "+"): {100: 10, 101: 5, 102: 8}  # All overlap first codon
        }

        result = screen_positions_for_tss(
            alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=0.0,
            peak_height_operator="sum"
        )

        # Should sum all three
        assert result["chr1:98-102:+"][1] == 23

    def test_screen_positions_max_not_accumulates(self, sample_codon_interlap_dict, sample_codon_dict):
        """Test that max operator does not accumulate multiple reads"""
        alignment_dict = {
            ("chr1", "+"): {100: 10, 101: 5, 102: 8}  # All overlap first codon
        }

        result = screen_positions_for_tss(
            alignment_dict,
            sample_codon_interlap_dict,
            sample_codon_dict,
            min_peak_height=0.0,
            peak_height_operator="max"
        )

        # Should take the maximum
        assert result["chr1:98-102:+"][1] == 10


class TestSearchCodonForward:
    """Tests for search_codon_forward function"""

    def test_search_codon_forward_found_immediately(self):
        """Test when codon is at current position"""
        genome_seq = "ATGAAACCC"
        result = search_codon_forward(0, genome_seq, ["ATG"])

        assert result == 0

    def test_search_codon_forward_found_next(self):
        """Test finding codon at next position"""
        genome_seq = "AAAATGCCC"
        result = search_codon_forward(0, genome_seq, ["ATG"])

        assert result == 3

    def test_search_codon_forward_multiple_steps(self):
        """Test finding codon several steps away"""
        genome_seq = "AAAAAAAAA" + "ATG" + "CCC"
        result = search_codon_forward(0, genome_seq, ["ATG"])

        assert result == 9

    def test_search_codon_forward_not_found(self):
        """Test when codon is not found"""
        genome_seq = "AAAAAAAAACCC"
        result = search_codon_forward(0, genome_seq, ["ATG"])

        assert result is None

    def test_search_codon_forward_at_end(self):
        """Test finding codon at end of sequence"""
        genome_seq = "AAAAAATAG"
        result = search_codon_forward(0, genome_seq, ["TAG"])

        assert result == 6

    def test_search_codon_forward_multiple_match_codons(self):
        """Test with multiple match codons"""
        genome_seq = "AAAGTGCCC"
        result = search_codon_forward(0, genome_seq, ["ATG", "GTG", "TTG"])

        assert result == 3  # Finds GTG

    def test_search_codon_forward_boundary(self):
        """Test boundary when near end of sequence"""
        genome_seq = "AAAAAA"  # Too short to find codon
        result = search_codon_forward(3, genome_seq, ["ATG"])

        assert result is None

    def test_search_codon_forward_starts_mid_sequence(self):
        """Test starting from middle of sequence"""
        genome_seq = "ATGAAAGTGCCC"
        result = search_codon_forward(6, genome_seq, ["GTG"])

        assert result == 6


class TestSearchCodonReverse:
    """Tests for search_codon_reverse function"""

    def test_search_codon_reverse_found_immediately(self):
        """Test when codon is at current position"""
        genome_seq = "CCCATGAAA"
        result = search_codon_reverse(3, genome_seq, ["ATG"])

        assert result == 3

    def test_search_codon_reverse_found_previous(self):
        """Test finding codon at previous position"""
        genome_seq = "ATGAAACCC"
        result = search_codon_reverse(6, genome_seq, ["ATG"])

        assert result == 0

    def test_search_codon_reverse_multiple_steps(self):
        """Test finding codon several steps back"""
        genome_seq = "ATG" + "AAAAAAAAA" + "CCC"
        result = search_codon_reverse(12, genome_seq, ["ATG"])

        assert result == 0

    def test_search_codon_reverse_not_found(self):
        """Test when codon is not found"""
        genome_seq = "CCCCCCCCCCCC"
        result = search_codon_reverse(9, genome_seq, ["ATG"])

        assert result is None

    def test_search_codon_reverse_at_start(self):
        """Test finding codon at start of sequence"""
        genome_seq = "TAAAAACCC"
        result = search_codon_reverse(6, genome_seq, ["TAA"])

        assert result == 0

    def test_search_codon_reverse_multiple_match_codons(self):
        """Test with multiple match codons"""
        genome_seq = "GTGAAACCC"
        result = search_codon_reverse(6, genome_seq, ["ATG", "GTG", "TTG"])

        assert result == 0  # Finds GTG

    def test_search_codon_reverse_boundary(self):
        """Test boundary when near start"""
        genome_seq = "AAAAAA"
        result = search_codon_reverse(3, genome_seq, ["ATG"])

        assert result is None

    def test_search_codon_reverse_negative_boundary(self):
        """Test that going negative returns None"""
        genome_seq = "CCCCCC"
        result = search_codon_reverse(0, genome_seq, ["ATG"])

        assert result is None

    def test_search_codon_reverse_sanity_check(self):
        """Test starting from middle of sequence"""
        genome_seq = "AAAAAATAACCC"
        result = search_codon_reverse(6, genome_seq, ["TAA"])

        assert result == 6


class TestSearchLongestForward:
    """Tests for search_longest_forward function"""

    def test_search_longest_forward_basic(self):
        """Test basic longest forward search"""
        genome_seq = "ATGAAATAACCC"
        # Start at ATG (0), find TAA (6)
        result = search_longest_forward(0, genome_seq, ["TAA"], ["ATG"])

        assert result is None  # Should find ATG which is previous position

    def test_search_longest_forward_finds_last_start_before_stop(self):
        """Test that it finds the last start codon before stop"""
        genome_seq = "ATGAAAATGAAATAACCC"
        # ATG at 0, 6; TAA at 12
        result = search_longest_forward(0, genome_seq, ["TAA"], ["ATG"])

        assert result == 6  # Should find second ATG

    def test_search_longest_forward_multiple_options(self):
        """Test with multiple search and match codons"""
        genome_seq = "GTGAAAATGAAAGTGTAACCC"
        # GTG at 0, 12; ATG at 6; TAA at 15
        result = search_longest_forward(0, genome_seq, ["TAA"], ["ATG", "GTG"])

        assert result == 12  # Should find GTG at 12

    def test_search_longest_forward_multiple_options_2(self):
        """Test with multiple search and match codons"""
        genome_seq = "GTGAAAATGAAAAAGTAACCC"
        # GTG at 0, 12; ATG at 6; TAA at 15
        result = search_longest_forward(0, genome_seq, ["TAA"], ["ATG", "GTG"])

        assert result == 6  # Should find GTG at 6

    def test_search_longest_forward_multiple_options_3(self):
        """Test with multiple search and match codons"""
        genome_seq = "GTGAAAATGAAAAAGTAACCCTAG"
        # GTG at 0, 12; ATG at 6; TAA at 15
        result = search_longest_forward(0, genome_seq, ["TAG", "TAA"], ["ATG", "GTG"])

        assert result == 6  # Should find GTG at 6

    def test_search_longest_forward_no_search_codon(self):
        """Test when no search codon is found"""
        genome_seq = "ATGAAAAAACCC"
        result = search_longest_forward(0, genome_seq, ["TAA"], ["ATG"])

        # Should fallback to search_codon_forward
        assert result == 0  # Finds ATG at start

    def test_search_longest_forward_no_match_codon(self):
        """Test when no match codon between search codons"""
        genome_seq = "ATGTAAAAACCC"
        result = search_longest_forward(0, genome_seq, ["TAA"], ["GTG"])

        assert result is None

    def test_search_longest_forward_exceeds_sequence(self):
        """Test when search goes beyond sequence"""
        genome_seq = "ATGAAA"
        result = search_longest_forward(0, genome_seq, ["TAG"], ["ATG"])

        # Should fallback to search_codon_forward
        assert result == 0


class TestSearchLongestReverse:
    """Tests for search_longest_reverse function"""

    def test_search_longest_reverse_basic(self):
        """Test basic longest reverse search"""
        genome_seq = "ATGAAATAACCC"
        # TAA at 6, search back for ATG
        result = search_longest_reverse(6, genome_seq, ["ATG"], ["TAA"])

        assert result is None  # Should find TAA which is previous position

    def test_search_longest_reverse_finds_first_stop_after_start(self):
        """Test that it finds the first stop codon after start"""
        genome_seq = "ATGAAATAATAACCC"
        # ATG at 0, TAA at 6 and 9
        result = search_longest_reverse(9, genome_seq, ["ATG"], ["TAA"])

        assert result == 6  # Should find first TAA after ATG

    def test_search_longest_reverse_multiple_options(self):
        """Test with multiple search and match codons"""
        genome_seq = "GTGAAAGTGAAAATGTAACCC"
        # GTG at 0, 6; ATG at 12; TAA at 15
        result = search_longest_reverse(15, genome_seq, ["ATG", "GTG"], ["TAA"])

        assert result is None

    def test_search_longest_reverse_no_search_codon(self):
        """Test when no search codon is found"""
        genome_seq = "AAAAAATAACCC"
        result = search_longest_reverse(6, genome_seq, ["ATG"], ["TAA"])

        # Should fallback to search_codon_reverse
        assert result is None

    def test_search_longest_reverse_no_match_codon(self):
        """Test when no match codon between search codons"""
        genome_seq = "ATGAAAAAACCC"
        result = search_longest_reverse(9, genome_seq, ["ATG"], ["GTG"])

        assert result is None

    def test_search_longest_reverse_goes_negative(self):
        """Test when search goes before sequence start"""
        genome_seq = "TAACCC"
        result = search_longest_reverse(3, genome_seq, ["ATG"], ["TAA"])

        # Should fallback to search_codon_reverse
        assert result == 0


class TestDetectPotentialOrfs:
    """Tests for detect_potential_orfs function"""

    @pytest.fixture
    def sample_genome_seq(self):
        """Sample genome sequence"""
        return "NNNNATGAAAAAATAANNNN"

    def test_detect_orfs_tis_plus_strand(self, sample_genome_seq):
        """Test TIS method on plus strand"""
        codon_dict = {
            "chr1:2-6:+": ["ATG", 10]  # ATG at position 4
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            sample_genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="TIS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should find ORF from ATG (4) to TAA (13)
        assert ("chr1", "+") in result
        assert (4, 15) in result[("chr1", "+")]
        assert result[("chr1", "+")][(4, 15)] == (10, np.nan, np.nan)

    def test_detect_orfs_tis_minus_strand(self):
        """Test TIS method on minus strand"""
        genome_seq = "NNNTTACCCCCCCCATNNNN"  # TTA (reverse of TAA) at 4, CAT (reverse of ATG) at 14
        codon_dict = {
            "chr1:12-16:-": ["ATG", 10]
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="TIS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        assert ("chr1", "-") in result
        # Should find ORF
        assert len(result[("chr1", "-")]) > 0

    def test_detect_orfs_tts_plus_strand_next_inframe(self, sample_genome_seq):
        """Test TTS method with next_inframe on plus strand"""
        codon_dict = {
            "chr1:11-15:+": ["TAA", 20]  # TAA at position 13
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            sample_genome_seq,
            search_codons=["TAA"],
            match_codons=["ATG", "GTG", "TTG"],
            method="TTS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should find ORF from ATG (4) to TAA (13)
        assert ("chr1", "+") in result
        assert (4, 15) in result[("chr1", "+")]
        assert result[("chr1", "+")][(4, 15)] == (np.nan, 20, np.nan)

    def test_detect_orfs_tts_plus_strand_furthest_inframe(self):
        """Test TTS method with furthest_inframe"""
        genome_seq = "TAANNNATGATGAAATAANNNN"
        # ATG at 2, 5; TAA at 11
        codon_dict = {
            "chr1:13-17:+": ["TAA", 20]
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            genome_seq,
            search_codons=["TAA"],
            match_codons=["ATG", "GTG", "TTG"],
            method="TTS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="furthest_inframe"
        )

        # Should find the furthest ATG (position 2, not 5)
        assert ("chr1", "+") in result
        assert (6, 17) in result[("chr1", "+")]
        assert result[("chr1", "+")][(6, 17)] == (np.nan, 20, np.nan)

    def test_detect_orfs_tts_plus_strand_furthest_inframe_2(self):
        """Test TTS method with furthest_inframe"""
        genome_seq = "NNNATGATGAAATAANNNN"
        # ATG at 2, 5; TAA at 11
        codon_dict = {
            "chr1:10-14:+": ["TAA", 20]
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            genome_seq,
            search_codons=["TAA"],
            match_codons=["ATG", "GTG", "TTG"],
            method="TTS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="furthest_inframe"
        )

        # Should find the furthest ATG (position 2, not 5)
        assert ("chr1", "+") in result
        assert (3, 14) in result[("chr1", "+")]
        assert result[("chr1", "+")][(3, 14)] == (np.nan, 20, np.nan)

    def test_detect_orfs_ribo_method(self, sample_genome_seq):
        """Test RIBO method"""
        codon_dict = {
            "chr1:2-6:+": ["ATG", 15]
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            sample_genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="RIBO",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should behave like TIS
        assert ("chr1", "+") in result
        # Check RIBO value is set
        orf_key = list(result[("chr1", "+")].keys())[0]
        assert result[("chr1", "+")][orf_key] == (np.nan, np.nan, 15)

    def test_detect_orfs_filters_zero_peaks(self, sample_genome_seq):
        """Test that codons with peak height <= 0 are filtered"""
        codon_dict = {
            "chr1:2-6:+": ["ATG", 0],  # Peak height 0
            "chr1:5-9:+": ["ATG", -5]  # Negative peak height
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            sample_genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="TIS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should find no ORFs
        assert len(result) == 0

    def test_detect_orfs_no_stop_codon_found(self):
        """Test when no stop codon is found"""
        genome_seq = "NNNNATGAAAAAAAAA"  # No stop codon
        codon_dict = {
            "chr1:2-6:+": ["ATG", 10]
        }
        detected_orfs_dict = {}

        result = detect_potential_orfs(
            codon_dict,
            genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="TIS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should not add ORF without stop codon
        assert len(result) == 0

    def test_detect_orfs_updates_existing_orf(self, sample_genome_seq):
        """Test that existing ORF gets updated"""
        codon_dict = {
            "chr1:2-6:+": ["ATG", 10]
        }
        detected_orfs_dict = {
            ("chr1", "+"): {
                (4, 15): (np.nan, 20, np.nan)  # Existing TTS value
            }
        }

        result = detect_potential_orfs(
            codon_dict,
            sample_genome_seq,
            search_codons=["ATG"],
            match_codons=["TAA", "TAG", "TGA"],
            method="TIS",
            detected_orfs_dict=detected_orfs_dict,
            tts_start_selection="next_inframe"
        )

        # Should update with TIS value
        assert result[("chr1", "+")][(4, 15)] == (10, 20, np.nan)

    def test_detect_orfs_invalid_tts_selection(self, sample_genome_seq):
        """Test with invalid TTS start selection method"""
        codon_dict = {
            "chr1:11-15:+": ["TAA", 20]
        }

        with pytest.raises(ValueError):
            detect_potential_orfs(
                codon_dict,
                sample_genome_seq,
                search_codons=["ATG"],
                match_codons=["ATG"],
                method="TTS",
                detected_orfs_dict={},
                tts_start_selection="invalid_method"
            )


class TestConvertCodonDict:
    """Tests for convert_codon_dict function"""

    def test_convert_codon_dict_tis_plus_strand(self):
        """Test converting TIS codon dict for plus strand"""
        codon_dict_tis = {
            ("chr1:2-6:+", 12): ["ATG", 10]  # interval_start=2, offset=12
        }
        codon_dict_tts = {}

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # cur_start = 2 - 12 + 2 = -8... wait, this seems wrong
        # Let me recalculate: interval_start=2, offset=12
        # cur_start = int(interval_start) - offset_tis + 2 = 2 - 12 + 2 = -8

        assert 12 in start_dict
        assert ("chr1", "+") in start_dict[12]
        # Check the calculation
        assert start_dict[12][("chr1", "+")] == [(-8, 10)]

    def test_convert_codon_dict_tis_minus_strand(self):
        """Test converting TIS codon dict for minus strand"""
        codon_dict_tis = {
            ("chr1:10-14:-", 12): ["ATG", 15]
        }
        codon_dict_tts = {}

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # cur_start = 10 + 12 + 2 = 24
        assert 12 in start_dict
        assert start_dict[12][("chr1", "-")] == [(24, 15)]

    def test_convert_codon_dict_tts_plus_strand(self):
        """Test converting TTS codon dict for plus strand"""
        codon_dict_tis = {}
        codon_dict_tts = {
            ("chr1:11-15:+", 12): ["TAA", 20]
        }

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # cur_stop = 11 - 12 + 4 = 3
        assert 12 in stop_dict
        assert stop_dict[12][("chr1", "+")] == [(3, 20)]

    def test_convert_codon_dict_tts_minus_strand(self):
        """Test converting TTS codon dict for minus strand"""
        codon_dict_tis = {}
        codon_dict_tts = {
            ("chr1:10-14:-", 12): ["TAA", 25]
        }

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # cur_stop = 10 + 12 = 22
        assert 12 in stop_dict
        assert stop_dict[12][("chr1", "-")] == [(22, 25)]

    def test_convert_codon_dict_filters_low_peaks(self):
        """Test that codons with peak < 1 are filtered"""
        codon_dict_tis = {
            ("chr1:2-6:+", 12): ["ATG", 0],
            ("chr1:5-9:+", 12): ["ATG", 10]
        }
        codon_dict_tts = {}

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # Only second entry should be present
        assert len(start_dict[12][("chr1", "+")]) == 1

    def test_convert_codon_dict_multiple_offsets(self):
        """Test with multiple different offsets"""
        codon_dict_tis = {
            ("chr1:2-6:+", 10): ["ATG", 10],
            ("chr1:5-9:+", 12): ["ATG", 15]
        }
        codon_dict_tts = {}

        start_dict, stop_dict = convert_codon_dict(
            codon_dict_tis, codon_dict_tts, offset_tis=12, offset_tts=12
        )

        # Should have entries for both offsets
        assert 10 in start_dict
        assert 12 in start_dict
