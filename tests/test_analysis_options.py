"""Scientific option checks using independently calculable sequences and libraries."""

import copy

import numpy as np
import pysam
import pytest
from Bio.Seq import Seq

from lib.alignment_reader import (
    LIBRARY_MINIMUM_KEY, IntervalReader, PositionReader, validate_normalization_scope,
)
from lib.expression import calculate_rpkm, normalization_read_count, retrieve_read_counts
from lib.genetic_code import resolve_genetic_code
from lib.misc import generate_result_dataframe, get_genome_information


def test_genetic_code_defaults_and_independent_codon_overrides():
    assert resolve_genetic_code() == (["ATG", "GTG", "TTG"], ["TAG", "TAA", "TGA"])
    starts, stops = resolve_genetic_code(4)
    assert "ATA" in starts
    assert set(stops) == {"TAA", "TAG"}
    assert resolve_genetic_code(4, " atg ,GTG,atg", ["tag"]) == (["ATG", "GTG"], ["TAG"])
    assert resolve_genetic_code(4, stop_codons=["TAA"])[0] == starts
    assert resolve_genetic_code(4, start_codons=["ATG"])[1] == stops


@pytest.mark.parametrize("code", [0, 999, -1, 11.0, "11", True, None])
def test_invalid_genetic_codes_fail_with_actionable_error(code):
    with pytest.raises(ValueError, match="genetic_code"):
        resolve_genetic_code(code)


@pytest.mark.parametrize("code", [27, 28, 31])
def test_context_dependent_tables_are_explicitly_unsupported(code):
    with pytest.raises(ValueError, match="context-dependent"):
        resolve_genetic_code(code)


@pytest.mark.parametrize("starts,stops,message", [
    ([], None, "DNA triplets"),
    (["AT"], None, "DNA triplets"),
    (["AUG"], None, "DNA triplets"),
    (["ANN"], None, "DNA triplets"),
    (["AAA"], None, "not start codons"),
    (["TAA"], None, "must not overlap"),
    (None, ["TGA"], "not stop codons"),
    (None, [], "DNA triplets"),
])
def test_custom_codons_must_agree_with_selected_table(starts, stops, message):
    with pytest.raises(ValueError, match=message):
        resolve_genetic_code(4, starts, stops)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("code,sequence,expected", [
    (11, "GTGGTGTAA", "MV*"),
    (11, "TTGTTGTAG", "ML*"),
    (4, "GTGTGATAA", "MW*"),
    (2, "ATTTGAAGA", "MW*"),
])
def test_translation_uses_selected_table_and_initiator_methionine(strand, code, sequence, expected):
    coding_region = sequence if strand == "+" else str(Seq(sequence).reverse_complement())
    genome = "C" * 15 + coding_region + "G" * 15
    nt, protein, upstream, start, stop = get_genome_information(15, 23, strand, genome, genetic_code=code)
    assert nt == sequence
    assert protein == expected
    assert upstream == "C" * 15
    assert (start, stop) == (sequence[:3], sequence[-3:])


def test_internal_stop_filter_uses_selected_genetic_code():
    args = ({("chr", "+"): {(0, 8): (10, np.nan, np.nan)}}, {}, {}, {},
            {"chr": "ATGTGATAA"}, {}, [], [], ("TIS", "TTS", "RIBO"))
    assert generate_result_dataframe(*args, genetic_code=11).empty
    selected = generate_result_dataframe(*args, genetic_code=4)
    assert selected["Amino_Acid_Seq"].tolist() == ["MW*"]


@pytest.fixture
def two_contig_alignment(tmp_path):
    path = tmp_path / "TIS-condition-1.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr", "LN": 200}, {"SN": "plasmid", "LN": 200}]}
    # Three accepted reads on chr (both strands), one on plasmid. A long read
    # is excluded by endpoint length filtering but optionally included in RPKM.
    records = [(0, 0, 10), (0, 0, 10), (0, 16, 10), (1, 0, 10), (1, 0, 20), (1, 256, 10)]
    with pysam.AlignmentFile(path, "w", header=header) as output:
        for number, (reference, flag, length) in enumerate(records):
            read = pysam.AlignedSegment(output.header)
            read.query_name = f"read{number}"
            read.flag = flag
            read.reference_id = reference
            read.reference_start = 20
            read.mapping_quality = 30
            read.cigarstring = f"{length}M"
            read.query_sequence = "A" * length
            output.write(read)
    return path


def position_reader(path):
    return PositionReader(path, {"default": ["10"]}, "fiveprime", {"default": {"default": 0}})


@pytest.mark.parametrize("scope,expected", [
    ("contig", {("chr", "+"): {20: 2000000 / 3}, ("chr", "-"): {29: 1000000 / 3},
                ("plasmid", "+"): {20: 1000000}}),
    ("library", {("chr", "+"): {20: 500000}, ("chr", "-"): {29: 250000},
                 ("plasmid", "+"): {20: 250000}}),
])
def test_rpm_scope_uses_accepted_alignments_on_both_strands(two_contig_alignment, scope, expected):
    reader = position_reader(two_contig_alignment)
    reader.normalize_read_counts("mil", {}, normalization_scope=scope)
    actual, accepted = reader.output()
    assert accepted == {"chr": 3, "plasmid": 1}
    assert actual.keys() == expected.keys()
    for key in expected:
        assert actual[key] == pytest.approx(expected[key])


def test_default_contig_scaling_and_raw_library_scope_preserve_behavior(two_contig_alignment):
    default = position_reader(two_contig_alignment)
    before = copy.deepcopy(default.output()[0])
    default.normalize_read_counts("raw", {}, normalization_scope="library")
    assert default.output()[0] == before
    default.normalize_read_counts("mil", {})
    assert default.output()[0][("plasmid", "+")][20] == 1000000


def test_invalid_scope_rejects_before_mutation(two_contig_alignment):
    reader = position_reader(two_contig_alignment)
    before = copy.deepcopy(reader.output())
    with pytest.raises(ValueError, match="normalization_scope"):
        reader.normalize_read_counts("raw", {}, normalization_scope="unknown")
    assert reader.output() == before


def test_library_min_uses_one_whole_library_target(two_contig_alignment):
    reader = position_reader(two_contig_alignment)
    reader.normalize_read_counts(
        "min", {LIBRARY_MINIMUM_KEY: 20}, normalization_scope="library",
    )
    positions, accepted = reader.output()
    assert accepted == {"chr": 3, "plasmid": 1}
    assert positions[("chr", "+")][20] == pytest.approx(10)
    assert positions[("chr", "-")][29] == pytest.approx(5)
    assert positions[("plasmid", "+")][20] == pytest.approx(5)


def test_library_min_requires_library_target_without_mutation(two_contig_alignment):
    reader = position_reader(two_contig_alignment)
    before = copy.deepcopy(reader.output())
    with pytest.raises(ValueError, match="minimum read count found for the library"):
        reader.normalize_read_counts("min", {"chr": 2, "plasmid": 3}, normalization_scope="library")
    assert reader.output() == before


def test_empty_library_scaling_has_no_division_by_zero(tmp_path):
    path = tmp_path / "empty.sam"
    path.write_text("@HD\tVN:1.6\n@SQ\tSN:chr\tLN:100\n", encoding="utf-8")
    reader = position_reader(path)
    reader.normalize_read_counts("mil", {}, normalization_scope="library")
    assert reader.output() == ({}, {})


@pytest.mark.parametrize("all_reads,total", [(False, 4), (True, 5)])
def test_library_rpkm_denominator_obeys_expression_length_selection(two_contig_alignment, all_reads, total):
    intervals, accepted = IntervalReader(two_contig_alignment, {"default": ["10"]}, all_reads).output()
    assert normalization_read_count(accepted, "chr", "library") == total
    assert normalization_read_count(accepted, "chr", "contig") == 3
    reads, accepted_list = retrieve_read_counts({("chr", 20, 28, "+"): []},
                                                [two_contig_alignment], {"default": ["10"]}, all_reads)
    df = generate_result_dataframe(
        {("chr", "+"): {(20, 28): (2, np.nan, np.nan)}}, {}, {}, {},
        {"chr": "C" * 20 + "ATGAAATAA"}, reads, accepted_list, ["TIS-condition-1"],
        ("TIS", "TTS", "RIBO"), normalization_scope="library")
    assert df["TIS-condition-1_rpkm"].tolist() == [round(2e9 / (total * 9), 2)]
    assert df["TIS_peak_height"].tolist() == [2]


def test_expression_missing_contigs_and_empty_library_are_defined():
    assert normalization_read_count({"chr": 3}, "missing", "contig") == 0
    assert normalization_read_count({"chr": 3}, "missing", "library") == 3
    assert calculate_rpkm(normalization_read_count({}, "missing", "library"), 0, 9) == 0.0
    with pytest.raises(ValueError, match="normalization_scope"):
        validate_normalization_scope("invalid")


def test_rpkm_and_te_share_each_samples_selected_scope():
    args = ({("chr", "+"): {(0, 8): (2, np.nan, np.nan)}}, {}, {}, {},
            {"chr": "ATGAAATAA"}, {("chr", 0, 8, "+"): [2, 4]},
            [{"chr": 10, "plasmid": 90}, {"chr": 80, "plasmid": 120}],
            ["TIS-condition-1", "RNATIS-condition-1"], ("TIS", "TTS", "RIBO"))
    library = generate_result_dataframe(*args, normalization_scope="library")
    contig = generate_result_dataframe(*args, normalization_scope="contig")
    assert library["TIS-condition-1_TE"].tolist() == [1.0]
    assert contig["TIS-condition-1_TE"].iloc[0] == pytest.approx(4.0)


def test_te_uses_unrounded_rates_when_displayed_rpkm_rounds_to_zero():
    frame = generate_result_dataframe(
        {("chr", "+"): {(0, 8): (2, np.nan, np.nan)}}, {}, {}, {},
        {"chr": "ATGAAATAA"}, {("chr", 0, 8, "+"): [1, 1]},
        [{"chr": 10**12}, {"chr": 2 * 10**12}],
        ["TIS-condition-1", "RNATIS-condition-1"],
        ("TIS", "TTS", "RIBO"),
    )

    assert frame["TIS-condition-1_rpkm"].tolist() == [0.0]
    assert frame["RNATIS-condition-1_rpkm"].tolist() == [0.0]
    assert frame["TIS-condition-1_TE"].tolist() == [2.0]
