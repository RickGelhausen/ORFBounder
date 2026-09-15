"""Regression tests exercising real alignments and biological coordinates."""

from dataclasses import FrozenInstanceError

import numpy as np
import pysam
import pytest
from Bio.Seq import Seq

from lib.alignment_reader import AlignmentPolicy, IntervalReader, PositionReader
from lib import expression, misc, predictions


def write_alignment(tmp_path, records, suffix=".sam"):
    """Write (name, flag, start, CIGAR, length, NH[, MAPQ]) records without an index."""
    path = tmp_path / f"TIS-control-1{suffix}"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 200}, {"SN": "chr2", "LN": 200}]}
    with pysam.AlignmentFile(path, "wb" if suffix == ".bam" else "w", header=header) as output:
        for record in records:
            name, flag, start, cigar, length, nh, *optional = record
            read = pysam.AlignedSegment(output.header)
            read.query_name = name
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = optional[0] if optional else 20
            read.cigarstring = cigar
            if length is not None:
                read.query_sequence = "A" * length
            if nh is not None:
                read.set_tag("NH", nh)
            output.write(read)
    return path


@pytest.mark.parametrize("arguments,message", [
    ({"min_mapq": -1}, "min_mapq"),
    ({"min_mapq": 256}, "min_mapq"),
    ({"min_mapq": 1.5}, "min_mapq"),
    ({"min_mapq": True}, "min_mapq"),
    ({"include_duplicates": 1}, "include_duplicates"),
    ({"multimappers": "fractional"}, "multimappers"),
])
def test_alignment_policy_validates_settings(arguments, message):
    with pytest.raises(ValueError, match=message):
        AlignmentPolicy(**arguments)


def test_alignment_policy_is_immutable_and_preserves_legacy_defaults():
    policy = AlignmentPolicy()
    assert (policy.min_mapq, policy.include_duplicates, policy.multimappers) == (0, True, "exclude")
    with pytest.raises(FrozenInstanceError):
        policy.min_mapq = 10


def test_alignment_policy_filters_records_and_reports_diagnostics(tmp_path):
    path = write_alignment(tmp_path, [
        ("low_mapq", 0, 10, "10M", 10, 1, 5),
        ("unique", 0, 30, "10M", 10, 1, 30),
        ("duplicate", 1024, 50, "10M", 10, 1, 30),
        ("multimapped", 0, 70, "10M", 10, 2, 30),
        ("without_nh", 0, 90, "10M", 10, None, 30),
        ("paired", 1, 110, "10M", 10, 1, 30),
    ])
    offsets = {"default": {"default": 0}}

    legacy = PositionReader(path, None, "threeprime", offsets)
    assert set(legacy.output()[0][("chr1", "+")]) == {19, 39, 59, 99, 119}
    assert legacy.alignment_diagnostics == {
        "policy": {"min_mapq": 0, "include_duplicates": True, "multimappers": "exclude"},
        "count_unit": "alignment_record",
        "total_records": 6,
        "accepted_records": 5,
        "excluded_records": 1,
        "records_without_nh": 1,
        "paired_records": 1,
        "duplicate_records": 1,
        "multimapped_records": 1,
        "exclusions": {
            "unmapped": 0,
            "secondary": 0,
            "supplementary": 0,
            "qc_fail": 0,
            "below_min_mapq": 0,
            "duplicate": 0,
            "multimapper": 1,
        },
    }

    policy = AlignmentPolicy(min_mapq=20, include_duplicates=False, multimappers="include")
    positions = PositionReader(path, None, "threeprime", offsets, policy)
    assert set(positions.output()[0][("chr1", "+")]) == {39, 79, 99, 119}
    assert positions.output()[1] == {"chr1": 4}
    assert positions.alignment_diagnostics["accepted_records"] == 4
    assert positions.alignment_diagnostics["excluded_records"] == 2
    assert positions.alignment_diagnostics["exclusions"]["below_min_mapq"] == 1
    assert positions.alignment_diagnostics["exclusions"]["duplicate"] == 1
    assert positions.alignment_diagnostics["exclusions"]["multimapper"] == 0

    intervals = IntervalReader(path, None, True, policy)
    assert intervals.output()[1] == {"chr1": 4}
    assert expression.count_reads("chr1", 0, 199, "+", intervals.output()[0]) == 4
    assert intervals.alignment_diagnostics == positions.alignment_diagnostics


@pytest.mark.parametrize("nh", [0, -1, "2"])
@pytest.mark.parametrize("reader", ["position", "interval"])
def test_malformed_nh_tag_fails_before_counting(tmp_path, nh, reader):
    path = write_alignment(tmp_path, [("bad_nh", 0, 30, "10M", 10, nh)])
    with pytest.raises(ValueError, match="invalid NH tag"):
        if reader == "position":
            PositionReader(path, None, "threeprime", {"default": {"default": 0}})
        else:
            IntervalReader(path, None, True)


@pytest.mark.parametrize("suffix", [".sam", ".bam"])
def test_streaming_optional_nh_and_primary_alignment_filter(tmp_path, suffix):
    path = write_alignment(tmp_path, [
        ("primary", 0, 20, "10M", 10, None),
        ("secondary", 256, 40, "10M", 10, None),
        ("supplementary", 2048, 60, "10M", 10, 1),
        ("failed_qc", 512, 80, "10M", 10, 1),
        ("multimapped", 0, 100, "10M", 10, 2),
    ], suffix)
    position_reader = PositionReader(path, None, "threeprime", {"default": {"default": 0}})
    positions, counts = position_reader.output()
    assert positions == {("chr1", "+"): {29: 1}}
    assert counts == {"chr1": 1}
    assert position_reader.alignment_diagnostics["total_records"] == 5
    assert position_reader.alignment_diagnostics["accepted_records"] == 1
    assert position_reader.alignment_diagnostics["excluded_records"] == 4
    assert position_reader.alignment_diagnostics["records_without_nh"] == 2
    assert position_reader.alignment_diagnostics["exclusions"] == {
        "unmapped": 0,
        "secondary": 1,
        "supplementary": 1,
        "qc_fail": 1,
        "below_min_mapq": 0,
        "duplicate": 0,
        "multimapper": 1,
    }
    interval_reader = IntervalReader(path, None, True)
    intervals, counts = interval_reader.output()
    assert expression.count_reads("chr1", 0, 199, "+", intervals) == 1
    assert counts == {"chr1": 1}
    assert interval_reader.alignment_diagnostics == position_reader.alignment_diagnostics


def test_read_length_json_offsets_use_cigar_aware_endpoints(tmp_path):
    path = write_alignment(tmp_path, [
        ("plus", 0, 40, "5S10M2D10M3S", 28, None),
        ("minus", 16, 80, "5S10M2D10M3S", 28, None),
    ])
    offsets = {"default": {"28": 3, "default": 12}}
    threeprime, _ = PositionReader(path, None, "threeprime", offsets).output()
    fiveprime, _ = PositionReader(path, None, "fiveprime", offsets).output()
    assert threeprime == {("chr1", "+"): {58: 1}, ("chr1", "-"): {83: 1}}
    assert fiveprime == {("chr1", "+"): {37: 1}, ("chr1", "-"): {104: 1}}


def test_missing_sam_sequence_uses_cigar_query_length_for_offsets_and_filtering(tmp_path):
    path = write_alignment(tmp_path, [
        ("sequence_omitted", 0, 40, "5S10M2D10M3S", None, None),
    ])
    lengths = {"default": ["28"]}

    positions, position_counts = PositionReader(
        path, lengths, "threeprime", {"default": {"28": 3}},
    ).output()
    intervals, interval_counts = IntervalReader(path, lengths, False).output()

    assert positions == {("chr1", "+"): {58: 1}}
    assert position_counts == {"chr1": 1}
    assert expression.count_reads("chr1", 40, 61, "+", intervals) == 1
    assert interval_counts == {"chr1": 1}


@pytest.mark.parametrize("reader", ["position", "interval"])
def test_mapped_alignment_without_cigar_is_rejected(tmp_path, reader):
    path = write_alignment(tmp_path, [
        ("missing_cigar", 0, 40, None, 10, None),
    ], suffix=".bam")

    with pytest.raises(
        ValueError,
        match="Alignment 'missing_cigar' is marked as mapped but has no CIGAR",
    ):
        if reader == "position":
            PositionReader(path, None, "threeprime", {"default": {"default": 0}})
        else:
            IntervalReader(path, None, True)


def test_cigar_insertions_do_not_extend_genomic_end(tmp_path):
    path = write_alignment(tmp_path, [("insertion", 0, 40, "10M2I10M", 22, None)])
    positions, _ = PositionReader(path, None, "threeprime", {"default": {"default": 0}}).output()
    assert positions == {("chr1", "+"): {59: 1}}


def test_cigar_gaps_have_no_coverage_and_each_alignment_counts_once(tmp_path):
    path = write_alignment(tmp_path, [
        ("first", 0, 20, "5M10N5M", 10, None),
        ("second", 0, 20, "5M10N5M", 10, None),
    ])
    positions, _ = PositionReader(path, None, "global", {"default": {"default": 0}}).output()
    assert set(positions[("chr1", "+")]) == set(range(20, 25)) | set(range(35, 40))
    intervals, _ = expression.create_interlap_dict(path)
    assert expression.count_reads("chr1", 25, 34, "+", intervals) == 0
    assert expression.count_reads("chr1", 24, 35, "+", intervals) == 2
    assert expression.count_reads("chr1", 0, 199, "-", intervals) == 0
    assert expression.count_reads("chr2", 0, 199, "+", intervals) == 0


def test_centered_offsets_are_strand_consistent_and_short_reads_only_lack_coverage(tmp_path):
    path = write_alignment(tmp_path, [
        ("plus", 0, 50, "25M", 25, None),
        ("minus", 16, 50, "25M", 25, None),
        ("too_short", 0, 100, "22M", 22, None),
    ])
    positions, counts = PositionReader(path, None, "centered", {"default": {"default": 2}}).output()
    assert set(positions[("chr1", "+")]) == {59, 60, 61}
    assert set(positions[("chr1", "-")]) == {63, 64, 65}
    assert sum(positions[("chr1", "+")].values()) == pytest.approx(1)
    assert sum(positions[("chr1", "-")].values()) == pytest.approx(1)
    # Selected primary records remain part of the normalization denominator;
    # 22-nt centered reads simply have no positions after 11+11 trimming.
    assert counts == {"chr1": 3}


def test_offsets_outside_reference_do_not_produce_invalid_wig_positions(tmp_path):
    path = write_alignment(tmp_path, [("minus", 16, 180, "20M", 20, None)])
    reader = PositionReader(path, None, "fiveprime", {"default": {"default": 12}})
    assert reader.output()[0] == {}


def test_missing_offset_reports_the_configuration_error(tmp_path):
    path = write_alignment(tmp_path, [("read", 0, 20, "10M", 10, None)])
    with pytest.raises(ValueError, match="Offsets and default value missing for .*readlength"):
        PositionReader(path, None, "threeprime", {"default": {"28": 12}})


def test_gff_density_uses_same_zero_based_coordinates_as_alignment(tmp_path):
    annotation = tmp_path / "genes.gff"
    annotation.write_text("chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds;locus_tag=gene\n")
    intervals, genes = misc.annotation_interlap(annotation)
    misc.calculate_density({("chr1", "+"): {0: 2, 8: 3, 9: 100}}, intervals, genes)
    assert genes["gene"] == ["chr1", 1, 9, "+", 5]


def test_gff_without_locus_tags_keeps_distinct_genes_and_ignores_embedded_fasta(tmp_path):
    annotation = tmp_path / "genes.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "chr1\t.\tregion\t1\t200\t.\t+\t.\t.\n"
        "chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=first\n"
        "chr1\t.\tCDS\t20\t28\t.\t+\t0\tID=second\n"
        "##FASTA\n>chr1\nATGCCCTAA\n"
    )
    intervals, genes = misc.annotation_interlap(annotation)
    misc.calculate_density({("chr1", "+"): {0: 2, 19: 3}}, intervals, genes)
    assert genes["first"][4] == 2
    assert genes["second"][4] == 3


def test_gff_without_cds_records_is_a_valid_empty_annotation(tmp_path):
    annotation = tmp_path / "genes.gff"
    annotation.write_text("##gff-version 3\n")
    assert misc.annotation_interlap(annotation) == ({}, {})


def test_gff_attribute_values_do_not_collide_with_attribute_names(tmp_path):
    annotation = tmp_path / "genes.gff"
    annotation.write_text("chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=name;Name=two%20words\n")
    parsed = misc.generate_annotation_dict(annotation)
    assert parsed["chr1:1-9:+"][0] == "name"
    assert parsed["chr1:1-9:+"][2] == "two words"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("method", ["TIS", "TTS"])
def test_terminal_orfs_are_retained_on_both_strands(strand, method):
    sequence = "ATGCCCTAA"
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    search_codons, match_codons = (["ATG"], ["TAA"]) if method == "TIS" else (["TAA"], ["ATG"])
    intervals, codons = misc.create_codon_interlaps("chr1", sequence, search_codons)
    position = (0 if method == "TIS" else 6) if strand == "+" else (8 if method == "TIS" else 2)
    predictions.screen_positions_for_tss({("chr1", strand): {position: 10}}, intervals, codons, 0, "max")
    result = predictions.detect_potential_orfs(codons, sequence, search_codons, match_codons, method, {}, "furthest_inframe")
    assert set(result[("chr1", strand)]) == {(0, 8)}
    assert misc.get_genome_information(0, 8, strand, sequence)[0] == "ATGCCCTAA"


def test_furthest_tts_search_handles_reference_boundaries():
    assert predictions.search_longest_reverse(0, "TAACCC", ["TAA"], ["ATG"]) is None
    assert predictions.search_longest_forward(0, "TAACCCATGATGCCC", ["TAA"], ["ATG"]) == 9
    assert predictions.search_longest_reverse(9, "ATGATGCCCTAA", ["TAA"], ["ATG"]) == 0


def test_multi_replicate_te_and_missing_contig_expression_form_valid_table():
    wildcards = ["RIBO-control-1", "RIBO-control-2", "RNA-control-1", "RNA-control-2"]
    frame = misc.generate_result_dataframe(
        {("chr1", "+"): {(0, 8): (10, np.nan, np.nan)}}, {}, {}, {},
        {"chr1": "ATGCCCTAA"}, {("chr1", 0, 8, "+"): [10, 20, 5, 0]},
        [{"chr1": 100}, {"chr1": 100}, {"chr1": 100}, {}], wildcards, ("TIS", "TTS", "RIBO"),
    )
    assert len(frame) == 1
    assert frame.loc[0, "RIBO-control-1_TE"] == pytest.approx(2)
    assert np.isnan(frame.loc[0, "RIBO-control-2_TE"])
    assert frame.loc[0, "RNA-control-2_rpkm"] == 0


def test_upstream_window_and_reverse_start_site_distances():
    assert misc.get_genome_information(3, 11, "+", "CCCATGCCCTAA")[2] == "CCC"
    assert misc.calculate_utr_distance(99, 269, "gene", {"gene": ("chr1", 100, 300, "-", 0)}) == (30, 170)


def test_te_ignores_unpaired_arbitrary_names_without_losing_named_pairs():
    names = ["my-sample", "my-sample-name", "TIS", "TIS--1", "TIS-control-1", "RNATIS-control-1"]
    assert expression.get_te_header(names) == ["TIS-control-1"]
    assert expression.calculate_te([100, 200, 300, 400, 10, 5], names) == [2]
    duplicate_names = ["TIS-control-1", "RNATIS-control-1", "TIS-control-1"]
    assert expression.get_te_header(duplicate_names) == ["TIS-control-1"]
    assert expression.calculate_te([10, 5, 10], duplicate_names) == [2]


def test_te_rejects_misaligned_counts_and_names():
    with pytest.raises(ValueError, match="same length"):
        expression.calculate_te([10], ["TIS-control-1", "RNATIS-control-1"])


def test_min_normalization_requires_every_covered_contig_before_changing_counts():
    reader = PositionReader.__new__(PositionReader)
    reader.reads_position_dict = {("chr1", "+"): {0: 2}, ("chr2", "+"): {0: 3}}
    reader.no_accepted_reads_dict = {"chr1": 2, "chr2": 3}
    with pytest.raises(ValueError, match="No minimum read count found for chr2"):
        reader.normalize_read_counts("min", {"chr1": 1})
    assert reader.reads_position_dict == {("chr1", "+"): {0: 2}, ("chr2", "+"): {0: 3}}
