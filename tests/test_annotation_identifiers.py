"""Regression coverage for annotation identity and coordinate identifiers."""

import numpy as np
import pysam
import pytest

from lib import io, misc, predictions
from lib.coordinates import format_coordinate_id, parse_coordinate_id


def _feature_lines(path):
    return [line for line in path.read_text().splitlines() if line and not line.startswith("#")]


def test_coordinate_identifier_round_trip_uses_right_hand_delimiters():
    identifier = format_coordinate_id("plasmid:alpha-beta", -2, 2, "+")

    assert identifier == "plasmid:alpha-beta:-2-2:+"
    assert parse_coordinate_id(identifier) == ("plasmid:alpha-beta", -2, 2, "+")


def test_coordinate_identifier_round_trip_accepts_a_negative_stop():
    identifier = format_coordinate_id("plasmid:alpha-beta", -5, -1, "-")

    assert identifier == "plasmid:alpha-beta:-5--1:-"
    assert parse_coordinate_id(identifier) == ("plasmid:alpha-beta", -5, -1, "-")


def test_delimited_url_encoded_seqid_flows_through_prediction_and_gff(tmp_path):
    seqid = "plasmid:alpha-beta,copy"
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "plasmid%3Aalpha-beta%2Ccopy\t.\tCDS\t1\t9\t.\t+\t0\t"
        "ID=cds%3Aone;locus_tag=gene-one\n"
    )
    genome = {seqid: "ATGCCCTAA"}

    io.validate_reference_inputs(genome, annotation, [])
    annotation_interlaps, genes = misc.annotation_interlap(annotation)
    codon_interlaps, codons = misc.create_codon_interlaps(seqid, genome[seqid], ["ATG"])
    predictions.screen_positions_for_tss(
        {(seqid, "+"): {0: 10}}, codon_interlaps, codons, 0, "max")
    detected = predictions.detect_potential_orfs(
        codons, genome[seqid], ["ATG"], ["TAA"], "TIS", {}, "next_inframe")

    assert (seqid, "+") in annotation_interlaps
    assert detected == {(seqid, "+"): {(0, 8): (10, np.nan, np.nan)}}

    result = misc.generate_result_dataframe(
        detected, genes, {}, {}, genome, {}, [], [], ("TIS", "TTS", "RIBO"))
    assert result.loc[0, "Identifier"] == f"{seqid}:1-9:+"
    assert result.loc[0, "Locus_tag"] == "gene-one"
    assert result.loc[0, "Type"] == "Annotated"

    io.write_codon_interval_gff(tmp_path, "codons.gff", codons)
    io.write_results_to_gff(result, tmp_path, "orfs", False)
    assert _feature_lines(tmp_path / "codons.gff")[0].split("\t", 1)[0] == "plasmid:alpha-beta%2Ccopy"
    assert _feature_lines(tmp_path / "orfs.gff")[0].split("\t", 1)[0] == "plasmid:alpha-beta%2Ccopy"


def test_repeated_display_labels_keep_distinct_loci_and_densities(tmp_path):
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=first;locus_tag=duplicate\n"
        "chr\t.\tCDS\t20\t28\t.\t+\t0\tID=second;locus_tag=duplicate\n"
    )

    annotation_interlaps, genes = misc.annotation_interlap(annotation)
    assert len(genes) == 2
    assert {record.display_label for record in genes.values()} == {"duplicate"}

    misc.calculate_density(
        {("chr", "+"): {0: 2, 19: 3}}, annotation_interlaps, genes)
    by_start = {record[1]: record for record in genes.values()}
    assert by_start[1][4] == 2
    assert by_start[20][4] == 3

    first_type, first_key = misc.get_gene_information("chr", 0, 8, "+", genes)
    second_type, second_key = misc.get_gene_information("chr", 19, 27, "+", genes)
    assert (first_type, second_type) == ("Annotated", "Annotated")
    assert first_key != second_key
    assert misc.annotation_display_label(first_key, genes) == "duplicate"
    assert misc.annotation_display_label(second_key, genes) == "duplicate"

    genome = {"chr": "ATGCCCTAA" + "C" * 10 + "ATGCCCTAA"}
    detected = {("chr", "+"): {
        (0, 8): (10, np.nan, np.nan),
        (19, 27): (10, np.nan, np.nan),
    }}
    result = misc.generate_result_dataframe(
        detected, genes, {}, {}, genome, {}, [], [], ("TIS", "TTS", "RIBO"))
    assert result["Type"].tolist() == ["Annotated", "Annotated"]
    assert result["Locus_tag"].tolist() == ["duplicate", "duplicate"]
    assert result["TIS_relative_density"].tolist() == [5, 10 / 3]


def test_multipart_cds_groups_phases_and_deduplicates_overlap(tmp_path):
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=joined;locus_tag=multi;"
        "exception=ribosomal%20slippage\n"
        "chr\t.\tCDS\t9\t18\t.\t+\t2\tID=joined;locus_tag=multi;"
        "exception=ribosomal%20slippage\n"
    )

    features = misc.parse_annotation_features(annotation)
    assert len(features) == 1
    assert features[0].segments == ((1, 9), (9, 18))
    assert features[0].phases == ("0", "2")
    assert features[0].exception == "ribosomal slippage"

    annotation_interlaps, genes = misc.annotation_interlap(annotation)
    assert list(genes) == ["multi"]
    assert genes["multi"] == ["chr", 1, 18, "+", 0]
    assert genes["multi"].is_complex is True

    misc.calculate_density(
        {("chr", "+"): {0: 2, 8: 7, 17: 3}}, annotation_interlaps, genes)
    assert genes["multi"][4] == 12

    category, gene_key = misc.get_gene_information("chr", 0, 17, "+", genes)
    assert (category, gene_key) == ("Annotated_Complex", "multi")

    genome = {"chr": "ATG" + "CCC" * 4 + "TAA"}
    detected = {("chr", "+"): {(0, 17): (10, np.nan, np.nan)}}
    result = misc.generate_result_dataframe(
        detected, genes, {}, {}, genome, {}, [], [], ("TIS", "TTS", "RIBO"))
    assert result.loc[0, "Type"] == "Annotated_Complex"
    assert result.loc[0, "Locus_tag"] == "multi"


def test_candidate_inside_multipart_cds_gap_remains_unannotated(tmp_path):
    annotation = tmp_path / "gapped-cds.gff"
    annotation.write_text(
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=joined;locus_tag=multi\n"
        "chr\t.\tCDS\t31\t39\t.\t+\t0\tID=joined;locus_tag=multi\n",
        encoding="utf-8",
    )
    _, genes = misc.annotation_interlap(annotation)

    # The logical CDS bounds are 1..39, but 13..21 does not touch either CDS
    # segment and therefore is not an internal ORF of this annotation.
    assert misc.get_gene_information("chr", 12, 20, "+", genes) == (
        "Unannotated", "chr:13-21:+",
    )
    # A same-frame candidate that really overlaps a segment remains internal.
    assert misc.get_gene_information("chr", 3, 8, "+", genes) == (
        "Internal_Inframe", "multi",
    )


def test_cds_inherits_gene_metadata_through_encoded_transcript_parent(tmp_path):
    annotation = tmp_path / "parent-chain.gff"
    annotation.write_text(
        "chr\t.\tgene\t1\t12\t.\t+\t.\tID=gene-one;locus_tag=LT1;Name=protein-one\n"
        "chr\t.\tmRNA\t1\t12\t.\t+\t.\tID=tx%2Cone;Parent=gene-one\n"
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=cds-one;Parent=tx%2Cone\n",
        encoding="utf-8",
    )

    feature = misc.parse_annotation_features(annotation)[0]

    assert feature.gene_id == "gene-one"
    assert feature.locus_tag == "LT1"
    assert feature.gene_name == "protein-one"
    assert feature.display_label == "LT1"


def test_invalid_multipart_row_cannot_be_hidden_by_valid_group_bounds(tmp_path):
    annotation = tmp_path / "invalid-multipart.gff"
    annotation.write_text(
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=joined\n"
        "chr\t.\tCDS\t20\t10\t.\t+\t0\tID=joined\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="invalid GFF coordinates on line 2"):
        misc.parse_annotation_features(annotation)


@pytest.mark.parametrize(
    "attributes",
    [
        "ID=first;ID=second;locus_tag=gene",
        "ID=first;locus_tag=one;LOCUS_TAG=two",
        "ID=first;Parent=gene-one;Parent=gene-two",
    ],
)
def test_duplicate_gff3_attribute_names_are_rejected(tmp_path, attributes):
    annotation = tmp_path / "duplicate-attributes.gff"
    annotation.write_text(
        f"chr\t.\tCDS\t1\t9\t.\t+\t0\t{attributes}\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="duplicate GFF3 attribute name .*comma-separated value list",
    ):
        misc.parse_annotation_features(annotation)


@pytest.mark.parametrize("attributes", ["=value", "   =value"])
def test_empty_gff3_attribute_name_is_rejected(tmp_path, attributes):
    annotation = tmp_path / "empty-attribute-name.gff"
    annotation.write_text(
        f"chr\t.\tCDS\t1\t9\t.\t+\t0\tID=cds;{attributes}\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="attribute names must not be empty"):
        misc.parse_annotation_features(annotation)


def test_gff3_parent_comma_list_remains_supported(tmp_path):
    annotation = tmp_path / "parent-list.gff"
    annotation.write_text(
        "chr\t.\tgene\t1\t9\t.\t+\t.\tID=gene-one;locus_tag=LT1\n"
        "chr\t.\tgene\t1\t9\t.\t+\t.\tID=gene-two;locus_tag=LT2\n"
        "chr\t.\tCDS\t1\t9\t.\t+\t0\tID=cds;Parent=gene-one,gene-two\n",
        encoding="utf-8",
    )

    feature = misc.parse_annotation_features(annotation)[0]

    assert feature.feature_id == "cds"
    assert feature.locus_tag == "LT1"


@pytest.mark.parametrize("suffix,mode", [(".sam", "w"), (".bam", "wb")])
def test_duplicate_alignment_reference_names_are_rejected(
    tmp_path, monkeypatch, suffix, mode,
):
    annotation = tmp_path / "annotation.gff"
    annotation.write_text("##gff-version 3\n", encoding="utf-8")
    alignment = tmp_path / f"duplicates{suffix}"
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr", "LN": 9}, {"SN": "chr", "LN": 9}],
    }
    if suffix == ".sam":
        alignment.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr\tLN:9\n@SQ\tSN:chr\tLN:9\n",
            encoding="utf-8",
        )
    else:
        with pysam.AlignmentFile(alignment, mode, header=header):
            pass

    def fail_if_pysam_parses_duplicate(*args, **kwargs):
        pytest.fail("duplicate reference preflight must run before pysam")

    monkeypatch.setattr(io.pysam, "AlignmentFile", fail_if_pysam_parses_duplicate)
    with pytest.raises(ValueError, match="duplicate reference name 'chr'"):
        io.validate_reference_inputs({"chr": "ATGCCCTAA"}, annotation, [alignment])


def test_duplicate_alignment_preflight_defers_corrupt_gzip(tmp_path):
    alignment = tmp_path / "corrupt.sam"
    corrupt_member = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xffBADDEFLATE"
    alignment.write_bytes(corrupt_member * 2)

    assert io._duplicate_alignment_reference_name(alignment) is None


def test_duplicate_alignment_preflight_has_a_scan_budget(tmp_path, monkeypatch):
    alignment = tmp_path / "oversized.sam"
    alignment.write_bytes(b"@CO\t" + b"x" * 64)
    monkeypatch.setattr(io, "_MAX_ALIGNMENT_HEADER_SCAN_BYTES", 32)

    assert io._duplicate_alignment_reference_name(alignment) is None
