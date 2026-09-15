"""Circular-reference regressions across calling, annotation, counts and output."""

import json
from pathlib import Path

import pandas as pd
import pysam
import pytest
from interlap import InterLap

import orfbounder
from helpers.final_to_gff import create_gff
from lib import expression, io, merging, misc, predictions, run_output
from lib.alignment_reader import PositionReader
from lib.coordinates import parse_coordinate_id
from lib.matched_statistics import ATTACH_FIELDS, attach_matched_to_orfs
from lib.statistics import analyze_peak_enrichment, attach_to_orfs


PLUS_SEQUENCE = "TGCCCTAACCCA"   # ATG at 11..1; in-frame TAA at 5..7.
MINUS_SEQUENCE = "ATCCCTTACCCC"  # reverse ATG at 11..1; reverse TAA at 5..7.


def _activated_codons(sequence, codons, strand):
    _, candidates = misc.create_codon_interlaps("circ", sequence, codons, circular=True)
    for identifier, value in candidates.items():
        value[1] = 7 if parse_coordinate_id(identifier)[3] == strand else 0
    return candidates


@pytest.mark.parametrize(
    "method,strand,sequence,search_codons,match_codons,expected",
    [
        ("TIS", "+", PLUS_SEQUENCE, ["ATG"], ["TAA"], (11, 19)),
        ("RIBO", "+", PLUS_SEQUENCE, ["ATG"], ["TAA"], (11, 19)),
        ("TIS", "-", MINUS_SEQUENCE, ["ATG"], ["TAA"], (5, 13)),
        ("TTS", "+", PLUS_SEQUENCE, ["TAA"], ["ATG"], (11, 19)),
        ("TTS", "-", MINUS_SEQUENCE, ["TAA"], ["ATG"], (5, 13)),
    ],
)
def test_circular_calling_finds_origin_spanning_orfs_on_both_strands(
        method, strand, sequence, search_codons, match_codons, expected):
    candidates = _activated_codons(sequence, search_codons, strand)
    result = predictions.detect_potential_orfs(
        candidates, sequence, search_codons, match_codons, method, {},
        "next_inframe", circular=True,
    )
    assert list(result[("circ", strand)]) == [expected]
    peak_index = {"TIS": 0, "TTS": 1, "RIBO": 2}[method]
    assert result[("circ", strand)][expected][peak_index] == 7


@pytest.mark.parametrize(
    "sequence,strand,next_expected,furthest_expected",
    [
        ("CCCATGTAACCCATG", "+", (3, 8), (12, 23)),
        ("CATGGGTTACATGGG", "-", (6, 11), (6, 17)),
    ],
)
def test_circular_furthest_tts_start_stops_before_one_revolution(
        sequence, strand, next_expected, furthest_expected):
    candidates = _activated_codons(sequence, ["TAA"], strand)
    closest = predictions.detect_potential_orfs(
        candidates, sequence, ["TAA"], ["ATG"], "TTS", {},
        "next_inframe", circular=True,
    )
    furthest = predictions.detect_potential_orfs(
        candidates, sequence, ["TAA"], ["ATG"], "TTS", {},
        "furthest_inframe", circular=True,
    )
    assert set(closest[("circ", strand)]) == {next_expected}
    assert set(furthest[("circ", strand)]) == {furthest_expected}
    assert furthest[("circ", strand)].keys() != closest[("circ", strand)].keys()


def test_boundary_spanning_codon_window_wraps_but_linear_mode_does_not():
    interlaps, candidates = misc.create_codon_interlaps(
        "circ", PLUS_SEQUENCE, ["ATG"], circular=True,
    )
    plus = next(key for key in candidates if parse_coordinate_id(key)[3] == "+")
    assert parse_coordinate_id(plus)[1:3] == (9, 13)
    predictions.screen_positions_for_tss(
        {("circ", "+"): {0: 4, 11: 3}}, interlaps, candidates, 0, "sum",
    )
    assert candidates[plus][1] == 7

    _, linear = misc.create_codon_interlaps("circ", PLUS_SEQUENCE, ["ATG"])
    assert not any(parse_coordinate_id(key)[3] == "+" for key in linear)


def test_circular_search_without_a_stop_is_bounded_on_non_triplet_length():
    sequence = "TGCCCCCCCA"  # length 10; wrapped ATG at position 9, no TAA.
    candidates = _activated_codons(sequence, ["ATG"], "+")
    result = predictions.detect_potential_orfs(
        candidates, sequence, ["ATG"], ["TAA"], "TIS", {},
        "furthest_inframe", circular=True,
    )
    assert result == {}


def _write_endpoint_sam(path: Path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "circ", "LN": 10}]}
    with pysam.AlignmentFile(path, "w", header=header) as handle:
        for name, position, flag in (("plus", 1, 0), ("minus", 8, 16)):
            read = pysam.AlignedSegment(handle.header)
            read.query_name = name
            read.reference_id = 0
            read.reference_start = position
            read.flag = flag
            read.mapping_quality = 60
            read.cigarstring = "1M"
            read.query_sequence = "A"
            handle.write(read)


def test_offsets_wrap_only_on_declared_circular_references(tmp_path):
    alignment = tmp_path / "TIS-x-1.sam"
    _write_endpoint_sam(alignment)
    offsets = {"default": {"default": 3}}
    circular = PositionReader(
        alignment, None, "fiveprime", offsets, circular_contigs={"circ"},
    )
    linear = PositionReader(alignment, None, "fiveprime", offsets)
    assert circular.reads_position_dict == {
        ("circ", "+"): {8: 1}, ("circ", "-"): {1: 1},
    }
    assert linear.reads_position_dict == {}
    assert circular.no_accepted_reads_dict == linear.no_accepted_reads_dict == {"circ": 2}


def test_centered_short_reads_remain_in_normalization_denominator(tmp_path):
    alignment = tmp_path / "TIS-x-1.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "circ", "LN": 30}]}
    with pysam.AlignmentFile(alignment, "w", header=header) as handle:
        read = pysam.AlignedSegment(handle.header)
        read.query_name = "short"
        read.reference_id = 0
        read.reference_start = 5
        read.mapping_quality = 60
        read.cigarstring = "10M"
        read.query_sequence = "A" * 10
        handle.write(read)
    reader = PositionReader(
        alignment, None, "centered", {"default": {"default": 0}},
    )
    assert reader.reads_position_dict == {}
    assert reader.no_accepted_reads_dict == {"circ": 1}


def _write_multipart_annotation(path: Path, seqid="circ", high_start=90):
    encoded = seqid.replace(":", "%3A")
    path.write_text(
        "##gff-version 3\n"
        f"##sequence-region {encoded} 1 100\n"
        f"{encoded}\ttest\tregion\t1\t100\t.\t.\t.\tID=region;Is_circular=true\n"
        f"{encoded}\ttest\tCDS\t{high_start}\t100\t.\t+\t0\tID=plus;locus_tag=repeat\n"
        f"{encoded}\ttest\tCDS\t1\t20\t.\t+\t1\tID=plus;locus_tag=repeat\n"
        f"{encoded}\ttest\tCDS\t{high_start}\t100\t.\t-\t2\tID=minus;locus_tag=repeat\n"
        f"{encoded}\ttest\tCDS\t1\t20\t.\t-\t0\tID=minus;locus_tag=repeat\n",
        encoding="utf-8",
    )


def test_circular_multipart_annotation_is_unwrapped_strand_aware_and_counted(tmp_path):
    annotation = tmp_path / "annotation.gff"
    _write_multipart_annotation(annotation, "circ:part-1")
    features = misc.parse_annotation_features(annotation)
    plus, minus = features
    assert (plus.chromosome, plus.start, plus.stop, plus.segments, plus.phases) == (
        "circ:part-1", 90, 120, ((90, 100), (1, 20)), ("0", "1"),
    )
    assert (minus.start, minus.stop, minus.segments, minus.phases) == (
        90, 120, ((1, 20), (90, 100)), ("0", "2"),
    )
    assert plus.crosses_origin and minus.crosses_origin

    annotation_interlaps, genes = misc.annotation_interlap(annotation)
    densities = misc.calculate_density(
        {
            ("circ:part-1", "+"): {89: 1, 0: 2, 19: 3},
            ("circ:part-1", "-"): {89: 1, 0: 2, 19: 3},
        },
        annotation_interlaps,
        genes,
    )
    assert sorted(record[4] for record in densities.values()) == [6, 6]
    for key, record in densities.items():
        category, matched = misc.get_gene_information(
            "circ:part-1", 89, 119, record[3], genes,
        )
        assert category == "Annotated_Complex"
        assert matched == key


@pytest.mark.parametrize(
    "strand,start,stop,expected",
    [
        ("+", 95, 119, "Annotated_Complex"),
        ("-", 95, 119, "Annotated_Complex"),
        ("+", 1, 19, "Near_Annotated"),
        ("-", 95, 113, "Near_Annotated"),
        ("+", 1, 16, "Internal_Inframe"),
        ("-", 0, 13, "Internal_Inframe"),
        ("+", 2, 16, "Internal_OutofFrame"),
        ("-", 0, 12, "Internal_OutofFrame"),
    ],
)
def test_circular_annotation_relationships_align_adjacent_virtual_copies(
        tmp_path, strand, start, stop, expected):
    annotation = tmp_path / "relationships.gff"
    _write_multipart_annotation(annotation, high_start=96)
    _, all_genes = misc.annotation_interlap(annotation)
    genes = {
        key: record for key, record in all_genes.items()
        if record[3] == strand
    }

    category, key = misc.get_gene_information(
        "circ", start, stop, strand, genes,
    )

    assert category == expected
    assert key in genes

    index = misc._build_annotation_candidate_index(genes)
    candidates = misc._overlapping_annotation_genes(
        index, genes, "circ", start, stop, strand,
    )
    assert misc.get_gene_information(
        "circ", start, stop, strand, candidates,
    ) == (category, key)


@pytest.mark.parametrize(
    "strand,start,stop",
    [
        ("+", 1, 19),   # Physical 2..20 is virtual 102..120 beside the CDS.
        ("-", 0, 13),   # Physical 1..14 is virtual 101..114 beside the CDS.
    ],
)
def test_circular_annotation_distances_use_the_matching_virtual_copy(
        tmp_path, strand, start, stop):
    annotation = tmp_path / "distances.gff"
    _write_multipart_annotation(annotation, high_start=96)
    _, all_genes = misc.annotation_interlap(annotation)
    key, record = next(
        (key, record) for key, record in all_genes.items()
        if record[3] == strand
    )

    assert misc.calculate_utr_distance(start, stop, key, {key: record}) == (6, 18)


def test_circular_extension_can_align_to_the_previous_virtual_copy():
    feature = misc.AnnotationFeature(
        feature_id="low", chromosome="circ", start=1, stop=20, strand="+",
        gene_id="low", locus_tag="low", name="", read_list=(),
        gene_name="", old_locus_tag="", segments=((1, 20),), phases=("0",),
        reference_length=100, is_circular=True,
    )
    genes = {"low": misc.AnnotationGene(feature)}

    category, key = misc.get_gene_information("circ", 89, 119, "+", genes)
    index = misc._build_annotation_candidate_index(genes)
    candidates = misc._overlapping_annotation_genes(
        index, genes, "circ", 89, 119, "+",
    )

    assert (category, key) == ("N-terminal_extension", "low")
    assert misc.get_gene_information(
        "circ", 89, 119, "+", candidates,
    ) == (category, key)
    assert misc.calculate_utr_distance(89, 119, key, genes) == (-11, 30)


def test_circular_annotation_topology_is_validated_against_fasta(tmp_path):
    annotation = tmp_path / "annotation.gff"
    _write_multipart_annotation(annotation)
    io.validate_reference_inputs({"circ": "A" * 100}, annotation, [])
    with pytest.raises(ValueError, match="length 100"):
        io.validate_reference_inputs({"circ": "A" * 99}, annotation, [])

    invalid = tmp_path / "virtual-input.gff"
    invalid.write_text(
        "##gff-version 3\n"
        "circ\ttest\tregion\t1\t12\t.\t.\t.\tID=r;Is_circular=true\n"
        "circ\ttest\tCDS\t12\t20\t.\t+\t0\tID=bad\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="physical coordinates"):
        misc.parse_annotation_features(invalid)


def test_conflicting_whole_reference_circularity_is_rejected(tmp_path):
    annotation = tmp_path / "conflicting-topology.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "##sequence-region circ 1 12\n"
        "circ\ttest\tregion\t1\t12\t.\t.\t.\tID=linear;Is_circular=false\n"
        "circ\ttest\tregion\t1\t12\t.\t.\t.\tID=circle;Is_circular=true\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Conflicting Is_circular declarations"):
        misc.parse_reference_topology(annotation)


def test_duplicate_is_circular_attribute_on_one_region_is_rejected(tmp_path):
    annotation = tmp_path / "duplicate-circular-attribute.gff"
    annotation.write_text(
        "##gff-version 3\n"
        "##sequence-region circ 1 12\n"
        "circ\ttest\tregion\t1\t12\t.\t.\t.\t"
        "ID=circle;Is_circular=true;is_circular=false\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate GFF3 attribute name"):
        misc.parse_reference_topology(annotation)


@pytest.mark.parametrize(
    "start,stop,strand,sequence,expected",
    [
        (11, 19, "+", PLUS_SEQUENCE, "ATGCCCTAA"),
        (5, 13, "-", MINUS_SEQUENCE, "ATGGGGTAA"),
    ],
)
def test_circular_sequence_and_upstream_window_use_physical_bases_once(
        start, stop, strand, sequence, expected):
    nt, _, upstream, start_codon, stop_codon = misc.get_genome_information(
        start, stop, strand, sequence, circular=True,
    )
    assert nt == expected
    assert (start_codon, stop_codon) == ("ATG", "TAA")
    assert len(upstream) == len(sequence)  # Reference is smaller than 15 nt.


def test_expression_query_wraps_and_deduplicates_alignment_records():
    intervals = {
        ("circ", "+"): InterLap([
            (10, 11, "spanning"), (0, 1, "spanning"), (2, 2, "other"),
        ])
    }
    assert expression.count_reads(
        "circ", 9, 14, "+", intervals,
        reference_lengths={"circ": 12}, circular_contigs={"circ"},
    ) == 2


def test_circular_local_background_uses_every_physical_base_at_most_once():
    counts = {("circ", "+"): {position: 1 for position in range(12)}}
    result = analyze_peak_enrichment(
        {"circ": PLUS_SEQUENCE}, counts, ["ATG"], flank_width=4,
        circular_contigs={"circ"},
    )
    plus = result.loc[(result.Strand == "+") & (result.Codon_start == 12)].iloc[0]
    assert (plus.Window_start, plus.Window_stop) == (10, 14)
    assert (plus.peak_count, plus.peak_width) == (5, 5)
    assert (plus.background_count, plus.background_width) == (7, 7)
    assert plus.peak_width + plus.background_width == 12

    tiny = analyze_peak_enrichment(
        {"tiny": "ATG"}, {("tiny", "+"): {0: 2, 1: 3, 2: 4}},
        ["ATG"], flank_width=20, circular_contigs={"tiny"},
    ).iloc[0]
    assert (tiny.peak_count, tiny.peak_width, tiny.background_count, tiny.background_width) == (9, 3, 0, 0)


@pytest.mark.parametrize(
    "method,strand,codon_start,start,stop",
    [
        ("TIS", "+", 12, 12, 20),
        ("TIS", "-", 12, 6, 14),
        ("TTS", "+", 6, 12, 20),
        ("TTS", "-", 6, 6, 14),
    ],
)
def test_local_and_matched_attachment_compare_circular_boundaries_modulo(
        method, strand, codon_start, start, stop):
    topology = {"Reference_length": 12, "Is_circular": True}
    orfs = pd.DataFrame([{
        "Genome": "circ", "Start": start, "Stop": stop, "Strand": strand, **topology,
    }])
    local = pd.DataFrame([{
        "Genome": "circ", "Strand": strand, "Codon_start": codon_start,
        "peak_count": 9, "background_count": 1, "p_value": 0.01,
        "q_value": 0.02, "significant": True, "correction": "by",
        "tests_in_family": 4, **topology,
    }])
    attached = attach_to_orfs(orfs, local, method, "sample")
    assert attached.loc[0, "sample_peak_count"] == 9

    matched_values = {field: pd.NA for field in ATTACH_FIELDS}
    matched_values.update({"matched_pair_count": 2, "q_value": 0.02, "significant": True})
    matched = pd.DataFrame([{
        "Genome": "circ", "Strand": strand, "Codon_start": codon_start,
        "assay": method, **topology, **matched_values,
    }])
    matched_attached = attach_matched_to_orfs(orfs, matched, method, "comparison")
    assert matched_attached.loc[0, "comparison_matched_pair_count"] == 2


def _minimal_result_frame():
    return pd.DataFrame([{
        "Type": "Annotated", "Identifier": "circ:part-1:12-20:+",
        "Genome": "circ:part-1", "Start": 12, "Stop": 20, "Strand": "+",
        "Locus_tag": "gene", "Codon_count": 3,
        "TIS_peak_height": 1, "TTS_peak_height": pd.NA, "RIBO_peak_height": pd.NA,
        "Start_codon": "ATG", "Stop_codon": "TAA", "15nt_window": "A" * 12,
        "Nucleotide_Seq": "ATGCCCTAA", "Amino_Acid_Seq": "MP*",
        "TIS_relative_density": 1.0, "TTS_relative_density": pd.NA,
        "RIBO_relative_density": pd.NA, "5'-distance": 0, "3'-distance": 8,
        "Reference_length": 12, "Is_circular": True,
    }])


def test_circular_topology_and_virtual_coordinates_survive_gff_and_merge(tmp_path):
    frame = _minimal_result_frame()
    io.write_results_to_gff(frame, tmp_path, "sample", False)
    sample_lines = (tmp_path / "sample.gff").read_text().splitlines()
    assert "##sequence-region circ:part-1 1 12" in sample_lines
    assert any("\tregion\t1\t12\t" in line and "Is_circular=true" in line for line in sample_lines)
    sample_cds = [line.split("\t") for line in sample_lines if "\tCDS\t" in line]
    assert [(fields[3], fields[4], fields[7]) for fields in sample_cds] == [
        ("12", "12", "0"), ("1", "8", "2"),
    ]

    metadata, dynamic = merging.extend_combined_dictionary(frame, {}, {})
    merged, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert merged.loc[0, ["Reference_length", "Is_circular"]].tolist() == [12, True]
    merging.write_merged_gff(merged, tmp_path / "merged.xlsx")
    merged_text = (tmp_path / "merged.gff").read_text()
    assert "##sequence-region circ:part-1 1 12" in merged_text
    assert "\tCDS\t12\t12\t" in merged_text
    assert "\tCDS\t1\t8\t" in merged_text


def test_split_gff_phases_follow_biological_strand_order():
    assert io.gff3_feature_segments(12, 20, "+", 12, True) == (
        (12, 12, "0"), (1, 8, "2"),
    )
    assert io.gff3_feature_segments(12, 20, "-", 12, True) == (
        (1, 8, "0"), (12, 12, "1"),
    )


def test_gff_writers_reject_conflicting_topology_on_one_contig(tmp_path):
    frame = pd.concat([_minimal_result_frame(), _minimal_result_frame()], ignore_index=True)
    frame.loc[1, "Is_circular"] = False
    with pytest.raises(ValueError, match="Conflicting reference topology"):
        io.write_results_to_gff(frame, tmp_path, "conflict", False)
    with pytest.raises(ValueError, match="Conflicting reference topology"):
        create_gff(frame, tmp_path / "conflict-helper.gff")


def test_empty_batch_merge_retains_circular_topology_for_gff(tmp_path):
    empty = _minimal_result_frame().iloc[:0].copy()
    empty.attrs["reference_lengths"] = {"circ:part-1": 12}
    empty.attrs["circular_contigs"] = ["circ:part-1"]
    metadata, dynamic = merging.extend_combined_dictionary(empty, {}, {})
    merged, _ = merging.build_merged_dataframe(metadata, dynamic)
    assert merged.empty
    assert merged.attrs["circular_contigs"] == ["circ:part-1"]
    merging.write_merged_gff(merged, tmp_path / "empty-merged.xlsx")
    text = (tmp_path / "empty-merged.gff").read_text()
    assert "##sequence-region circ:part-1 1 12" in text
    assert "Is_circular=true" in text


def _write_circular_run_inputs(tmp_path):
    seqid = "circ:part-1"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">{seqid}\n{PLUS_SEQUENCE}\n", encoding="utf-8")
    annotation = tmp_path / "annotation.gff"
    annotation.write_text(
        "##gff-version 3\n"
        f"##sequence-region {seqid} 1 12\n"
        f"{seqid}\ttest\tregion\t1\t12\t.\t.\t.\tID=region;Is_circular=true\n"
        f"{seqid}\ttest\tCDS\t12\t12\t.\t+\t0\tID=gene;locus_tag=gene\n"
        f"{seqid}\ttest\tCDS\t1\t8\t.\t+\t2\tID=gene;locus_tag=gene\n",
        encoding="utf-8",
    )
    alignment = tmp_path / "TIS-circ-1.sam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": seqid, "LN": 12}]}
    with pysam.AlignmentFile(alignment, "w", header=header) as handle:
        read = pysam.AlignedSegment(handle.header)
        read.query_name = "peak"
        read.reference_id = 0
        read.reference_start = 11
        read.mapping_quality = 60
        read.cigarstring = "1M"
        read.query_sequence = "A"
        handle.write(read)
    offsets = tmp_path / "offsets.json"
    offsets.write_text('{"default":{"default":0}}', encoding="utf-8")
    return seqid, fasta, annotation, alignment, offsets


def test_direct_circular_run_records_topology_and_empty_exports_keep_it(tmp_path):
    seqid, fasta, annotation, alignment, offsets = _write_circular_run_inputs(tmp_path)
    output = tmp_path / "result"
    frame, candidates = orfbounder.run_orfbounder(
        alignment, None, None, None, "raw", "fiveprime", annotation, fasta,
        ["ATG"], ["TAA"], output, "circ-run", offsets,
        "furthest_inframe", 0, "sum", True, None, None, statistics="local",
    )
    assert frame[["Start", "Stop", "Nucleotide_Seq"]].iloc[0].tolist() == [12, 20, "ATGCCCTAA"]
    assert frame.loc[0, "Is_circular"]
    assert set(candidates["Reference_length"]) == {12}
    assert candidates["Is_circular"].all()
    manifest = json.loads((output / "circ-run.run.json").read_text())
    assert manifest["reference_lengths"] == {seqid: 12}
    assert manifest["circular_contigs"] == [seqid]
    assert "scipy" in manifest["dependencies"]

    run_output.export_sample_results(frame, output, "circ-run")
    gff_text = (output / "gff_per_sample/circ-run.gff").read_text()
    assert "\tCDS\t12\t12\t" in gff_text
    assert "\tCDS\t1\t8\t" in gff_text

    empty = frame.iloc[:0].copy()
    io.write_results_to_gff(empty, output, "empty-circular", False)
    empty_text = (output / "empty-circular.gff").read_text()
    assert f"##sequence-region {seqid} 1 12" in empty_text
    assert "Is_circular=true" in empty_text
    assert "\tCDS\t" not in empty_text
