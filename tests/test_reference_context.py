"""Batch reference-context reuse without sharing mutable analysis state."""

from dataclasses import FrozenInstanceError

import pandas as pd
import pytest

import orfbounder
from helpers.toy_example_generation import ToyExampleGenerator
from lib import config, io, misc, run_output
from orfbounder_batch import call_orfbounder


@pytest.fixture
def example(tmp_path):
    generated = ToyExampleGenerator(tmp_path / "example")
    generated.generate_all()
    return generated


def _run_sample(example, output_path, *, reference_cache=None):
    return orfbounder.run_orfbounder(
        example.bam_dir / "TIS-WT-1.bam",
        None,
        example.bam_dir / "RIBO-WT-1.bam",
        example.config_dir / "read_lengths.json",
        "raw",
        "fiveprime",
        example.data_dir / "annotation.gff",
        example.data_dir / "genome.fa",
        ["ATG"],
        ["TAG", "TAA", "TGA"],
        output_path,
        "WT-1",
        example.config_dir / "offsets.json",
        "furthest_inframe",
        1,
        "sum",
        True,
        None,
        None,
        _reference_cache=reference_cache,
    )[0]


def test_cached_reference_matches_uncached_results_and_densities_do_not_leak(
        example, tmp_path):
    expected = _run_sample(example, tmp_path / "uncached")
    reference_cache = {}
    first = _run_sample(example, tmp_path / "cached-first", reference_cache=reference_cache)
    second = _run_sample(example, tmp_path / "cached-second", reference_cache=reference_cache)

    pd.testing.assert_frame_equal(first, expected)
    pd.testing.assert_frame_equal(second, expected)
    assert len(reference_cache) == 1
    assert first.filter(like="relative_density").notna().any().any()

    context = next(iter(reference_cache.values()))
    assert isinstance(context.annotation_features, tuple)
    assert isinstance(context.circular_contigs, frozenset)
    with pytest.raises(FrozenInstanceError):
        context.circular_contigs = frozenset({"changed"})

    first_interlaps, first_genes = misc.annotation_interlap(
        example.data_dir / "annotation.gff",
        annotation_features=context.annotation_features,
    )
    _, fresh_genes = misc.annotation_interlap(
        example.data_dir / "annotation.gff",
        annotation_features=context.annotation_features,
    )
    misc.calculate_density(
        {(example.chrom, "+"): {99: 7}}, first_interlaps, first_genes,
    )
    assert any(gene[4] > 0 for gene in first_genes.values())
    assert all(gene[4] == 0 for gene in fresh_genes.values())
    assert all(first_genes[key] is not fresh_genes[key] for key in first_genes)


def test_batch_grid_parses_one_reference_for_two_samples_and_two_cells(
        example, tmp_path, monkeypatch):
    config_path = example.config_dir / "config.tsv"
    frame = pd.read_csv(config_path, sep="\t", keep_default_na=False)
    frame["normalization_method"] = "raw,mil"
    # Reference caching is independent of expression counting; omit those extra
    # alignment passes to keep this regression focused and fast.
    frame["alignment_folder_path"] = ""
    frame.to_csv(config_path, sep="\t", index=False)
    validated = config.check_config_sheet(config_path)

    calls = {"genome": 0, "topology": 0, "features": 0}

    def counted(name, function):
        def wrapper(*args, **kwargs):
            calls[name] += 1
            return function(*args, **kwargs)
        return wrapper

    monkeypatch.setattr(
        io, "generate_genome_dict", counted("genome", io.generate_genome_dict),
    )
    monkeypatch.setattr(
        misc, "parse_reference_topology",
        counted("topology", misc.parse_reference_topology),
    )
    monkeypatch.setattr(
        misc, "parse_annotation_features",
        counted("features", misc.parse_annotation_features),
    )

    output = tmp_path / "results"
    call_orfbounder(validated, output)

    assert calls == {"genome": 1, "topology": 1, "features": 1}
    for normalization in ("raw", "mil"):
        destination = output / "toy_example" / "fiveprime" / normalization
        complete = pd.read_json(destination / "complete.json", typ="series")
        assert complete["samples"] == ["WT-1", "WT-2"]
        merged = pd.read_csv(destination / "toy_example_final.csv", sep="\t")
        assert len(merged) == 3


def test_config_validation_parses_a_shared_reference_once(
        example, monkeypatch):
    config_path = example.config_dir / "config.tsv"
    first = pd.read_csv(config_path, sep="\t", dtype=str, keep_default_na=False)
    second = first.copy()
    second["experiment_name"] = "second_experiment"
    pd.concat([first, second], ignore_index=True).to_csv(
        config_path, sep="\t", index=False,
    )

    calls = {"genome": 0, "topology": 0, "features": 0}

    def counted(name, function):
        def wrapper(*args, **kwargs):
            calls[name] += 1
            return function(*args, **kwargs)
        return wrapper

    monkeypatch.setattr(
        io, "generate_genome_dict", counted("genome", io.generate_genome_dict),
    )
    monkeypatch.setattr(
        misc, "parse_reference_topology",
        counted("topology", misc.parse_reference_topology),
    )
    monkeypatch.setattr(
        misc, "parse_annotation_features",
        counted("features", misc.parse_annotation_features),
    )

    validated = config.check_config_sheet(config_path)

    assert len(validated) == 2
    assert calls == {"genome": 1, "topology": 1, "features": 1}


def test_reference_context_cache_invalidates_on_sha_or_path_change(example, tmp_path):
    genome = example.data_dir / "genome.fa"
    annotation = example.data_dir / "annotation.gff"
    fingerprint_cache = {}
    reference_cache = {}

    records = run_output.fingerprint_inputs(
        [genome, annotation], cache=fingerprint_cache,
    )
    original = io.load_reference_context(
        genome, annotation, records, cache=reference_cache,
    )
    assert io.load_reference_context(
        genome, annotation, records, cache=reference_cache,
    ) is original

    annotation.write_text(
        annotation.read_text(encoding="utf-8") + "# content-key change\n",
        encoding="utf-8",
    )
    changed_records = run_output.fingerprint_inputs(
        [genome, annotation], cache=fingerprint_cache,
    )
    changed = io.load_reference_context(
        genome, annotation, changed_records, cache=reference_cache,
    )
    assert changed is not original
    assert changed.cache_key[1][1] != original.cache_key[1][1]

    copied_annotation = tmp_path / "copied-annotation.gff"
    copied_annotation.write_bytes(annotation.read_bytes())
    copied_records = run_output.fingerprint_inputs(
        [genome, copied_annotation], cache=fingerprint_cache,
    )
    copied = io.load_reference_context(
        genome, copied_annotation, copied_records, cache=reference_cache,
    )
    assert copied is not changed
    assert copied.cache_key[1][0] != changed.cache_key[1][0]
    assert copied.cache_key[1][1] == changed.cache_key[1][1]
    assert len(reference_cache) == 3
