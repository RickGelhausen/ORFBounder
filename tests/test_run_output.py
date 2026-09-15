"""Failure recovery and portable sequence exports for complete analysis runs."""

import hashlib
import errno
import json
import os
from pathlib import Path
from copy import deepcopy

import pandas as pd
import pytest

from lib import run_output


def test_staging_publishes_completion_last(tmp_path, monkeypatch):
    published = []
    original_link = run_output.os.link

    def record_link(source, target, *args, **kwargs):
        published.append(Path(target).name)
        original_link(source, target, *args, **kwargs)

    monkeypatch.setattr(run_output.os, "link", record_link)
    destination = tmp_path / "results"
    with run_output.staged_output(destination) as stage:
        (stage / "tables").mkdir()
        (stage / "tables/example.tsv").write_text("some results")
        (stage / "complete.json").write_text('{"status":"complete"}')
        assert not (destination / "tables/example.tsv").exists()
    assert published[-1] == "complete.json"
    assert (destination / "tables/example.tsv").read_text() == "some results"
    assert not list(destination.glob(".orfbounder-stage-*"))


def test_export_exception_leaves_no_result_files(tmp_path):
    with pytest.raises(OSError, match="disk full"):
        with run_output.staged_output(tmp_path) as stage:
            (stage / "table.tsv").write_text("partial")
            raise OSError("disk full")
    assert list(tmp_path.iterdir()) == []


def test_failed_staging_removes_only_its_new_empty_destination(tmp_path):
    new_destination = tmp_path / "new"
    with pytest.raises(OSError, match="analysis failed"):
        with run_output.staged_output(new_destination):
            raise OSError("analysis failed")
    assert not new_destination.exists()

    existing_destination = tmp_path / "existing"
    existing_destination.mkdir()
    with pytest.raises(OSError, match="analysis failed"):
        with run_output.staged_output(existing_destination):
            raise OSError("analysis failed")
    assert existing_destination.is_dir()
    assert not list(existing_destination.iterdir())


def test_failed_staging_preserves_concurrently_added_content(tmp_path):
    destination = tmp_path / "results"
    with pytest.raises(OSError, match="analysis failed"):
        with run_output.staged_output(destination):
            (destination / "user-note.txt").write_text("keep")
            raise OSError("analysis failed")
    assert (destination / "user-note.txt").read_text() == "keep"
    assert {path.name for path in destination.iterdir()} == {"user-note.txt"}


def test_failed_staging_preserves_replaced_destination(tmp_path):
    destination = tmp_path / "results"
    moved = tmp_path / "moved"
    with pytest.raises(OSError, match="analysis failed"):
        with run_output.staged_output(destination):
            destination.rename(moved)
            destination.mkdir()
            raise OSError("analysis failed")
    assert destination.is_dir()
    assert moved.is_dir()


def test_publication_failure_rolls_back_own_files(tmp_path, monkeypatch):
    (tmp_path / "keep.txt").write_text("user content")
    original_link = run_output.os.link
    count = 0

    def failing_link(source, target, *args, **kwargs):
        nonlocal count
        count += 1
        if count == 2:
            raise OSError("publication failed")
        original_link(source, target, *args, **kwargs)

    monkeypatch.setattr(run_output.os, "link", failing_link)
    with pytest.raises(OSError, match="publication failed"):
        with run_output.staged_output(tmp_path) as stage:
            (stage / "new").mkdir()
            (stage / "new/a.tsv").write_text("a")
            (stage / "new/b.tsv").write_text("b")
            (stage / "complete.json").write_text("complete")
    assert {path.name for path in tmp_path.iterdir()} == {"keep.txt"}
    assert (tmp_path / "keep.txt").read_text() == "user content"


def test_publication_falls_back_to_exclusive_copy_when_links_are_unsupported(tmp_path, monkeypatch):
    def unsupported_link(*args, **kwargs):
        raise OSError(errno.EOPNOTSUPP, "hard links unsupported")

    monkeypatch.setattr(run_output.os, "link", unsupported_link)
    with run_output.staged_output(tmp_path) as stage:
        (stage / "table.tsv").write_text("complete results")
        (stage / "complete.json").write_text("complete")

    assert (tmp_path / "table.tsv").read_text() == "complete results"
    assert (tmp_path / "complete.json").read_text() == "complete"
    assert not list(tmp_path.glob(".orfbounder-stage-*"))


def test_copy_fallback_opens_staged_source_without_blocking(
        tmp_path, monkeypatch):
    source = tmp_path / "source.tsv"
    source.write_text("results")
    destination = tmp_path / "destination"
    destination.mkdir()
    original_open = run_output.os.open
    checked_source = False

    def unsupported_link(*args, **kwargs):
        raise OSError(errno.EOPNOTSUPP, "hard links unsupported")

    def checked_open(path, flags, *args, **kwargs):
        nonlocal checked_source
        if Path(path) == source:
            checked_source = True
            assert flags & os.O_NONBLOCK
        return original_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(run_output.os, "link", unsupported_link)
    monkeypatch.setattr(run_output.os, "open", checked_open)

    run_output._publish_exclusive(source, destination, Path("result.tsv"))

    assert checked_source
    assert (destination / "result.tsv").read_text() == "results"


def test_fallback_publication_journals_each_copied_target(tmp_path, monkeypatch):
    """Every copied file needs an identity so interrupted runs are recoverable."""
    def unsupported_link(*args, **kwargs):
        raise OSError(errno.EOPNOTSUPP, "hard links unsupported")

    journals = []
    original_write_journal = run_output._write_transaction_journal

    def record_journal(stage, entries):
        journals.append(deepcopy(entries))
        original_write_journal(stage, entries)

    monkeypatch.setattr(run_output.os, "link", unsupported_link)
    monkeypatch.setattr(run_output, "_write_transaction_journal", record_journal)
    with run_output.staged_output(tmp_path) as stage:
        (stage / "a.tsv").write_text("a")
        (stage / "b.tsv").write_text("b")
        (stage / "complete.json").write_text("complete")

    assert {entry["relative"] for entry in journals[-1] if "target_identity" in entry} == {
        "a.tsv", "b.tsv", "complete.json",
    }
    assert all(entry["publication_complete"] for entry in journals[-1])


def test_existing_output_prevents_any_publication(tmp_path):
    (tmp_path / "z.tsv").write_text("old")
    with pytest.raises(FileExistsError):
        with run_output.staged_output(tmp_path) as stage:
            (stage / "a.tsv").write_text("new")
            (stage / "z.tsv").write_text("replace")
    assert not (tmp_path / "a.tsv").exists()
    assert (tmp_path / "z.tsv").read_text() == "old"


def _stale_stage(destination, files):
    stage = destination / ".orfbounder-stage-interrupted"
    stage.mkdir()
    (stage / ".transaction.lock").touch()
    entries = []
    for relative, contents, completion in files:
        source = stage / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(contents)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        os.link(source, target)
        entries.append({"relative": relative, "completion": completion})
    run_output._write_transaction_journal(stage, entries)
    return stage


def test_next_run_recovers_provably_owned_files_from_interrupted_publication(tmp_path):
    stage = _stale_stage(tmp_path, [("old/table.tsv", "partial", False)])

    with run_output.staged_output(tmp_path) as fresh:
        assert not stage.exists()
        assert not (tmp_path / "old/table.tsv").exists()
        (fresh / "new.tsv").write_text("new")

    assert (tmp_path / "new.tsv").read_text() == "new"


def test_sample_preflight_recovers_its_interrupted_partial_output(tmp_path):
    stage = _stale_stage(
        tmp_path,
        [("table_per_sample/WT-1.csv", "partial", False)],
    )

    run_output.assert_sample_output_available(tmp_path, "WT-1")

    assert not stage.exists()
    assert not (tmp_path / "table_per_sample/WT-1.csv").exists()


def test_batch_preflight_recovers_partial_files_and_allows_empty_directories(tmp_path):
    stage = _stale_stage(tmp_path, [("tables/partial.tsv", "partial", False)])

    run_output.assert_output_directory_available(tmp_path)

    assert not stage.exists()
    assert (tmp_path / "tables").is_dir()


def test_batch_preflight_never_reuses_a_directory_containing_a_file(tmp_path):
    (tmp_path / "tables").mkdir()
    (tmp_path / "tables/existing.tsv").write_text("keep")
    with pytest.raises(FileExistsError, match="not empty"):
        run_output.assert_output_directory_available(tmp_path)


def test_recovery_preserves_target_whose_inode_no_longer_matches(tmp_path):
    stage = _stale_stage(tmp_path, [("old.tsv", "partial", False)])
    (tmp_path / "old.tsv").unlink()
    (tmp_path / "old.tsv").write_text("user replacement")

    with run_output.staged_output(tmp_path) as fresh:
        (fresh / "new.tsv").write_text("new")

    assert not stage.exists()
    assert (tmp_path / "old.tsv").read_text() == "user replacement"


def test_recovery_retains_journal_when_owned_target_cannot_be_unlinked(
        tmp_path, monkeypatch):
    stage = _stale_stage(tmp_path, [("old.tsv", "partial", False)])
    original_unlink = run_output._unlink_if_owned
    monkeypatch.setattr(run_output, "_unlink_if_owned", lambda *args: False)

    run_output._recover_stale_transactions(tmp_path)

    assert stage.is_dir()
    assert (tmp_path / "old.tsv").read_text() == "partial"

    monkeypatch.setattr(run_output, "_unlink_if_owned", original_unlink)
    run_output._recover_stale_transactions(tmp_path)
    assert not stage.exists()
    assert not (tmp_path / "old.tsv").exists()


def test_recovery_retains_journal_when_target_cannot_be_inspected(
        tmp_path, monkeypatch):
    stage = _stale_stage(tmp_path, [("old.tsv", "partial", False)])

    def inspection_failure(*args):
        raise PermissionError("temporarily inaccessible")

    monkeypatch.setattr(
        run_output, "_output_entry_identity", inspection_failure,
    )
    run_output._recover_stale_transactions(tmp_path)

    assert stage.is_dir()
    assert (tmp_path / "old.tsv").read_text() == "partial"


def test_nonregular_link_target_is_journaled_and_recovered(
        tmp_path, monkeypatch):
    original_cleanup = run_output._unlink_at_if_owned

    def publish_fifo(source, name, *args, dst_dir_fd=None, **kwargs):
        os.mkfifo(name, dir_fd=dst_dir_fd)

    monkeypatch.setattr(run_output.os, "link", publish_fifo)
    monkeypatch.setattr(
        run_output, "_unlink_at_if_owned", lambda *args: False,
    )

    with pytest.raises(OSError, match="not a regular file"):
        with run_output.staged_output(tmp_path) as stage:
            (stage / "bad.tsv").write_text("data")

    stale = next(tmp_path.glob(".orfbounder-stage-*"))
    journal = json.loads((stale / ".transaction.json").read_text())
    assert journal["targets"][0]["target_identity"] == [
        (tmp_path / "bad.tsv").lstat().st_dev,
        (tmp_path / "bad.tsv").lstat().st_ino,
    ]

    monkeypatch.setattr(run_output, "_unlink_at_if_owned", original_cleanup)
    run_output._recover_stale_transactions(tmp_path)
    assert not (tmp_path / "bad.tsv").exists()
    assert not stale.exists()


def test_recovery_never_follows_an_unsafe_journal_target(tmp_path):
    outside = tmp_path.parent / f"{tmp_path.name}-outside.tsv"
    outside.write_text("keep")
    stage = tmp_path / ".orfbounder-stage-untrusted"
    stage.mkdir()
    run_output.write_json(stage / ".transaction.json", {
        "version": 1,
        "targets": [{
            "relative": f"../{outside.name}",
            "completion": False,
            "target_identity": [outside.stat().st_dev, outside.stat().st_ino],
        }],
    })

    with run_output.staged_output(tmp_path) as fresh:
        (fresh / "new.tsv").write_text("new")

    assert outside.read_text() == "keep"
    assert stage.exists(), "invalid journals are retained for manual inspection"


def test_recovery_retains_an_unknown_journal_version(tmp_path):
    stage = tmp_path / ".orfbounder-stage-newer"
    stage.mkdir()
    (stage / ".transaction.lock").touch()
    run_output.write_json(stage / ".transaction.json", {
        "version": 3,
        "targets": [],
    })

    run_output._recover_stale_transactions(tmp_path)

    assert stage.exists(), "unknown journal semantics must not be guessed"


@pytest.mark.parametrize("control_file", [".transaction.lock", ".transaction.json"])
def test_recovery_does_not_follow_a_linked_control_file(tmp_path, control_file):
    outside = tmp_path / "outside.txt"
    outside_contents = (
        '{"version": 2, "targets": []}'
        if control_file == ".transaction.json" else "keep"
    )
    outside.write_text(outside_contents)
    stage = tmp_path / ".orfbounder-stage-untrusted"
    stage.mkdir()
    (stage / control_file).symlink_to(outside)
    if control_file == ".transaction.lock":
        run_output.write_json(stage / ".transaction.json", {
            "version": 2, "targets": [],
        })
    else:
        (stage / ".transaction.lock").touch()

    run_output._recover_stale_transactions(tmp_path)

    assert stage.is_dir(), "untrusted control files must be retained for inspection"
    assert (stage / control_file).is_symlink()
    assert outside.read_text() == outside_contents


def test_recovery_does_not_block_on_a_fifo_lock(tmp_path):
    stage = tmp_path / ".orfbounder-stage-untrusted"
    stage.mkdir()
    os.mkfifo(stage / ".transaction.lock")
    run_output.write_json(stage / ".transaction.json", {
        "version": 2, "targets": [],
    })

    run_output._recover_stale_transactions(tmp_path)

    assert stage.is_dir()


def test_recovery_preserves_unjournaled_stage_with_unknown_files(tmp_path):
    stage = tmp_path / ".orfbounder-stage-user-data"
    stage.mkdir()
    (stage / ".transaction.lock").touch()
    note = stage / "important.txt"
    note.write_text("keep")

    run_output._recover_stale_transactions(tmp_path)

    assert note.read_text() == "keep"


def test_recovery_does_not_create_a_lock_in_an_unowned_lookalike(tmp_path):
    stage = tmp_path / ".orfbounder-stage-user-data"
    stage.mkdir()
    note = stage / "important.txt"
    note.write_text("keep")

    run_output._recover_stale_transactions(tmp_path)

    assert note.read_text() == "keep"
    assert not (stage / ".transaction.lock").exists()


def test_recovery_removes_a_partially_copied_completion_marker(tmp_path):
    stage = tmp_path / ".orfbounder-stage-interrupted-copy"
    stage.mkdir()
    (stage / ".transaction.lock").touch()
    source = stage / "complete.json"
    source.write_text('{"status":"complete"}\n')
    target = tmp_path / "complete.json"
    target.write_text("partial")
    run_output._write_transaction_journal(stage, [{
        "relative": "complete.json",
        "completion": True,
        "publication_complete": False,
        "target_identity": [target.stat().st_dev, target.stat().st_ino],
    }])

    run_output._recover_stale_transactions(tmp_path)

    assert not target.exists()
    assert not stage.exists()


def test_recovery_never_follows_an_intermediate_target_symlink(tmp_path):
    outside = tmp_path.parent / f"{tmp_path.name}-outside"
    outside.mkdir()
    victim = outside / "victim.tsv"
    victim.write_text("keep")
    (tmp_path / "linked").symlink_to(outside, target_is_directory=True)
    stage = tmp_path / ".orfbounder-stage-untrusted"
    stage.mkdir()
    run_output.write_json(stage / ".transaction.json", {
        "version": 1,
        "targets": [{
            "relative": "linked/victim.tsv",
            "completion": False,
            "target_identity": [victim.stat().st_dev, victim.stat().st_ino],
        }],
    })

    run_output._recover_stale_transactions(tmp_path)

    assert victim.read_text() == "keep"
    assert stage.exists(), "unsafe journals are retained for manual inspection"


def test_recovery_keeps_committed_outputs_and_removes_only_private_stage(tmp_path):
    stage = _stale_stage(tmp_path, [
        ("table.tsv", "results", False), ("complete.json", "complete", True),
    ])

    with run_output.staged_output(tmp_path) as fresh:
        (fresh / "another.tsv").write_text("another")

    assert not stage.exists()
    assert (tmp_path / "table.tsv").read_text() == "results"
    assert (tmp_path / "complete.json").read_text() == "complete"


def test_output_directory_symlink_is_rejected(tmp_path):
    destination = tmp_path / "out"
    destination.mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (destination / "tables").symlink_to(elsewhere, target_is_directory=True)
    with pytest.raises(FileExistsError, match="regular directory"):
        with run_output.staged_output(destination) as stage:
            (stage / "tables").mkdir()
            (stage / "tables/new.tsv").write_text("new")
    assert list(elsewhere.iterdir()) == []


def test_failed_publication_never_rolls_back_through_a_replaced_parent_symlink(
        tmp_path, monkeypatch):
    outside = tmp_path / "outside"
    outside.mkdir()
    destination = tmp_path / "results"
    original_publish = run_output._publish_exclusive
    publications = 0

    def replace_parent_after_first_publish(
            source, root, relative, on_created=None):
        nonlocal publications
        identity = original_publish(source, root, relative, on_created)
        publications += 1
        if publications == 1:
            target = root / relative
            parent = target.parent
            parent.rename(destination / "moved")
            parent.symlink_to(outside, target_is_directory=True)
            os.link(source, outside / target.name)
        return identity

    monkeypatch.setattr(
        run_output, "_publish_exclusive", replace_parent_after_first_publish,
    )

    with pytest.raises(FileExistsError, match="regular directory"):
        with run_output.staged_output(destination) as stage:
            (stage / "tables").mkdir()
            (stage / "tables/a.tsv").write_text("a")
            (stage / "tables/b.tsv").write_text("b")

    assert (outside / "a.tsv").read_text() == "a"
    assert list(destination.glob(".orfbounder-stage-*")), (
        "an incompletely rolled-back transaction must retain its journal"
    )


def test_staging_rejects_destination_symlink_without_touching_target(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    destination = tmp_path / "results"
    destination.symlink_to(elsewhere, target_is_directory=True)

    with pytest.raises(FileExistsError, match="regular directory"):
        with run_output.staged_output(destination) as stage:
            (stage / "result.tsv").write_text("must not be published")

    assert list(elsewhere.iterdir()) == []


@pytest.mark.parametrize(
    "operation",
    [
        lambda path: run_output.assert_sample_output_available(path, "sample"),
        run_output.assert_output_directory_available,
        run_output.recover_stale_transactions,
    ],
)
def test_output_preflights_reject_symlinked_ancestor_without_touching_target(
        tmp_path, operation):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(elsewhere, target_is_directory=True)
    destination = alias / "nested" / "results"

    with pytest.raises(FileExistsError, match="contains a symbolic link"):
        operation(destination)

    assert list(elsewhere.iterdir()) == []


def test_staging_rejects_symlinked_ancestor_without_touching_target(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(elsewhere, target_is_directory=True)
    destination = alias / "nested" / "results"

    with pytest.raises(FileExistsError, match="contains a symbolic link"):
        with run_output.staged_output(destination) as stage:
            (stage / "result.tsv").write_text("must not be published")

    assert list(elsewhere.iterdir()) == []


def test_staging_rejects_link_then_parent_output_alias(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    nested = elsewhere / "nested"
    nested.mkdir(parents=True)
    alias = tmp_path / "alias"
    alias.symlink_to(nested, target_is_directory=True)
    destination = alias / ".." / "results"

    with pytest.raises(ValueError, match="must not contain"):
        with run_output.staged_output(destination) as stage:
            (stage / "result.tsv").write_text("must not be published")

    assert not (elsewhere / "results").exists()
    assert not (tmp_path / "results").exists()


def test_fingerprints_hash_bytes_and_detect_input_changes(tmp_path):
    path = tmp_path / "input.fa"
    data = b">chr\nATGTAA\n"
    path.write_bytes(data)
    cache = {}
    records = run_output.fingerprint_inputs([path, path, None], cache)
    assert len(records) == 1
    assert records[0]["sha256"] == hashlib.sha256(data).hexdigest()
    assert records[0]["size_bytes"] == len(data)
    assert run_output.fingerprint_inputs([path], cache) == records
    run_output.verify_inputs_unchanged(records)
    before = path.stat()
    path.write_bytes(b">chr\nGTGTAA\n")
    # Same-size writes can retain timestamps within a filesystem clock tick.
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
    with pytest.raises(OSError, match="changed during analysis"):
        run_output.verify_inputs_unchanged(records)
    assert run_output.fingerprint_inputs([path], cache)[0]["sha256"] != records[0]["sha256"]


def test_input_snapshots_merge_without_replacing_conflicting_evidence(tmp_path):
    path = tmp_path / "comparison.tsv"
    path.write_text("comparison\n")
    snapshot = run_output.fingerprint_inputs([path])
    merged = run_output.merge_input_records(snapshot, [dict(snapshot[0])])
    assert merged == snapshot

    conflicting = dict(snapshot[0])
    conflicting["sha256"] = "0" * 64
    with pytest.raises(OSError, match="Conflicting input fingerprints"):
        run_output.merge_input_records(snapshot, [conflicting])


def test_retargeted_input_alias_is_detected_even_when_old_target_is_unchanged(tmp_path):
    first, second = tmp_path / "first.bam", tmp_path / "second.bam"
    first.write_bytes(b"first")
    second.write_bytes(b"second")
    alias = tmp_path / "TIS-WT-1.bam"
    alias.symlink_to(first)
    snapshot = run_output.fingerprint_inputs([alias, first])
    assert len(snapshot) == 1
    assert set(snapshot[0]["source_paths"]) == {str(alias), str(first)}
    alias.unlink()
    alias.symlink_to(second)
    with pytest.raises(OSError, match="Input link changed"):
        run_output.verify_inputs_unchanged(snapshot)


def test_final_content_verification_rehashes_and_detects_digest_mismatch(tmp_path):
    source = tmp_path / "input.fa"
    source.write_bytes(b">chr\nATGTAA\n")
    snapshot = run_output.fingerprint_inputs([source])
    snapshot[0]["sha256"] = "0" * 64

    # The ordinary metadata verification still succeeds, while the final
    # content verification independently checks the recorded digest.
    run_output.verify_inputs_unchanged(snapshot)
    with pytest.raises(OSError, match="Input content changed"):
        run_output.verify_inputs_unchanged(snapshot, verify_content=True)


def test_final_content_verification_hashes_duplicate_manifest_inputs_once(tmp_path, monkeypatch):
    source = tmp_path / "input.fa"
    source.write_text("ATGTAA")
    record = run_output.fingerprint_inputs([source])[0]
    calls = 0
    original = run_output._stable_sha256

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(run_output, "_stable_sha256", counted)
    run_output.verify_inputs_unchanged([record, dict(record)], verify_content=True)
    assert calls == 1


def test_existing_track_folder_is_rejected_before_processing(tmp_path):
    (tmp_path / "coverage_files/WT-1").mkdir(parents=True)
    with pytest.raises(FileExistsError, match="already exists"):
        run_output.assert_sample_output_available(tmp_path, "WT-1")
    run_output.assert_sample_output_available(tmp_path, "WT-2")


def test_sample_preflight_rejects_linked_output_subdirectory(tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (tmp_path / "table_per_sample").symlink_to(
        elsewhere, target_is_directory=True,
    )

    with pytest.raises(FileExistsError, match="contains a symbolic link"):
        run_output.assert_sample_output_available(tmp_path, "WT-1")

    assert list(elsewhere.iterdir()) == []


def test_completion_inventory_requires_unchanged_inputs(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text("ATGTAA")
    output = tmp_path / "out"
    output.mkdir()
    manifest = {"output_basename": "WT-1", "input_files": run_output.fingerprint_inputs([source]),
                "status": "analysis_complete"}
    run_output.write_json(output / "WT-1.run.json", manifest)
    (output / "table.tsv").write_text("results")
    final = run_output.complete_run(output, "WT-1", [output / "table.tsv"])
    assert final["status"] == "complete"
    run_output.write_completion(output, "WT-1.complete.json", [final])
    completion = json.loads((output / "WT-1.complete.json").read_text())
    assert completion["samples"] == ["WT-1"]
    assert completion["outputs"] == ["WT-1.run.json", "table.tsv"]
    before = source.stat()
    source.write_text("GTGTAA")
    os.utime(source, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
    with pytest.raises(OSError):
        run_output.write_completion(output, "other.complete.json", [final])
    assert not (output / "other.complete.json").exists()


def test_sequence_exports_have_stable_ids_and_no_terminal_protein_stop(tmp_path):
    frame = pd.DataFrame([{
        "Identifier": "chr:1-9:-", "Locus_tag": "gene\nlabel", "Nucleotide_Seq": "GTGAAATAA",
        "Amino_Acid_Seq": "MK*",
    }])
    outputs = run_output.write_sequences(frame, tmp_path, "WT-1")
    assert outputs[".fna"].read_text() == ">chr:1-9:- gene label\nGTGAAATAA\n"
    assert outputs[".faa"].read_text() == ">chr:1-9:- gene label\nMK\n"
    run_output.write_sequences(frame.iloc[0:0], tmp_path, "empty")
    assert (tmp_path / "empty.fna").read_text() == ""


def test_sequence_exports_omit_missing_annotation_description(tmp_path):
    frame = pd.DataFrame([{
        "Identifier": "chr:1-9:+", "Locus_tag": pd.NA,
        "Nucleotide_Seq": "ATGAAATAA", "Amino_Acid_Seq": "MK*",
    }])

    outputs = run_output.write_sequences(frame, tmp_path, "unannotated")

    assert outputs[".fna"].read_text() == ">chr:1-9:+\nATGAAATAA\n"
    assert outputs[".faa"].read_text() == ">chr:1-9:+\nMK\n"


@pytest.mark.parametrize("column,missing", [
    ("Nucleotide_Seq", pd.NA),
    ("Amino_Acid_Seq", float("nan")),
])
def test_sequence_exports_reject_missing_sequences(tmp_path, column, missing):
    frame = pd.DataFrame([{
        "Identifier": "chr:1-9:+", "Locus_tag": "gene",
        "Nucleotide_Seq": "ATGAAATAA", "Amino_Acid_Seq": "MK*",
    }])
    frame.loc[0, column] = missing

    with pytest.raises(ValueError, match=rf"{column} is missing"):
        run_output.write_sequences(frame, tmp_path, "incomplete")


def test_file_index_escapes_labels_and_link_paths(tmp_path):
    folder = tmp_path / "coverage_files"
    folder.mkdir()
    file = folder / 'sample<&>.wig'
    file.write_text("coverage")
    output = run_output.write_file_index(tmp_path, "files.html", "<unsafe title>", [file])
    html = output.read_text()
    assert "&lt;unsafe title&gt;" in html
    assert "sample%3C%26%3E.wig" in html
    assert "sample&lt;&amp;&gt;.wig" in html
