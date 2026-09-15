"""Input provenance and safe publication of a complete set of analysis outputs."""

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import errno
import json
import os
from pathlib import Path
import re
import shutil
import stat
import tempfile
from urllib.parse import quote
from html import escape

import pandas as pd

try:  # POSIX advisory locks are available in the supported container/Linux path.
    import fcntl
except ImportError:  # pragma: no cover - conservative fallback on non-POSIX hosts
    fcntl = None


_STAGE_PREFIX = ".orfbounder-stage-"
_LOCK_NAME = ".transaction.lock"
_JOURNAL_NAME = ".transaction.json"


def _is_transaction_file(path):
    return Path(path).name in {_LOCK_NAME, _JOURNAL_NAME, _JOURNAL_NAME + ".next"}


def _identity(stat):
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def _stable_sha256(path, expected_identity):
    """Hash a file while proving that its filesystem identity stayed stable."""
    path = Path(path)
    if _identity(path.stat()) != expected_identity:
        raise OSError(f"Input changed during analysis: {path}. Run again with stable input files.")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        if _identity(os.fstat(handle.fileno())) != expected_identity:
            raise OSError(f"Input changed while opening it: {path}. Run again with stable input files.")
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
        if _identity(os.fstat(handle.fileno())) != expected_identity:
            raise OSError(f"Input changed while hashing it: {path}. Run again with stable input files.")
    if _identity(path.stat()) != expected_identity:
        raise OSError(f"Input changed while hashing it: {path}. Run again with stable input files.")
    return digest.hexdigest()


def fingerprint_inputs(paths, cache=None):
    """Hash each distinct input once, reusing only unchanged file snapshots."""
    cache = cache if cache is not None else {}
    records, aliases = [], {}
    for source in paths:
        if source is not None:
            source = Path(source).absolute()
            aliases.setdefault(source.resolve(strict=True), set()).add(str(source))
    for path in sorted(aliases):
        before = path.stat()
        identity = _identity(before)
        cached = cache.get(str(path))
        if cached is None or cached[0] != identity:
            digest = _stable_sha256(path, identity)
            record = {"path": str(path), "size_bytes": before.st_size, "sha256": digest,
                      "mtime_ns": before.st_mtime_ns, "ctime_ns": before.st_ctime_ns,
                      "device": before.st_dev, "inode": before.st_ino}
            cache[str(path)] = (identity, record)
        records.append({**cache[str(path)][1], "source_paths": sorted(aliases[path])})
    verify_inputs_unchanged(records)
    return records


def verify_inputs_unchanged(records, *, verify_content=False):
    """Verify aliases and metadata, optionally rehashing each canonical file once."""
    verified = {}
    for record in records:
        path = Path(record["path"])
        for source in record.get("source_paths", [str(path)]):
            if Path(source).resolve(strict=True) != path:
                raise OSError(f"Input link changed during analysis: {source}. Run again with stable input files.")
        stat = path.stat()
        expected = (record["device"], record["inode"], record["size_bytes"],
                    record["mtime_ns"], record["ctime_ns"])
        if _identity(stat) != expected:
            raise OSError(f"Input changed during analysis: {path}. Run again with stable input files.")
        previous = verified.get(str(path))
        if previous is not None and previous != (expected, record["sha256"]):
            raise OSError(f"Conflicting input fingerprints recorded for {path}.")
        if previous is None and verify_content:
            actual = _stable_sha256(path, expected)
            if actual != record["sha256"]:
                raise OSError(f"Input content changed during analysis: {path}. Run again with stable input files.")
        verified[str(path)] = (expected, record["sha256"])


def merge_input_records(*record_groups):
    """Combine verified fingerprint snapshots without replacing older evidence."""
    records = [record for group in record_groups for record in group]
    verify_inputs_unchanged(records)
    merged = {}
    for record in records:
        path = record["path"]
        if path not in merged:
            merged[path] = {**record, "source_paths": list(record.get("source_paths", [path]))}
            continue
        aliases = set(merged[path].get("source_paths", [path]))
        aliases.update(record.get("source_paths", [path]))
        merged[path]["source_paths"] = sorted(aliases)
    return [merged[path] for path in sorted(merged)]


def _absolute_output_path(path):
    """Return a lexical absolute path without resolving symbolic links."""
    if ".." in Path(path).parts:
        # abspath would collapse ``linked/..`` before a no-follow inspection,
        # while the kernel traverses the link before processing ``..``.
        raise ValueError(f"Output path must not contain '..' components: {path}")
    return Path(os.path.abspath(os.fspath(path)))


def _assert_no_output_symlink(path):
    """Reject a pre-existing symbolic link anywhere in an output path.

    Output aliases are ambiguous publication boundaries: a path that appears
    to be below the requested result directory can otherwise write through a
    linked ancestor. Inputs intentionally support aliases and record them in
    their fingerprints, but output paths use a stricter no-link contract.
    """
    absolute = _absolute_output_path(path)
    current = Path(absolute.anchor)
    for component in absolute.parts[1:]:
        current /= component
        try:
            current_status = current.lstat()
        except FileNotFoundError:
            # No later component can exist before this one is created.
            return
        if stat.S_ISLNK(current_status.st_mode):
            raise FileExistsError(
                f"Output path is not a regular directory; contains a symbolic link: {current}"
            )


@contextmanager
def _open_directory_path(path):
    """Open an absolute directory path one component at a time, without links."""
    required = ("O_DIRECTORY", "O_NOFOLLOW")
    if any(not hasattr(os, name) for name in required):
        raise OSError(errno.ENOTSUP, "No-follow directory operations are unavailable")
    flags = (getattr(os, "O_PATH", os.O_RDONLY) | os.O_DIRECTORY
             | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0))
    absolute = _absolute_output_path(path)
    descriptors = []
    try:
        descriptor = os.open(absolute.anchor, flags)
        descriptors.append(descriptor)
        for component in absolute.parts[1:]:
            descriptor = os.open(component, flags, dir_fd=descriptor)
            descriptors.append(descriptor)
        yield descriptors[-1]
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


@contextmanager
def _open_parent_directory(root, relative):
    """Open a relative parent one component at a time without following links."""
    relative = Path(relative)
    if (relative.is_absolute() or not relative.parts or ".." in relative.parts
            or relative.name in {"", "."}):
        raise ValueError(f"Unsafe relative output path: {relative}")
    flags = (getattr(os, "O_PATH", os.O_RDONLY) | os.O_DIRECTORY
             | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0))
    descriptors = []
    with _open_directory_path(root) as root_descriptor:
        try:
            descriptor = root_descriptor
            for component in relative.parts[:-1]:
                descriptor = os.open(component, flags, dir_fd=descriptor)
                descriptors.append(descriptor)
            yield descriptor, relative.name
        finally:
            for descriptor in reversed(descriptors):
                os.close(descriptor)


def _parent_chain_is_safe(root, relative):
    try:
        with _open_parent_directory(root, relative):
            return True
    except FileNotFoundError:
        # A missing parent means that no target currently exists to follow.
        return True
    except (OSError, ValueError):
        return False


def _regular_file_identity(root, relative):
    try:
        with _open_parent_directory(root, relative) as (parent, name):
            current = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(current.st_mode):
        return None
    return current.st_dev, current.st_ino


def _output_entry_identity(root, relative):
    """Return the identity of a non-directory entry without following links."""
    try:
        with _open_parent_directory(root, relative) as (parent, name):
            current = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except FileNotFoundError:
        return None
    if stat.S_ISDIR(current.st_mode):
        return None
    return current.st_dev, current.st_ino


def _unlink_at_if_owned(parent, name, expected_identity):
    """Unlink a non-directory entry at a stable parent iff it is still ours."""
    try:
        current = os.stat(name, dir_fd=parent, follow_symlinks=False)
        if (stat.S_ISDIR(current.st_mode)
                or (current.st_dev, current.st_ino) != tuple(expected_identity)):
            return False
        os.unlink(name, dir_fd=parent)
        return True
    except FileNotFoundError:
        return True
    except OSError:
        return False


def _unlink_if_owned(root, relative, expected_identity):
    """Unlink through a stable parent descriptor only when identity still matches."""
    try:
        with _open_parent_directory(root, relative) as (parent, name):
            return _unlink_at_if_owned(parent, name, expected_identity)
    except FileNotFoundError:
        return True
    except (OSError, ValueError):
        return False


def _remove_empty_directory(root, relative):
    try:
        with _open_parent_directory(root, relative) as (parent, name):
            os.rmdir(name, dir_fd=parent)
        return True
    except FileNotFoundError:
        return True
    except (OSError, ValueError):
        return False


def _remove_empty_directory_if_owned(root, relative, expected_identity):
    """Remove a newly created empty directory only while its identity matches."""
    try:
        with _open_parent_directory(root, relative) as (parent, name):
            current = os.stat(name, dir_fd=parent, follow_symlinks=False)
            if (not stat.S_ISDIR(current.st_mode)
                    or (current.st_dev, current.st_ino) != expected_identity):
                return False
            os.rmdir(name, dir_fd=parent)
        return True
    except (OSError, ValueError):
        return False


def _ensure_parent_directories(root, relative):
    """Create missing output parents without traversing a substituted link."""
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ValueError(f"Unsafe relative output path: {relative}")
    flags = (getattr(os, "O_PATH", os.O_RDONLY) | os.O_DIRECTORY
             | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0))
    descriptors, created, prefix = [], [], []
    with _open_directory_path(root) as root_descriptor:
        try:
            descriptor = root_descriptor
            for component in relative.parts[:-1]:
                prefix.append(component)
                try:
                    child = os.open(component, flags, dir_fd=descriptor)
                except FileNotFoundError:
                    try:
                        os.mkdir(component, dir_fd=descriptor)
                        created.append(Path(*prefix).as_posix())
                    except FileExistsError:
                        pass
                    child = os.open(component, flags, dir_fd=descriptor)
                descriptors.append(child)
                descriptor = child
            return created
        finally:
            for descriptor in reversed(descriptors):
                os.close(descriptor)


def _publish_exclusive(source, destination, relative, on_created=None):
    """Publish exclusively and report ownership as soon as a target exists."""
    source = Path(source)
    with _open_parent_directory(destination, relative) as (parent, name):
        source_status = os.stat(source, follow_symlinks=False)
        if not stat.S_ISREG(source_status.st_mode):
            raise OSError(f"Analysis output is not a regular file: {relative}")
        published_identity = None
        try:
            try:
                os.link(source, name, dst_dir_fd=parent, follow_symlinks=False)
            except OSError as exc:
                unsupported = {errno.EXDEV, errno.EPERM, errno.EOPNOTSUPP}
                if hasattr(errno, "ENOTSUP"):
                    unsupported.add(errno.ENOTSUP)
                if exc.errno not in unsupported:
                    raise
                source_descriptor = os.open(
                    source,
                    os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
                    | getattr(os, "O_CLOEXEC", 0)
                    | getattr(os, "O_NONBLOCK", 0),
                )
                descriptor = -1
                try:
                    opened_source = os.fstat(source_descriptor)
                    if not stat.S_ISREG(opened_source.st_mode):
                        raise OSError(
                            f"Analysis output is not a regular file: {relative}"
                        )
                    descriptor = os.open(
                        name, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                        opened_source.st_mode & 0o777, dir_fd=parent,
                    )
                    published = os.fstat(descriptor)
                    published_identity = (published.st_dev, published.st_ino)
                    if on_created is not None:
                        on_created(published_identity)
                    with os.fdopen(source_descriptor, "rb") as source_handle:
                        source_descriptor = -1
                        with os.fdopen(descriptor, "wb") as target_handle:
                            descriptor = -1
                            shutil.copyfileobj(
                                source_handle, target_handle,
                                length=1024 * 1024,
                            )
                            target_handle.flush()
                            os.fsync(target_handle.fileno())
                finally:
                    if source_descriptor >= 0:
                        os.close(source_descriptor)
                    if descriptor >= 0:
                        os.close(descriptor)
            else:
                published = os.stat(name, dir_fd=parent, follow_symlinks=False)
                published_identity = (published.st_dev, published.st_ino)
                if on_created is not None:
                    on_created(published_identity)
                if not stat.S_ISREG(published.st_mode):
                    raise OSError(
                        f"Published output is not a regular file: {relative}"
                    )
            return published_identity
        except BaseException:
            if published_identity is not None:
                _unlink_at_if_owned(parent, name, published_identity)
            raise


def _write_transaction_journal(stage, entries):
    temporary = stage / (_JOURNAL_NAME + ".next")
    write_json(temporary, {"version": 2, "targets": entries})
    os.replace(temporary, stage / _JOURNAL_NAME)


def _validate_transaction_entries(raw_entries):
    """Reject journals that could name files outside their own destination."""
    if not isinstance(raw_entries, list):
        raise ValueError("Transaction targets must be a list.")
    entries, seen = [], set()
    for entry in raw_entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("relative"), str):
            raise ValueError("Invalid transaction target entry.")
        relative = Path(entry["relative"])
        if (not entry["relative"] or not relative.parts or relative.name in {"", "."}
                or relative.is_absolute() or ".." in relative.parts
                or (relative.parts and relative.parts[0].startswith(_STAGE_PREFIX))):
            raise ValueError("Unsafe transaction target path.")
        normalized = relative.as_posix()
        if normalized in seen:
            raise ValueError("Duplicate transaction target path.")
        seen.add(normalized)
        if not isinstance(entry.get("completion"), bool):
            raise ValueError("Invalid transaction completion marker.")
        publication_complete = entry.get("publication_complete")
        if (publication_complete is not None
                and not isinstance(publication_complete, bool)):
            raise ValueError("Invalid publication completion state.")
        identity = entry.get("target_identity")
        if identity is not None and (
                not isinstance(identity, list) or len(identity) != 2
                or any(isinstance(value, bool) or not isinstance(value, int) for value in identity)):
            raise ValueError("Invalid transaction target identity.")
        entries.append({**entry, "relative": normalized})
    return entries


def _open_regular_stage_file(stage, name, flags, mode=0o600):
    """Open a transaction control file through a stable, no-follow parent."""
    with _open_directory_path(stage) as directory:
        descriptor = os.open(
            name,
            flags | os.O_NOFOLLOW | getattr(os, "O_NONBLOCK", 0)
            | getattr(os, "O_CLOEXEC", 0),
            mode,
            dir_fd=directory,
        )
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise OSError(errno.EINVAL, f"Transaction control file is not regular: {stage / name}")
    return descriptor


def _recover_stale_transactions(destination):
    """Roll back files provably owned by interrupted publication attempts.

    Active stages remain locked and are never touched. A completed transaction
    is recognized by its completion marker and only has its private stage
    removed. For an incomplete transaction, targets are removed only when their
    inode matches the staged source (hard link) or the recorded copied target.
    """
    if fcntl is None:
        return
    for stage in sorted(destination.glob(_STAGE_PREFIX + "*")):
        if stage.is_symlink() or not stage.is_dir():
            continue
        try:
            lock = os.fdopen(
                _open_regular_stage_file(stage, _LOCK_NAME, os.O_RDWR),
                "r+b",
            )
        except OSError:
            continue
        try:
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except (BlockingIOError, OSError):
                continue
            try:
                journal_descriptor = _open_regular_stage_file(
                    stage, _JOURNAL_NAME, os.O_RDONLY,
                )
            except FileNotFoundError:
                # Without a journal no target ownership can be proved. Only
                # an otherwise empty stage can be discarded automatically;
                # preserve any staged files or user-created lookalike for
                # inspection instead of deleting unknown data.
                if not any(child.name != _LOCK_NAME for child in stage.iterdir()):
                    shutil.rmtree(stage)
                continue
            except OSError:
                # A linked or non-regular journal cannot establish ownership.
                continue
            try:
                with os.fdopen(journal_descriptor, "r", encoding="utf-8") as handle:
                    journal = json.load(handle)
                journal_version = journal.get("version")
                if journal_version not in {1, 2}:
                    raise ValueError("Unsupported transaction journal version.")
                entries = _validate_transaction_entries(journal["targets"])
            except (OSError, ValueError, KeyError, TypeError):
                # A corrupt journal cannot establish ownership, so preserve it
                # for manual inspection instead of guessing what to delete.
                continue

            if any(
                    not _parent_chain_is_safe(destination, entry["relative"])
                    or not _parent_chain_is_safe(stage, entry["relative"])
                    for entry in entries):
                # Never resolve a journal target through a link or another
                # non-directory entry. Preserve the stage for inspection just
                # like a syntactically unsafe journal.
                continue

            def owned_target_identity(entry, *, require_regular=False):
                copied = entry.get("target_identity")
                identity_reader = (
                    _regular_file_identity if require_regular
                    else _output_entry_identity
                )
                target = identity_reader(destination, entry["relative"])
                if copied and target == tuple(copied):
                    return target
                source = _regular_file_identity(stage, entry["relative"])
                return source if target is not None and target == source else None

            def publication_is_complete(entry):
                target = _regular_file_identity(
                    destination, entry["relative"],
                )
                if target is None:
                    return False
                source = _regular_file_identity(stage, entry["relative"])
                if source is not None and target == source:
                    # A hard link appears atomically with its complete inode.
                    return True
                if journal_version == 1:
                    # Version 1 recorded copied identities only after the copy
                    # and fsync completed.
                    return target == tuple(entry.get("target_identity") or ())
                return (
                    entry.get("publication_complete") is True
                    and target == tuple(entry.get("target_identity") or ())
                )

            rollback_complete = True
            try:
                completion_entries = [
                    entry for entry in entries if entry.get("completion")
                ]
                committed = bool(completion_entries) and all(
                    owned_target_identity(entry, require_regular=True)
                    is not None
                    and publication_is_complete(entry)
                    for entry in completion_entries
                )
                if not committed:
                    for entry in reversed(entries):
                        identity = owned_target_identity(entry)
                        if identity is not None and not _unlink_if_owned(
                                destination, entry["relative"], identity):
                            rollback_complete = False
            except (OSError, ValueError):
                # A path that became uninspectable cannot be classified as
                # absent, replaced, or owned. Keep its journal as evidence.
                rollback_complete = False
            if rollback_complete:
                shutil.rmtree(stage)
        finally:
            lock.close()


def assert_sample_output_available(destination, basename):
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", basename) is None:
        raise ValueError("output_basename must be a filename without path separators, starting with a letter or digit.")
    destination = Path(destination)
    _assert_no_output_symlink(destination)
    if destination.is_symlink() or (destination.exists() and not destination.is_dir()):
        raise FileExistsError(f"Output path is not a regular directory: {destination}")
    if destination.exists():
        _recover_stale_transactions(destination)
    expected = [destination / f"{basename}{suffix}" for suffix in
                (".run.json", ".complete.json", ".report.html", ".files.html")]
    expected += [destination / directory / f"{basename}{suffix}" for directory, suffix in
                 (("table_per_sample", ".csv"), ("table_per_sample", ".xlsx"),
                  ("gff_per_sample", ".gff"), ("sequence_per_sample", ".fna"),
                  ("sequence_per_sample", ".faa"), ("coverage_files", ""))]
    for path in expected:
        _assert_no_output_symlink(path)
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Output already exists: {path}. Choose a new output path or basename.")


def recover_stale_transactions(destination):
    """Safely recover interrupted publications without requiring an empty root."""
    destination = Path(destination)
    _assert_no_output_symlink(destination)
    if destination.is_symlink() or (
            destination.exists() and not destination.is_dir()):
        raise FileExistsError(
            f"Output path is not a regular directory: {destination}"
        )
    if destination.exists():
        _recover_stale_transactions(destination)


def assert_output_directory_available(destination):
    """Recover interrupted work, then require no existing result files.

    Empty directories are safe to reuse because publication remains exclusive
    at the file level. Files, links, and other entries are never overwritten.
    """
    destination = Path(destination)
    _assert_no_output_symlink(destination)
    if destination.is_symlink() or (destination.exists() and not destination.is_dir()):
        raise FileExistsError(f"Result folder is not a regular directory: {destination}")
    if not destination.exists():
        return
    _recover_stale_transactions(destination)
    occupied = next(
        (path for path in destination.rglob("*")
         if path.is_symlink() or not path.is_dir()),
        None,
    )
    if occupied is not None:
        raise FileExistsError(
            f"Result folder is not empty: {destination}. "
            "Choose a new result path to preserve existing results."
        )


def write_json(path, contents):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(contents, indent=2, allow_nan=False) + "\n", encoding="utf-8")


@contextmanager
def staged_output(destination):
    """Prepare files privately, then publish without replacing existing files.

    A completion record is published last. If normal publication fails, remove
    only files created by this transaction. A later run recovers an interrupted
    transaction when ownership can be proved from its journal; consumers must
    still require the completion file. Staging lives on the output filesystem
    so exclusive hard links are atomic where supported. A newly created output
    directory is removed on failure only if it is still the same empty directory.
    """
    destination = Path(destination)
    _assert_no_output_symlink(destination)
    if destination.is_symlink() or (destination.exists() and not destination.is_dir()):
        raise FileExistsError(f"Output path is not a regular directory: {destination}")
    created_destination_identity = None
    try:
        destination.mkdir(parents=True)
    except FileExistsError:
        # A pre-existing directory belongs to its owner, even if this run
        # later fails and leaves it empty.
        pass
    else:
        created = destination.lstat()
        created_destination_identity = (created.st_dev, created.st_ino)
    _assert_no_output_symlink(destination)
    # Repeat the check after creation so a root replaced between the initial
    # validation and mkdir cannot redirect staging through a symbolic link.
    if destination.is_symlink() or not destination.is_dir():
        raise FileExistsError(f"Output path is not a regular directory: {destination}")
    _recover_stale_transactions(destination)
    staging = Path(tempfile.mkdtemp(prefix=_STAGE_PREFIX, dir=destination))
    lock = (staging / _LOCK_NAME).open("a+b")
    if fcntl is not None:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
    published, created_dirs = [], []
    remove_staging = True
    failed = False
    try:
        yield staging
        sources = sorted(path for path in staging.rglob("*") if path.is_file()
                         and not _is_transaction_file(path))
        if any(path.is_symlink() for path in staging.rglob("*")):
            raise ValueError("Analysis outputs must not contain symbolic links.")
        # Completion records are only visible after all other files exist.
        sources.sort(key=lambda path: path.name == "complete.json" or path.name.endswith(".complete.json"))
        journal = [{
            "relative": source.relative_to(staging).as_posix(),
            "completion": source.name == "complete.json" or source.name.endswith(".complete.json"),
            "publication_complete": False,
        } for source in sources]
        _write_transaction_journal(staging, journal)
        for index, source in enumerate(sources):
            target = destination / source.relative_to(staging)
            if target.exists() or target.is_symlink():
                raise FileExistsError(f"Output already exists: {target}. Choose a new output path or basename.")
            parent = target.parent
            while parent != destination:
                if parent.is_symlink() or (parent.exists() and not parent.is_dir()):
                    raise FileExistsError(f"Output folder is not a regular directory: {parent}")
                parent = parent.parent
        for index, source in enumerate(sources):
            relative = source.relative_to(staging)
            relative_text = relative.as_posix()

            def record_publication(identity):
                device, inode = identity
                published.append((relative_text, device, inode))
                journal[index]["target_identity"] = [device, inode]
                _write_transaction_journal(staging, journal)

            try:
                created_dirs.extend(
                    _ensure_parent_directories(destination, relative)
                )
                _publish_exclusive(
                    source, destination, relative, record_publication,
                )
                journal[index]["publication_complete"] = True
                _write_transaction_journal(staging, journal)
            except OSError as exc:
                if exc.errno in {errno.ELOOP, errno.ENOTDIR}:
                    raise FileExistsError(
                        f"Output folder is not a regular directory: "
                        f"{(destination / relative).parent}"
                    ) from exc
                raise
    except BaseException:
        failed = True
        rollback_complete = True
        for relative, device, inode in reversed(published):
            if not _unlink_if_owned(destination, relative, (device, inode)):
                rollback_complete = False
        for relative in reversed(created_dirs):
            if not _remove_empty_directory(destination, relative):
                rollback_complete = False
        # If anything could not be safely rolled back, retain the locked stage
        # and journal so a later run or a human can establish ownership.
        remove_staging = rollback_complete
        raise
    finally:
        lock.close()
        if remove_staging:
            shutil.rmtree(staging, ignore_errors=True)
            if failed and created_destination_identity is not None:
                _remove_empty_directory_if_owned(
                    destination.parent, destination.name,
                    created_destination_identity,
                )


def write_sequences(frame, directory, basename):
    """Export coordinate-identified nucleotide and translated protein FASTA."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {suffix: directory / f"{basename}{suffix}" for suffix in (".fna", ".faa")}

    def write_record(handle, identifier, description, sequence):
        header = f">{identifier}"
        if description:
            header += f" {description}"
        handle.write(header + "\n")
        for start in range(0, len(sequence), 60):
            handle.write(sequence[start:start + 60] + "\n")

    def required_sequence(value, column, identifier):
        if pd.isna(value) or not isinstance(value, str) or not value:
            raise ValueError(
                f"Cannot export FASTA record {identifier!r}: {column} is missing."
            )
        return value

    columns = ["Identifier", "Locus_tag", "Nucleotide_Seq", "Amino_Acid_Seq"]
    with outputs[".fna"].open("w", encoding="utf-8") as nucleotide, \
            outputs[".faa"].open("w", encoding="utf-8") as protein:
        for identifier, description, nt_sequence, aa_sequence in frame[columns].itertuples(index=False, name=None):
            # Identifiers are coordinate-based and contain no whitespace.
            identifier = "_".join(str(identifier).split())
            description = "" if pd.isna(description) else " ".join(str(description).split())
            nt_sequence = required_sequence(nt_sequence, "Nucleotide_Seq", identifier)
            aa_sequence = required_sequence(aa_sequence, "Amino_Acid_Seq", identifier)
            write_record(nucleotide, identifier, description, nt_sequence)
            write_record(protein, identifier, description, aa_sequence.removesuffix("*"))
    return outputs


def export_sample_results(frame, directory, basename, split_gff=False):
    from lib import io
    directory = Path(directory)
    io.write_results_to_gff(frame, directory / "gff_per_sample", basename, split_gff)
    io.write_results_to_table(frame, directory / "table_per_sample", basename)
    write_sequences(frame, directory / "sequence_per_sample", basename)


def report_links(directory, basename, merged=False):
    """Relative URLs remain valid after staging is published or results moved."""
    paths = {
        "Spreadsheet": f"{basename}.xlsx" if merged else f"table_per_sample/{basename}.xlsx",
        "Tab-separated table": f"{basename}.csv" if merged else f"table_per_sample/{basename}.csv",
        "Genome annotation": f"{basename}.gff" if merged else f"gff_per_sample/{basename}.gff",
        "Nucleotide sequences": f"{basename}.fna" if merged else f"sequence_per_sample/{basename}.fna",
        "Protein sequences": f"{basename}.faa" if merged else f"sequence_per_sample/{basename}.faa",
    }
    return {label: quote(path, safe="/") for label, path in paths.items() if (Path(directory) / path).is_file()}


def write_file_index(directory, name, title, paths):
    """List downloadable tracks/statistics for an offline browser report."""
    directory = Path(directory)
    items = []
    for path in sorted(set(map(Path, paths))):
        relative = path.relative_to(directory).as_posix()
        items.append(f'<li><a href="{quote(relative, safe="/")}">{escape(relative)}</a></li>')
    output = directory / name
    output.write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">'
        f'<title>{escape(title)}</title><style>body{{font:17px system-ui;max-width:65rem;margin:3rem auto;padding:0 1rem;line-height:1.6}}'
        'a{color:#165b82}li{overflow-wrap:anywhere}</style>'
        f'<h1>{escape(title)}</h1><p>Open coverage tracks with the matching reference genome in your genome browser. '
        'Statistical tables are tab-separated files you can open in a spreadsheet.</p>'
        + ("<ul>" + "".join(items) + "</ul>" if items else "<p>No coverage tracks or statistical tables were produced.</p>")
        + "</html>", encoding="utf-8")
    return output


def complete_run(directory, basename, output_files):
    """Mark a sample manifest complete after its requested files are written."""
    path = Path(directory) / f"{basename}.run.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    verify_inputs_unchanged(record["input_files"])
    record["status"] = "complete"
    record["completed_at"] = datetime.now(timezone.utc).isoformat()
    record["outputs"] = sorted(
        str(Path(file).relative_to(directory)) for file in output_files
        if not _is_transaction_file(file)
    )
    write_json(path, record)
    return record


def write_completion(directory, filename, manifests):
    """Write a final inventory after all exports; publication puts it last."""
    directory = Path(directory)
    # Re-read every distinct input at the final commit point. Metadata checks
    # during analysis catch ordinary edits cheaply; this closes the same-size,
    # timestamp-preserving gap before the completion marker is created.
    verify_inputs_unchanged(
        [record for manifest in manifests for record in manifest["input_files"]],
        verify_content=True,
    )
    outputs = sorted(path.relative_to(directory).as_posix() for path in directory.rglob("*")
                     if path.is_file() and not _is_transaction_file(path))
    write_json(directory / filename, {"status": "complete", "completed_at": datetime.now(timezone.utc).isoformat(),
                                      "samples": [record["output_basename"] for record in manifests],
                                      "outputs": outputs})
