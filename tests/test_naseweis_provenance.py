"""Acquisition fingerprints identify bytes and reject changes during hashing."""

import hashlib
import os

import pytest

from src import naseweis_provenance as provenance


def test_fingerprint_records_canonical_path_stat_and_sha256(tmp_path):
    path = tmp_path / "acquisition.mzML"
    content = b"mzML acquisition bytes\x00" * 100000
    path.write_bytes(content)
    record = provenance.acquisition_fingerprint(path.parent / "." / path.name)
    assert record == {
        "path": str(path.resolve()),
        "size": len(content),
        "mtime_ns": path.stat().st_mtime_ns,
        "ctime_ns": path.stat().st_ctime_ns,
        "inode": path.stat().st_ino,
        "device": path.stat().st_dev,
        "sha256": hashlib.sha256(content).hexdigest(),
    }
    assert provenance.assert_acquisition_unchanged(path, record) is None


def test_hash_reads_bounded_chunks_without_parsing_mzml(tmp_path, monkeypatch):
    path = tmp_path / "arbitrary.mzML"
    path.write_bytes(b"not XML" * 400000)
    native_hash = hashlib.sha256
    sizes = []

    class RecordingHash:
        def __init__(self):
            self.digest = native_hash()

        def update(self, chunk):
            sizes.append(len(chunk))
            self.digest.update(chunk)

        def hexdigest(self):
            return self.digest.hexdigest()

    monkeypatch.setattr(provenance.hashlib, "sha256", RecordingHash)
    assert (
        provenance.acquisition_fingerprint(path)["sha256"]
        == native_hash(path.read_bytes()).hexdigest()
    )
    assert len(sizes) > 1
    assert max(sizes) <= 1024 * 1024


def test_fingerprint_rejects_a_file_changed_during_hashing(tmp_path, monkeypatch):
    path = tmp_path / "changing.mzML"
    path.write_bytes(b"original")
    native_hash = hashlib.sha256

    class ChangingHash:
        def __init__(self):
            self.digest = native_hash()

        def update(self, chunk):
            self.digest.update(chunk)
            path.write_bytes(b"changed acquisition")

    monkeypatch.setattr(provenance.hashlib, "sha256", ChangingHash)
    with pytest.raises(ValueError, match="changed.*rerun"):
        provenance.acquisition_fingerprint(path)


def test_missing_files_and_directories_are_rejected(tmp_path):
    with pytest.raises(OSError, match="Cannot fingerprint acquisition.*Check the file"):
        provenance.acquisition_fingerprint(tmp_path / "missing.mzML")
    with pytest.raises(ValueError, match="not a regular file"):
        provenance.acquisition_fingerprint(tmp_path)


@pytest.mark.parametrize("change", ["bytes", "mtime", "missing", "path"])
def test_post_search_check_rejects_changed_acquisition(tmp_path, change):
    path = tmp_path / "original.mzML"
    path.write_bytes(b"acquisition")
    record = provenance.acquisition_fingerprint(path)
    if change == "bytes":
        path.write_bytes(b"replacement bytes")
    elif change == "mtime":
        os.utime(path, ns=(record["mtime_ns"] + 1000000, record["mtime_ns"] + 1000000))
    elif change == "missing":
        path.unlink()
    else:
        copy = tmp_path / "copy.mzML"
        copy.write_bytes(path.read_bytes())
        path = copy
    with pytest.raises(ValueError, match="changed during the search.*rerun"):
        provenance.assert_acquisition_unchanged(path, record)


@pytest.mark.parametrize("record", [None, [], {}, {"path": "wrong"}])
def test_post_search_check_rejects_missing_or_invalid_records(tmp_path, record):
    path = tmp_path / "acquisition.mzML"
    path.touch()
    with pytest.raises(ValueError, match="rerun"):
        provenance.assert_acquisition_unchanged(path, record)


def test_identical_restored_bytes_keep_the_same_content_identity(tmp_path):
    path = tmp_path / "acquisition.mzML"
    path.write_bytes(b"original")
    before = provenance.acquisition_fingerprint(path)
    path.write_bytes(b"replacement")
    path.write_bytes(b"original")
    os.utime(path, ns=(before["mtime_ns"] + 1000000, before["mtime_ns"] + 1000000))
    restored = provenance.acquisition_fingerprint(path)
    assert (restored["sha256"], restored["size"]) == (before["sha256"], before["size"])
    assert restored["mtime_ns"] != before["mtime_ns"]


@pytest.mark.parametrize("replacement", ["overwrite", "replace"])
def test_same_size_replacement_with_restored_mtime_is_rejected(tmp_path, replacement):
    path = tmp_path / "acquisition.mzML"
    path.write_bytes(b"original")
    before = provenance.acquisition_fingerprint(path)
    if replacement == "overwrite":
        path.write_bytes(b"replaced")
    else:
        other = tmp_path / "replacement.mzML"
        other.write_bytes(b"replaced")
        other.replace(path)
    os.utime(path, ns=(before["mtime_ns"], before["mtime_ns"]))
    assert path.stat().st_size == before["size"]
    assert path.stat().st_mtime_ns == before["mtime_ns"]
    with pytest.raises(ValueError, match="changed during the search"):
        provenance.assert_acquisition_unchanged(path, before)


def test_windows_post_search_check_hashes_even_when_all_stat_fields_match(
    tmp_path, monkeypatch
):
    path = tmp_path / "acquisition.mzML"
    path.write_bytes(b"original")
    original = provenance.acquisition_fingerprint(path)
    path.write_bytes(b"replaced")
    os.utime(path, ns=(original["mtime_ns"], original["mtime_ns"]))
    # Model creation-time ctime semantics: all current stat fields can match the
    # captured record even though bytes changed without changing size or mtime.
    record = provenance.acquisition_fingerprint(path) | {"sha256": original["sha256"]}
    monkeypatch.setattr(provenance, "requires_content_verification", lambda: True)
    with pytest.raises(ValueError, match="changed during the search"):
        provenance.assert_acquisition_unchanged(path, record)


def test_windows_unchanged_acquisition_passes_fresh_content_verification(
    tmp_path, monkeypatch
):
    path = tmp_path / "acquisition.mzML"
    path.write_bytes(b"original")
    record = provenance.acquisition_fingerprint(path)
    hasher = provenance.acquisition_fingerprint
    calls = []

    def fingerprint(source):
        calls.append(source)
        return hasher(source)

    monkeypatch.setattr(provenance, "requires_content_verification", lambda: True)
    monkeypatch.setattr(provenance, "acquisition_fingerprint", fingerprint)
    assert provenance.assert_acquisition_unchanged(path, record) is None
    assert len(calls) == 1
