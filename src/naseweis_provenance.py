"""Bind a saved search to the acquisition bytes that were actually searched."""

import hashlib
import os
from pathlib import Path


def _identity(stat):
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def requires_content_verification():
    """Windows ctime may denote creation, so stat identity alone is insufficient."""
    return os.name == "nt"


def _record_identity(path, stat):
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "inode": stat.st_ino,
        "device": stat.st_dev,
    }


def acquisition_fingerprint(path):
    """Hash one acquisition in bounded chunks, rejecting concurrent changes."""
    try:
        source = Path(path).resolve(strict=True)
        if not source.is_file():
            raise ValueError(f"Acquisition is not a regular file: {source}")
        before = source.stat()
        digest = hashlib.sha256()
        with source.open("rb") as acquisition:
            if _identity(before) != _identity(os.fstat(acquisition.fileno())):
                raise ValueError(
                    "The acquisition changed while being fingerprinted. Select the file again and rerun the search."
                )
            while chunk := acquisition.read(1024 * 1024):
                digest.update(chunk)
            after_read = os.fstat(acquisition.fileno())
        after = source.stat()
    except OSError as exc:
        raise OSError(
            f"Cannot fingerprint acquisition {path}: {exc}. Check the file and rerun the search."
        ) from exc
    if _identity(before) != _identity(after_read) or _identity(before) != _identity(
        after
    ):
        raise ValueError(
            "The acquisition changed while being fingerprinted. Select the file again and rerun the search."
        )
    return _record_identity(source, after) | {"sha256": digest.hexdigest()}


def assert_acquisition_unchanged(path, record):
    """Reject a file changed or replaced while a search was running."""
    message = "The acquisition changed during the search or is no longer available. Select the original file and rerun the search before viewing spectra or running quantification."
    try:
        source = Path(path).resolve(strict=True)
        info = source.stat()
        unchanged = (
            source.is_file()
            and isinstance(record, dict)
            and all(
                record.get(key) == value
                for key, value in _record_identity(source, info).items()
            )
        )
        if unchanged and requires_content_verification():
            current = acquisition_fingerprint(source)
            unchanged = all(record.get(key) == value for key, value in current.items())
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError(message) from exc
    if not unchanged:
        raise ValueError(message)
