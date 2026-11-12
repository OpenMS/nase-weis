"""Check live local workflows without mistaking orphaned PID files for workers."""

from pathlib import Path

import psutil


def pid_markers(pid_dir: Path) -> list[Path]:
    """Snapshot positive numeric regular markers; a removed directory is idle.

    Never treat directories, symlinks, process-group ID zero or unrelated files
    as owned processes. Other I/O errors remain visible to the caller.
    """
    try:
        markers = list(Path(pid_dir).iterdir())
    except FileNotFoundError:
        return []
    return [
        marker for marker in markers
        if marker.name.isascii() and marker.name.isdecimal()
        and int(marker.name) > 0
        and not marker.is_symlink() and marker.is_file()
    ]


def local_workflow_running(pid_dir: Path) -> bool:
    """Inspect controlled numeric PID markers and prune only dead-process files.

    Other files, directories and symlinks are left alone. A process whose state
    cannot be inspected is kept conservatively; an unknown state is not proof
    that a second workflow can safely start.
    """
    try:
        markers = pid_markers(pid_dir)
    except OSError:
        return True
    running = False
    for marker in markers:
        # Inspect every marker, so live workers do not leave other dead
        # markers behind through short-circuit evaluation.
        running = prune_exited_marker(marker) or running
    return running


def prune_exited_marker(marker: Path) -> bool:
    """Prune a dead owned PID marker; return True for a live/unknown process.

    Signal delivery is asynchronous. Cancellation must retain tracking until
    the process actually exits, including processes that ignore SIGTERM.
    Callers select positive regular markers with ``pid_markers`` first.
    """
    try:
        process = psutil.Process(int(marker.name))
        alive = process.is_running() and process.status() not in {
            psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD
        }
    except (psutil.NoSuchProcess, ValueError, OverflowError):
        alive = False
    except psutil.AccessDenied:
        alive = True
    if not alive:
        try:
            marker.unlink(missing_ok=True)
        except OSError:
            # An unremovable marker does not make a proven-dead process live.
            pass
    return alive
