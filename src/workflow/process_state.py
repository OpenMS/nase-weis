"""Check live local workflows without mistaking orphaned PID files for workers."""

from pathlib import Path

import psutil


def local_workflow_running(pid_dir: Path) -> bool:
    """Inspect controlled numeric PID markers and prune only dead-process files.

    Other files, directories and symlinks are left alone. A process whose state
    cannot be inspected is kept conservatively; an unknown state is not proof
    that a second workflow can safely start.
    """
    try:
        markers = list(Path(pid_dir).iterdir())
    except FileNotFoundError:
        return False
    except OSError:
        return True
    running = False
    for marker in markers:
        if not marker.name.isdecimal() or marker.is_symlink() or not marker.is_file():
            continue
        pid = int(marker.name)
        try:
            process = psutil.Process(pid)
            alive = process.is_running() and process.status() not in {
                psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD
            }
        except (psutil.NoSuchProcess, ValueError, OverflowError):
            alive = False
        except psutil.AccessDenied:
            alive = True
        if alive:
            running = True
        else:
            try:
                marker.unlink(missing_ok=True)
            except OSError:
                # The UI can still report that no process is alive if this
                # stale marker cannot be removed from the shared filesystem.
                pass
    return running
