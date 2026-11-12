"""
Keeps every run of a workflow as its own named result set.

A workflow that sets ``multi_run = True`` writes each run into
``<workflow_dir>/runs/<run_id>/`` instead of overwriting
``<workflow_dir>/results/``:

    runs/<run_id>/run.json      name, status and timestamps
    runs/<run_id>/params.json   the parameters the run was started with
    runs/<run_id>/results/      everything execution() wrote
    runs/<run_id>/logs/         the run's logs, copied when it ends

This module must stay importable without Streamlit: the RQ worker
(``tasks.py``) uses it to mark a queued run finished or failed.
"""

import json
import re
import shutil
import stat
import time
import uuid
import zipfile
import zlib
from datetime import datetime
from pathlib import Path
from typing import Optional

from .file_lock import file_lock

RUN_FILE = "run.json"
MAX_NAME_LENGTH = 100
# Largest unpacked size import_zip() accepts by default
MAX_IMPORT_BYTES = 20 * 1024**3
# Most files a run ZIP may add, and the largest run.json / params.json read
MAX_IMPORT_MEMBERS = 100_000
MAX_RUN_FILE_BYTES = 1024**2
MAX_PARAMS_FILE_BYTES = 64 * 1024**2


class RunImportError(ValueError):
    """A ZIP that can not be imported as a run; the message says why."""

# Run states written to run.json. A run whose process died without reaching
# finish() stays "running"; the UI shows it as interrupted once the workflow
# is no longer running.
RUNNING = "running"
FINISHED = "finished"
FAILED = "failed"
CANCELLED = "cancelled"


class RunManager:
    """Creates, lists, renames, deletes and archives the runs of one workflow."""

    def __init__(self, workflow_dir: Path) -> None:
        self.workflow_dir = Path(workflow_dir)
        self.runs_dir = Path(self.workflow_dir, "runs")

    # Paths

    def run_dir(self, run_id: str) -> Path:
        return Path(self.runs_dir, run_id)

    def results_dir(self, run_id: str) -> Path:
        return Path(self.runs_dir, run_id, "results")

    def params_file(self, run_id: str) -> Path:
        return Path(self.runs_dir, run_id, "params.json")

    # Reading

    def get(self, run_id: Optional[str]) -> Optional[dict]:
        """The run's metadata, or None if there is no such run."""
        if not run_id or not _is_safe_id(run_id):
            return None
        try:
            with open(Path(self.run_dir(run_id), RUN_FILE), "r", encoding="utf-8") as f:
                run = json.load(f)
        except (OSError, ValueError):
            return None
        if not isinstance(run, dict):
            return None
        run["id"] = run_id
        return run

    def list_runs(self, status: Optional[str] = None) -> list[dict]:
        """All runs, newest first. Pass ``status`` to keep only runs in that state."""
        try:
            directories = list(self.runs_dir.iterdir())
        except (FileNotFoundError, NotADirectoryError):
            return []
        runs = [self.get(d.name) for d in directories if not d.name.startswith(".") and d.is_dir()]
        runs = [r for r in runs if r is not None]
        if status is not None:
            runs = [r for r in runs if r.get("status") == status]
        return sorted(runs, key=lambda r: (r.get("created", ""), r["id"]), reverse=True)

    def latest(self, status: Optional[str] = FINISHED) -> Optional[str]:
        """Id of the newest run (in ``status``, if given), or None."""
        runs = self.list_runs(status)
        return runs[0]["id"] if runs else None

    def params(self, run_id: str) -> dict:
        """The parameters a run was started with (empty if not recorded)."""
        try:
            with open(self.params_file(run_id), "r", encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    # Writing

    def new_run(self, name: str = "", params_file: Optional[Path] = None) -> str:
        """
        Create the directory for a new run and return its id.

        Args:
            name: Display name. Defaults to "Run <n>".
            params_file: The workflow's params.json, copied into the run so the
                results can be read with the parameters that produced them.
        """
        with self._lock():
            run_id = self._reserve_run()
            try:
                self.results_dir(run_id).mkdir()
                name = _clean_name(name) or self.default_name()
                if params_file is not None:
                    shutil.copy2(params_file, self.params_file(run_id))
                self._write(run_id, {"name": name, "created": _now(), "status": RUNNING})
            except BaseException:
                shutil.rmtree(self.run_dir(run_id), ignore_errors=True)
                raise
        return run_id

    def _reserve_run(self) -> str:
        """Claim a new ID with mkdir, including against older uncoordinated writers."""
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        base = datetime.now().strftime("%Y%m%d-%H%M%S")
        run_id, n = base, 1
        while True:
            try:
                self.run_dir(run_id).mkdir()
                return run_id
            except FileExistsError:
                n += 1
                run_id = f"{base}-{n}"

    def finish(self, run_id: Optional[str], status: str) -> None:
        """Record how a run ended and keep a copy of its logs."""
        if not run_id or not _is_safe_id(run_id):
            return
        with self._lock():
            run = self.get(run_id)
            # Check terminal state under the same lock as the update. A stale
            # cancellation request cannot cancel a run that already finished.
            if run is None or run.get("status") == CANCELLED or (
                status == CANCELLED and run.get("status") != RUNNING
            ):
                return
            run.update(status=status, finished=_now())
            self._write(run_id, run)
        self._copy_logs(run_id, status)

    def rename(self, run_id: str, name: str) -> bool:
        """Set a run's display name. False if the name is empty or too long."""
        name = _clean_name(name)
        if not _is_safe_id(run_id) or not name or len(name) > MAX_NAME_LENGTH:
            return False
        with self._lock():
            run = self.get(run_id)
            if run is None:
                return False
            run["name"] = name
            self._write(run_id, run)
        return True

    def delete(self, run_id: str) -> bool:
        """Delete a run and all its files."""
        if not _is_safe_id(run_id):
            return False
        removed = self.workflow_dir / f".run-delete-{uuid.uuid4().hex}"
        with self._lock():
            if self.get(run_id) is None:
                return False
            self.run_dir(run_id).rename(removed)
            Path(self.workflow_dir, "downloads", f"{run_id}.zip").unlink(missing_ok=True)
        shutil.rmtree(removed, ignore_errors=True)
        return True

    def zip(self, run_id: str, exclude: Optional[list[str]] = None) -> Path:
        """
        Pack a run's results, parameters and logs into
        ``<workflow_dir>/downloads/<run_id>.zip`` and return the path.

        Args:
            exclude: File or folder names to leave out (glob patterns allowed),
                e.g. caches such as ``["insight_cache"]``.
        """
        from fnmatch import fnmatch

        exclude = list(exclude or [])
        run_dir = self.run_dir(run_id)
        zip_path = Path(self.workflow_dir, "downloads", f"{run_id}.zip")
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = zip_path.with_name(f".{run_id}-{uuid.uuid4().hex}.zip")
        try:
            if self.get(run_id) is None:
                raise FileNotFoundError("The selected run is no longer available.")
            # A second tab can prepare a ZIP too. Readers retain the previous
            # complete archive until this private replacement is ready.
            with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED) as zf:
                for path in sorted(run_dir.rglob("*")):
                    rel = path.relative_to(run_dir)
                    if path.is_file() and not any(
                        fnmatch(part, p) for part in rel.parts for p in exclude
                    ):
                        zf.write(path, rel)
            with self._lock():
                if self.get(run_id) is None:
                    raise FileNotFoundError("The run was deleted while preparing its archive.")
                temporary.replace(zip_path)
        finally:
            temporary.unlink(missing_ok=True)
        return zip_path

    def import_zip(
        self, source, name: str = "", max_bytes: int = MAX_IMPORT_BYTES
    ) -> str:
        """
        Add a run from a ZIP made by zip() and return its new id.

        The ZIP holds ``results/`` and, as zip() writes them, ``run.json``,
        ``params.json`` and ``logs/``; a ZIP of the run folder itself (all of
        that under one top folder) works too. The run keeps its name, status,
        times and parameters, and gets a new id, so importing the same ZIP
        twice gives two runs.

        Every member is checked before anything is written: absolute paths,
        ``..``, drive letters, backslashes and symbolic links are refused, as
        as are names Windows reserves (``NUL``, ``CON``, ...), a ZIP that
        unpacks to more than ``max_bytes`` or holds more than
        MAX_IMPORT_MEMBERS files, and an oversized run.json or params.json.
        A run that was still running when it was zipped is imported as failed. Unpacking goes to a
        temporary folder that becomes the run only when it is complete.

        Args:
            source: Path or binary file object (e.g. a Streamlit upload).
            name: Display name; defaults to the name in run.json, then
                "Imported run".
            max_bytes: Largest unpacked size accepted.

        Raises:
            RunImportError: The file is not such a ZIP, or is unsafe.
        """
        try:
            archive = zipfile.ZipFile(source)
        except (zipfile.BadZipFile, OSError) as e:
            raise RunImportError("This is not a ZIP file.") from e
        try:
            with archive:
                return self._import_archive(archive, name, max_bytes)
        except (zipfile.BadZipFile, zlib.error, EOFError) as e:
            raise RunImportError(f"The ZIP is damaged ({e}).") from e
        except RuntimeError as e:
            # Encrypted members, an unsupported compression method, or JSON
            # nested too deeply to parse
            raise RunImportError(f"The ZIP can not be read ({e}).") from e
        except OSError as e:
            raise RunImportError(f"The ZIP could not be unpacked ({e}).") from e

    def _import_archive(self, archive: zipfile.ZipFile, name: str, max_bytes: int) -> str:
        members = [m for m in archive.infolist() if not m.is_dir()]
        for m in archive.infolist():
            _check_member(m)
        if len(members) > MAX_IMPORT_MEMBERS:
            raise RunImportError(
                f"The ZIP holds {len(members)} files, more than the "
                f"{MAX_IMPORT_MEMBERS} allowed."
            )
        total = sum(m.file_size for m in members)
        if total > max_bytes:
            raise RunImportError(
                f"The ZIP unpacks to {total / 1024**3:.1f} GB, more than "
                f"the {max_bytes / 1024**3:.1f} GB allowed."
            )
        prefix = _run_prefix([m.filename for m in members])
        if prefix is None:
            raise RunImportError(
                "This ZIP holds no run: a run's ZIP has a results/ "
                "folder or a run.json."
            )

        for m in members:
            limit = {
                f"{prefix}{RUN_FILE}": MAX_RUN_FILE_BYTES,
                f"{prefix}params.json": MAX_PARAMS_FILE_BYTES,
            }.get(m.filename)
            if limit is not None and m.file_size > limit:
                raise RunImportError(
                    f"Refused: {m.filename} is larger than "
                    f"{limit / 1024**2:.0f} MB."
                )

        meta = {}
        if f"{prefix}{RUN_FILE}" in archive.namelist():
            try:
                meta = json.loads(archive.read(f"{prefix}{RUN_FILE}"))
            except ValueError:
                meta = {}
        if not isinstance(meta, dict):
            meta = {}
        status = meta.get("status")
        if status == RUNNING:
            # Zipped while it ran, or after its process died: its results
            # may be incomplete.
            status = FAILED
        elif status not in (FINISHED, FAILED, CANCELLED):
            status = FINISHED

        tmp = Path(self.workflow_dir, f".run-import-{uuid.uuid4().hex[:8]}")
        tmp.mkdir(parents=True)
        try:
            for m in members:
                rel = m.filename[len(prefix):]
                if not rel or rel == RUN_FILE:
                    continue
                target = Path(tmp, *rel.split("/"))
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(m) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            params = Path(tmp, "params.json")
            if params.exists():
                try:
                    if not isinstance(json.loads(params.read_text(encoding="utf-8")), dict):
                        raise ValueError
                except ValueError:
                    params.unlink()
            Path(tmp, "results").mkdir(exist_ok=True)

            run = {
                "name": _clean_name(name) or _clean_name(meta.get("name")) or "Imported run",
                "created": _iso_time(meta.get("created")) or _now(),
                "status": status,
                "imported": _now(),
            }
            if _iso_time(meta.get("finished")):
                run["finished"] = _iso_time(meta["finished"])
            with self._lock():
                run["name"] = self._unique_name(run["name"])
                self.runs_dir.mkdir(parents=True, exist_ok=True)
                base = datetime.now().strftime("%Y%m%d-%H%M%S")
                run_id, n = base, 1
                while self.run_dir(run_id).exists():
                    n += 1
                    run_id = f"{base}-{n}"
                with open(Path(tmp, RUN_FILE), "w", encoding="utf-8") as f:
                    json.dump(run, f, indent=2)
                tmp.rename(self.run_dir(run_id))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        return run_id

    def _unique_name(self, name: str) -> str:
        """``name``, or ``name (n)`` if a run has it, within MAX_NAME_LENGTH."""
        taken = {r.get("name") for r in self.list_runs()}
        name = name[:MAX_NAME_LENGTH]
        if name not in taken:
            return name
        n = 2
        while True:
            suffix = f" ({n})"
            candidate = name[: MAX_NAME_LENGTH - len(suffix)].rstrip() + suffix
            if candidate not in taken:
                return candidate
            n += 1

    def adopt_legacy_results(self) -> Optional[str]:
        """
        Turn a ``<workflow_dir>/results`` left by a single-run version of the
        workflow into a run named "Earlier run", so switching a workflow to
        multi-run keeps what users already computed. Returns the new run id,
        or None if there was nothing to adopt.
        """
        legacy = Path(self.workflow_dir, "results")
        # Empty workspaces are checked before opening a lock so merely viewing
        # an unused workflow does not create run state.
        try:
            if not any(legacy.iterdir()):
                return None
        except (FileNotFoundError, NotADirectoryError):
            return None
        with self._lock():
            try:
                if not any(legacy.iterdir()) or self.list_runs():
                    return None
                created = datetime.fromtimestamp(legacy.stat().st_mtime)
            except (FileNotFoundError, NotADirectoryError):
                return None
            run_id = self._reserve_run()
            try:
                shutil.move(str(legacy), str(self.results_dir(run_id)))
            except OSError:
                shutil.rmtree(self.run_dir(run_id), ignore_errors=True)
                return None
            self._write(
                run_id,
                {
                    "name": "Earlier run",
                    "created": created.isoformat(timespec="seconds"),
                    "finished": created.isoformat(timespec="seconds"),
                    "status": FINISHED,
                },
            )
            params = Path(self.workflow_dir, "params.json")
            try:
                shutil.copy2(params, self.params_file(run_id))
            except FileNotFoundError:
                pass
            return run_id

    # Internals

    def default_name(self) -> str:
        return f"Run {len(self.list_runs()) + 1}"

    def _lock(self):
        # Stable across RunManager instances and worker processes; deletion of
        # a run must not remove the lock another process is still holding.
        return file_lock(self.runs_dir / ".locks" / "metadata.lock")

    def _copy_logs(self, run_id: str, status: str) -> None:
        """Copy logs off-lock, then publish only if this run still exists."""
        staged = self.workflow_dir / f".run-logs-{uuid.uuid4().hex}"
        retired = self.workflow_dir / f".run-logs-old-{uuid.uuid4().hex}"
        try:
            try:
                shutil.copytree(self.workflow_dir / "logs", staged)
            except FileNotFoundError:
                # A new run may have cleared the shared log folder. Its
                # disappearance cannot turn the completed tool into a failure.
                return
            except shutil.Error as error:
                # copytree aggregates per-file failures and drops errno. A
                # missing source after the copy attempt identifies the same
                # log-clearing race; existing sources retain genuine errors.
                failures = error.args[0]
                if not isinstance(failures, list) or not failures or any(
                    not isinstance(item, (tuple, list)) or len(item) != 3
                    or Path(item[0]).exists()
                    for item in failures
                ):
                    raise
                return
            with self._lock():
                current = self.get(run_id)
                if current is None or current.get("status") != status:
                    return
                target = self.run_dir(run_id) / "logs"
                if target.exists():
                    target.rename(retired)
                try:
                    staged.rename(target)
                except OSError:
                    if retired.exists() and not target.exists():
                        retired.rename(target)
                    raise
        finally:
            shutil.rmtree(staged, ignore_errors=True)
            shutil.rmtree(retired, ignore_errors=True)

    def _write(self, run_id: str, run: dict) -> None:
        run = {k: v for k, v in run.items() if k != "id"}
        path = Path(self.run_dir(run_id), RUN_FILE)
        tmp = path.with_name(f".{RUN_FILE}-{uuid.uuid4().hex}.tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(run, f, indent=2)
            # Windows can briefly deny replacement while a status reader has
            # run.json open. Keep the previous complete JSON during retries.
            for attempt in range(20):
                try:
                    tmp.replace(path)
                    break
                except PermissionError:
                    if attempt == 19:
                        raise
                    time.sleep(0.05)
        finally:
            tmp.unlink(missing_ok=True)


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _clean_name(name: Optional[str]) -> str:
    return " ".join(str(name or "").split())


def _is_safe_id(run_id: str) -> bool:
    # Run ids come from session state and query parameters; never let one
    # leave the runs directory.
    return bool(run_id) and "/" not in run_id and "\\" not in run_id and run_id not in (".", "..")


def _iso_time(value) -> Optional[str]:
    """``value`` if it is an ISO timestamp string, else None."""
    if not isinstance(value, str):
        return None
    try:
        datetime.fromisoformat(value)
    except ValueError:
        return None
    return value


# Names Windows maps to devices, with or without an extension (NUL, nul.txt)
_WINDOWS_RESERVED = re.compile(r"^(con|prn|aux|nul|com[1-9]|lpt[1-9])(\..*)?$", re.I)


def _check_member(member: zipfile.ZipInfo) -> None:
    """Refuse a ZIP member that could land outside the folder it is unpacked to."""
    name = member.filename
    parts = name.replace("\\", "/").split("/")
    if (
        "\\" in name
        or name.startswith("/")
        or ".." in parts
        or (len(name) > 1 and name[1] == ":")
        or "\x00" in name
        or any(_WINDOWS_RESERVED.match(p.rstrip(" .")) for p in parts if p)
    ):
        raise RunImportError(f"Refused: the ZIP has an unsafe path ({name!r}).")
    mode = member.external_attr >> 16
    if stat.S_ISLNK(mode):
        raise RunImportError(f"Refused: the ZIP holds a symbolic link ({name!r}).")


def _run_prefix(names: list[str]) -> Optional[str]:
    """
    Where the run sits in the ZIP: "" when results/ is at the top, "<dir>/"
    when everything is in one folder, None when it holds no run. A run is
    a results/ folder or a run.json: zip() leaves out an empty results/.
    """
    def is_run(prefix: str) -> bool:
        return any(
            n.startswith(f"{prefix}results/") or n == f"{prefix}{RUN_FILE}"
            for n in names
        )

    if is_run(""):
        return ""
    tops = {n.split("/", 1)[0] for n in names}
    if len(tops) == 1:
        top = tops.pop() + "/"
        if is_run(top):
            return top
    return None


def format_time(iso: Optional[str]) -> str:
    """'2026-10-30T09:12:05' -> '2026-10-30 09:12' (empty for None)."""
    if not iso:
        return ""
    try:
        return datetime.fromisoformat(iso).strftime("%Y-%m-%d %H:%M")
    except ValueError:
        return iso
