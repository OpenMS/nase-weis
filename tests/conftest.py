"""Opt-in requirements for the complete NASEWEIS regression gate."""

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def pytest_addoption(parser):
    parser.addoption(
        "--strict-regressions",
        action="store_true",
        help="Require TOPP tools and fail if any regression check skips or xfails.",
    )
    parser.addoption(
        "--naseweis-baseline",
        action="store_true",
        help="Compare every NME1 match/mapping with the pinned nightly reference.",
    )


def pytest_sessionstart(session):
    config = session.config
    if not (
        config.getoption("--strict-regressions")
        or config.getoption("--naseweis-baseline")
    ):
        return
    missing = [
        name
        for name in ["DecoyDatabase", "NucleicAcidSearchEngine"]
        if not shutil.which(name)
    ]
    if missing:
        raise pytest.UsageError(
            "Missing required TOPP tools: "
            + ", ".join(missing)
            + ". Run: python scripts/check_regressions.py --container"
        )
    if not config.getoption("--naseweis-baseline"):
        return
    manifest = json.loads(
        (ROOT / "tests/data/naseweis/nightly/manifest.json").read_text()
    )
    help_text = subprocess.run(
        ["NucleicAcidSearchEngine", "--help"],
        capture_output=True,
        text=True,
        check=True,
        timeout=30,
    )
    version = re.search(r"Version: (\S+)", help_text.stdout + help_text.stderr)
    if not version or version.group(1) != manifest["openms_version"]:
        raise pytest.UsageError(
            "Reference results require the pinned OpenMS nightly runtime. "
            "Run: python scripts/check_regressions.py --container"
        )
    if (
        manifest.get("revision")
        and f"Revision: {manifest['revision']}"
        not in help_text.stdout + help_text.stderr
    ):
        raise pytest.UsageError(
            "OpenMS source revision does not match the nightly reference."
        )
    import pyopenms

    if pyopenms.__version__ != manifest["pyopenms_version"]:
        raise pytest.UsageError(
            "pyOpenMS version does not match the reference manifest."
        )
    for relative, expected in manifest["inputs"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise pytest.UsageError(
                f"Reference input changed: {relative}. Review and explicitly regenerate its baseline."
            )


def pytest_sessionfinish(session, exitstatus):
    if session.config.getoption("--strict-regressions"):
        reporter = session.config.pluginmanager.get_plugin("terminalreporter")
        incomplete = sum(
            len(reporter.stats.get(status, []))
            for status in ["skipped", "xfailed", "xpassed"]
        )
        if incomplete:
            reporter.write_sep(
                "!", f"Regression gate failed: {incomplete} skipped/xfail checks"
            )
            session.exitstatus = pytest.ExitCode.TESTS_FAILED
