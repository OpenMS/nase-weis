"""Keep future app pages/modules inside the application coverage gate."""

import configparser
import fnmatch
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_all_app_pages_and_owned_analysis_modules_have_coverage():
    spec = yaml.safe_load((ROOT / "app.yaml").read_text())
    config = configparser.ConfigParser()
    config.read(ROOT / ".coveragerc")
    includes = config["report"]["include"].split()
    required = {"app.py", *spec["pages"]}
    required.update(
        path
        for path in spec["sync"]["app_owned"]
        if path.startswith("src/") and path.endswith(".py") and "*" not in path
    )
    missing = sorted(
        path
        for path in required
        if not any(fnmatch.fnmatch(path, pattern) for pattern in includes)
    )
    assert not missing, (
        f"New application code is outside the regression coverage gate: {missing}"
    )


def test_main_ci_calls_the_strict_gate_without_duplicate_event_triggers():
    # BaseLoader preserves GitHub's "on" key instead of YAML 1.1 treating it
    # as a boolean. Check the event graph as well as the script invocation.
    ci = yaml.load(
        (ROOT / ".github/workflows/ci.yml").read_text(), Loader=yaml.BaseLoader
    )
    job = ci["jobs"]["regressions"]
    assert "pull_request" in ci["on"]
    assert "main" in ci["on"]["push"]["branches"]
    assert job["uses"] == "./.github/workflows/naseweis-regressions.yml"
    reusable = yaml.load((ROOT / job["uses"]).read_text(), Loader=yaml.BaseLoader)
    assert "workflow_call" in reusable["on"]
    assert "workflow_dispatch" in reusable["on"]
    assert not {"push", "pull_request"}.intersection(reusable["on"])
    steps = reusable["jobs"]["pinned-openms-regressions"]["steps"]
    assert any(
        step.get("run") == "python3 scripts/check_regressions.py --container"
        for step in steps
    )
    reports = next(
        step
        for step in steps
        if step.get("uses", "").startswith("actions/upload-artifact@")
    )
    assert reports["if"] == "always()"
    assert reports["with"]["path"] == ".local/regressions/*.xml"
