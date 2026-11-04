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
