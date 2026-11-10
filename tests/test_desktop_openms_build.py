"""Exercise desktop build routing against modern and legacy OpenMS checkouts."""

import json
import os
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = yaml.safe_load((ROOT / ".github/workflows/openms.yml").read_text())
ROUTE = next(
    step["run"]
    for step in WORKFLOW["jobs"]["build"]["steps"]
    if step.get("id") == "route"
)


@pytest.mark.parametrize(
    "platform,preset,thirdparty",
    [
        ("Linux", "linux-x64-release", "Linux/x86_64"),
        ("macOS", "macos-arm64-release", "MacOS/arm64"),
        ("Windows", "windows-x64-release", "Windows/x86_64"),
    ],
)
def test_current_openms_uses_its_platform_preset_without_removed_ctest_scripts(
    tmp_path, platform, preset, thirdparty
):
    (tmp_path / "CMakePresets.json").write_text(
        json.dumps({"configurePresets": [{"name": preset}]})
    )
    (tmp_path / "vcpkg.json").write_text("{}")
    result, output = run_route(tmp_path, platform, thirdparty)
    assert result.returncode == 0, result.stdout + result.stderr
    assert output == {"route": "vcpkg", "preset": preset}


@pytest.mark.parametrize("platform", ["Linux", "macOS", "Windows"])
def test_legacy_openms_retains_its_supported_build_route(tmp_path, platform):
    if platform != "Windows":
        scripts = tmp_path / "tools/ci"
        scripts.mkdir(parents=True)
        dependency_script = "deps-macos.sh" if platform == "macOS" else "deps-ubuntu.sh"
        for name in [dependency_script, "cibuild.cmake", "cipackage.cmake"]:
            (scripts / name).touch()
    result, output = run_route(tmp_path, platform, "fixture")
    assert result.returncode == 0, result.stdout + result.stderr
    assert output["route"] == ("contrib" if platform == "Windows" else "cibuild")


def test_missing_supported_build_route_fails_before_installing_dependencies(tmp_path):
    result, output = run_route(tmp_path, "Linux", "fixture")
    assert result.returncode != 0
    assert output == {}
    assert "neither a linux-x64-release preset nor" in result.stdout


def run_route(directory, platform, thirdparty):
    # Route discovery should inspect checkout files, without contacting GitHub
    # or actually fetching the third-party tools after it selects the route.
    binaries = directory / "bin"
    binaries.mkdir()
    git = binaries / "git"
    git.write_text("#!/bin/sh\nexit 0\n")
    git.chmod(0o755)
    for name in [thirdparty, "All"]:
        tools = directory / "THIRDPARTY" / name
        tools.mkdir(parents=True, exist_ok=True)
        (tools / "tool").touch()
    output = directory / "output"
    result = subprocess.run(
        ["bash", "-e", "-o", "pipefail"],
        input=ROUTE.replace("${{ matrix.thirdparty }}", thirdparty),
        cwd=directory,
        env=os.environ
        | {
            "PATH": f"{binaries}{os.pathsep}{os.environ['PATH']}",
            "RUNNER_OS": platform,
            "REPO": "OpenMS/OpenMS",
            "SHA": "test-source",
            "GITHUB_OUTPUT": str(output),
        },
        capture_output=True,
        text=True,
        check=False,
    )
    values = (
        dict(line.split("=", 1) for line in output.read_text().splitlines())
        if output.exists()
        else {}
    )
    return result, values
