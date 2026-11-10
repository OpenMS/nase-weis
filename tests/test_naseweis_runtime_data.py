"""Keep generated native runtime data when source/build trees are removed."""

import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("dockerfile", ["Dockerfile", "Dockerfile.arm"])
@pytest.mark.parametrize("with_managed_reader", [True, False])
def test_installed_runtime_assets_and_source_examples_survive_staging(
    tmp_path, dockerfile, with_managed_reader
):
    source = tmp_path / "source share"
    source.mkdir()
    (source / "examples").mkdir()
    (source / "examples/rna.fasta").write_text(">rna\nACGU\n")
    (source / "schema.xsd").write_text("source schema")
    build = tmp_path / "build"
    installed = build / "installed-share"
    (installed / "TOOLS").mkdir(parents=True)
    (installed / "TOOLS/OpenMS.tsv").write_text("NucleicAcidSearchEngine\tRNA\n")
    expected = {"TOOLS/OpenMS.tsv": "NucleicAcidSearchEngine\tRNA\n"}
    # Older releases also expose the generated registry directly in the build
    # tree; preserving it alone must not stand in for installing all runtime data.
    registry = build / "share/OpenMS/TOOLS"
    registry.mkdir(parents=True)
    (registry / "OpenMS.tsv").write_text(expected["TOOLS/OpenMS.tsv"])
    if with_managed_reader:
        managed = installed / "openms_thermo_bridge/managed"
        managed.mkdir(parents=True)
        for filename in [
            "ThermoWrapperManaged.dll",
            "ThermoWrapperManaged.runtimeconfig.json",
            "ThermoFisher.CommonCore.RawFileReader.dll",
            "THIRD-PARTY-NOTICES.txt",
        ]:
            (managed / filename).write_text(filename)
            expected[f"openms_thermo_bridge/managed/{filename}"] = filename
        (installed / "LICENSES").mkdir()
        (installed / "LICENSES/ThermoRawFileReader-License.doc").write_text("license")
        expected["LICENSES/ThermoRawFileReader-License.doc"] = "license"

    # Model CMake's share-component contract; the Docker step itself performs
    # the source/install merge. Native validation separately exercises CMake
    # against the pinned OpenMS build and opens a real RAW file with its bridge.
    binaries = tmp_path / "bin"
    binaries.mkdir()
    cmake = binaries / "cmake"
    cmake.write_text(
        f"#!{sys.executable}\n"
        "import pathlib, shutil, sys\n"
        "assert sys.argv[1:2] == ['--install']\n"
        "assert sys.argv[3:6] == ['--component', 'share', '--prefix']\n"
        "source = pathlib.Path(sys.argv[2]) / 'installed-share'\n"
        "destination = pathlib.Path(sys.argv[6]) / 'share/OpenMS'\n"
        "shutil.copytree(source, destination)\n"
    )
    cmake.chmod(0o755)
    destination = tmp_path / "runtime share"
    command = staging_command(ROOT / dockerfile)
    command = command.replace("OpenMS/share/OpenMS", shlex.quote(str(source)))
    command = command.replace(
        "openms-build/share/OpenMS/TOOLS", shlex.quote(str(registry))
    )
    command = command.replace("/openms-build", shlex.quote(str(build)))
    command = command.replace("/openms/share", shlex.quote(str(destination)))
    subprocess.run(
        ["bash", "-e", "-o", "pipefail", "-c", command],
        env=os.environ | {"PATH": f"{binaries}{os.pathsep}{os.environ['PATH']}"},
        check=True,
        capture_output=True,
        text=True,
    )
    # Data must be self-contained after the original source/build are gone.
    shutil.rmtree(source)
    shutil.rmtree(build)
    assert (destination / "examples/rna.fasta").read_text() == ">rna\nACGU\n"
    assert (destination / "schema.xsd").read_text() == "source schema"
    for relative, content in expected.items():
        assert (destination / relative).read_text() == content


def staging_command(dockerfile):
    lines = dockerfile.read_text().splitlines()
    start = next(
        index
        for index, line in enumerate(lines)
        if line.startswith("RUN cp -r OpenMS/share/OpenMS /openms/share")
    )
    end = next(
        index
        for index in range(start, len(lines))
        if lines[index] == "RUN rm -rf OpenMS"
    )
    return "\n".join(
        line.removeprefix("RUN ").removesuffix("\\")
        for line in lines[start:end]
        if not line.startswith("#")
    )
