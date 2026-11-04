"""Check the registry repair used by the container's real-tool build smoke test."""

import pytest

from docker.restore_tool_registry import restore_registry


def test_registry_matches_the_installed_tools(tmp_path):
    source = tmp_path / "executables.cmake"
    source.write_text("""# tool declarations
openms_topp_tool(DecoyDatabase "Database processing")
openms_topp_tool(NucleicAcidSearchEngine "RNA search")
openms_topp_tool(NotInstalled "Other")
""")
    binaries = tmp_path / "bin"
    binaries.mkdir()
    for name in ["DecoyDatabase", "NucleicAcidSearchEngine"]:
        (binaries / name).touch()
    output = tmp_path / "share/TOOLS/OpenMS.tsv"
    assert restore_registry(source, binaries, output) == 2
    assert "NucleicAcidSearchEngine\tRNA search" in output.read_text()
    assert "NotInstalled" not in output.read_text()


def test_registry_repair_refuses_a_mismatching_release(tmp_path):
    source = tmp_path / "executables.cmake"
    source.write_text("set(TOPP_TOOLS OldTool)")
    with pytest.raises(ValueError, match="do not contain both"):
        restore_registry(source, tmp_path, tmp_path / "registry.tsv")
