"""The MS data panels say why they draw nothing for a missing or damaged mzML."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from src.viewer import panels


class _Streamlit:
    def __init__(self):
        self.shown = []

    def markdown(self, text, **_):
        self.shown.append(("markdown", text))

    def info(self, text):
        self.shown.append(("info", text))

    def error(self, text):
        self.shown.append(("error", text))


@pytest.fixture
def st(monkeypatch):
    fake = _Streamlit()
    monkeypatch.setattr(panels, "st", fake)
    return fake


@pytest.mark.parametrize("panel", ["file_scans", "file_spectrum", "file_chromatogram", "file_peak_map"])
def test_damaged_mzml_shows_an_error(st, tmp_path, panel):
    damaged = tmp_path / "broken.mzML"
    damaged.write_text("<mzML>not really")
    src = SimpleNamespace(mzml=lambda: damaged, notes=[])
    getattr(panels, panel)(SimpleNamespace(), src)
    kinds = [kind for kind, _ in st.shown]
    assert kinds == ["markdown", "error"]
    assert "broken.mzML" in st.shown[1][1] and "damaged or not mzML" in st.shown[1][1]


def test_missing_mzml_says_why(st):
    src = SimpleNamespace(mzml=lambda: None, notes=["the file was removed", "no copy kept"])
    panels.file_scans(SimpleNamespace(), src)
    assert st.shown == [
        ("markdown", "**Scans**"),
        ("info", "The spectra of this acquisition are not available: the file was removed; no copy kept"),
    ]
