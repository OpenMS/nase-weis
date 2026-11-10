"""Viewer provenance, file boundaries and immutable revision-aware caches."""

import json

import pandas as pd
import pytest

from src import naseweis_viewers as viewers


def test_snapshot_uses_its_exact_saved_values(tmp_path):
    path = tmp_path / "params.json"
    params = {
        "mzML-files": "/saved/original.mzML",
        "FeatureFinderMetaboIdent": {"extract:mz_window": 12.5},
    }
    path.write_text(json.dumps(params))
    assert viewers.read_run_snapshot(path) == params


@pytest.mark.parametrize("content", [None, "invalid JSON", "[]", "null"])
def test_missing_or_invalid_snapshot_has_an_actionable_error(tmp_path, content):
    path = tmp_path / "params.json"
    if content is not None:
        path.write_text(content)
    with pytest.raises(ValueError, match="snapshot"):
        viewers.read_run_snapshot(path)


def test_acquisition_boundaries_and_file_revisions(tmp_path):
    workflow = tmp_path / "workspace" / "naseweis"
    workflow.mkdir(parents=True)
    inside = workflow / "spectra.mzML"
    inside.touch()
    first = viewers.acquisition_revision(str(inside), workflow, True)
    stat = inside.stat()
    assert first == (
        str(inside.resolve()),
        stat.st_mtime_ns,
        0,
        stat.st_ctime_ns,
        stat.st_ino,
        stat.st_dev,
    )
    inside.write_bytes(b"changed")
    assert viewers.acquisition_revision(str(inside), workflow, True) != first
    outside = tmp_path / "outside.mzML"
    outside.touch()
    assert viewers.acquisition_revision(str(outside), workflow)[0] == str(outside)
    with pytest.raises(ValueError, match="outside this online workspace"):
        viewers.acquisition_revision(str(outside), workflow, True)
    symlink = workflow / "linked.mzML"
    symlink.symlink_to(outside)
    with pytest.raises(ValueError, match="outside this online workspace"):
        viewers.acquisition_revision(str(symlink), workflow, True)


@pytest.mark.parametrize("selected", [None, "", [], ["/path.mzML"]])
def test_missing_or_ambiguous_saved_acquisition_is_rejected(tmp_path, selected):
    with pytest.raises(ValueError, match="saved acquisition"):
        viewers.acquisition_revision(selected, tmp_path)


def test_discover_registered_and_legacy_inputs_without_changing_selection(tmp_path):
    workflow = tmp_path / "workspace" / "naseweis"
    uploaded = workflow / "input-files" / "mzML-files"
    legacy = workflow.parent / "mzML-files"
    uploaded.mkdir(parents=True)
    legacy.mkdir()
    files = [uploaded / "current.mzML", legacy / "old.MZML", uploaded / "thermo.raw"]
    for path in files:
        path.touch()
    vendor = uploaded / "bruker.d"
    vendor.mkdir()
    external = tmp_path / "external.mzML"
    external.touch()
    (uploaded / "external_files.txt").write_text(
        f"{external}\n\n{files[0]}\n{tmp_path / 'missing.mzML'}\n"
    )
    (uploaded / "unrelated.txt").touch()
    (uploaded / "external-link.mzML").symlink_to(external)
    params = {"mzML-files": str(files[0])}
    original = params.copy()
    assert set(viewers.input_acquisitions(workflow, params)) == {
        str(path.resolve()) for path in files + [external, vendor]
    }
    assert set(viewers.input_acquisitions(workflow, params, online=True)) == {
        str(path.resolve()) for path in files + [vendor]
    }
    assert params == original
    assert viewers.input_acquisitions(tmp_path / "empty" / "naseweis", {}) == []


def test_widget_namespaces_do_not_reuse_another_file_or_run(tmp_path):
    assert viewers.viewer_key(tmp_path / "run-a") == viewers.viewer_key(
        str(tmp_path / "run-a")
    )
    assert viewers.viewer_key(tmp_path / "run-a") != viewers.viewer_key(
        tmp_path / "run-b"
    )


def test_caches_follow_file_revisions_and_do_not_expose_mutable_tables(monkeypatch):
    calls = []

    def loader(path):
        calls.append(path)
        return pd.DataFrame({"index": [0], "native_id": [path]})

    viewers.cached_scan_index.clear()
    monkeypatch.setattr(viewers, "scan_index", loader)
    first = viewers.cached_scan_index(("file-a.mzML", 1, 10))
    first.loc[0, "native_id"] = "changed in one view"
    assert (
        viewers.cached_scan_index(("file-a.mzML", 1, 10)).loc[0, "native_id"]
        == "file-a.mzML"
    )
    viewers.cached_scan_index(("file-a.mzML", 2, 10))
    viewers.cached_scan_index(("file-b.mzML", 2, 10))
    assert calls == ["file-a.mzML", "file-a.mzML", "file-b.mzML"]
    viewers.cached_scan_index.clear()


@pytest.mark.parametrize(
    "loader,cache,arguments,result",
    [
        (
            "ms1_overview",
            "cached_overview",
            (),
            (pd.DataFrame({"rt": [1]}), pd.DataFrame({"trace": ["TIC"]})),
        ),
        ("spectrum_peaks", "cached_spectrum", (3,), pd.DataFrame({"mz": [100.0]})),
        (
            "xic",
            "cached_xic",
            (100.0, 10.0, "ppm", (0.0, 5.0)),
            pd.DataFrame({"intensity": [0.0]}),
        ),
        ("read_annotations", "cached_annotations", (), pd.DataFrame({"ion": ["w1"]})),
    ],
)
def test_each_cached_loader_passes_exact_arguments(
    monkeypatch, loader, cache, arguments, result
):
    from unittest.mock import Mock

    function = Mock(return_value=result)
    monkeypatch.setattr(viewers, loader, function)
    cached = getattr(viewers, cache)
    cached.clear()
    revision = ("selected-file.mzML", 10, 20)
    value = cached(revision, *arguments)
    function.assert_called_once_with("selected-file.mzML", *arguments)
    if isinstance(result, tuple):
        for actual, expected in zip(value, result):
            pd.testing.assert_frame_equal(actual, expected)
    else:
        pd.testing.assert_frame_equal(value, result)
    cached.clear()


@pytest.mark.parametrize(
    "feature,expected,boundaries",
    [
        ({"leftWidth": "10", "rightWidth": "20", "rt": 15}, (5, 25), (10, 20)),
        ({"rt_start": 10, "rt_end": 40, "rt": 20}, (-5, 55), (10, 40)),
        ({"leftWidth": "bad", "rightWidth": 20, "rt": 15}, (-105, 135), None),
        ({"leftWidth": 30, "rightWidth": 20, "rt": 15}, (-105, 135), None),
        ({"leftWidth": float("nan"), "rightWidth": 20, "rt": 15}, (-105, 135), None),
    ],
)
def test_feature_context_uses_detected_boundaries(feature, expected, boundaries):
    assert viewers._feature_range(feature) == (expected, boundaries)


def test_fingerprints_are_checked_only_when_evidence_is_opened(tmp_path, monkeypatch):
    from src.naseweis_provenance import acquisition_fingerprint

    source = tmp_path / "raw.mzML"
    source.write_bytes(b"original")
    params = {
        "mzML-files": str(source),
        "_acquisition_fingerprint": acquisition_fingerprint(source),
    }
    monkeypatch.setattr(
        viewers, "cached_fingerprint", lambda _: pytest.fail("Eager file hashing")
    )
    revision = viewers.saved_acquisition_revision(
        params, tmp_path, verify_content=False
    )
    assert revision[0] == str(source)
    monkeypatch.setattr(
        viewers, "cached_fingerprint", lambda _: acquisition_fingerprint(source)
    )
    assert viewers.saved_acquisition_revision(params, tmp_path) == (
        *revision,
        params["_acquisition_fingerprint"]["sha256"],
    )
    source.write_bytes(b"replaced")
    with pytest.raises(ValueError, match="has changed"):
        viewers.saved_acquisition_revision(params, tmp_path)


@pytest.mark.parametrize(
    "record",
    [
        {"sha256": "0" * 64, "size": True},
        {"sha256": "0" * 64, "size": -1},
        {"sha256": 123, "size": 0},
        {"sha256": "0" * 64, "size": "0"},
    ],
)
def test_invalid_fingerprint_fields_are_not_treated_as_legacy(tmp_path, record):
    path = tmp_path / "raw.mzML"
    path.touch()
    with pytest.raises(ValueError, match="fingerprint is invalid"):
        viewers.saved_acquisition_revision(
            {"mzML-files": str(path), "_acquisition_fingerprint": record}, tmp_path
        )


@pytest.mark.parametrize("recorded", [False, True])
def test_windows_saved_views_verify_content_without_reusing_a_stat_only_hash(
    tmp_path, monkeypatch, recorded
):
    import os

    from src import naseweis_provenance as provenance

    source = tmp_path / "raw.mzML"
    source.write_bytes(b"original")
    original_stat = source.stat()
    params = {"mzML-files": str(source)}
    expected = provenance.acquisition_fingerprint(source)
    if recorded:
        params["_acquisition_fingerprint"] = expected
    monkeypatch.setattr(provenance, "requires_content_verification", lambda: True)
    monkeypatch.setattr(
        viewers,
        "cached_fingerprint",
        lambda _: pytest.fail("Reused cached Windows hash"),
    )
    assert (
        viewers.saved_acquisition_revision(params, tmp_path)[-1] == expected["sha256"]
    )
    source.write_bytes(b"replaced")
    os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    if recorded:
        with pytest.raises(ValueError, match="has changed"):
            viewers.saved_acquisition_revision(params, tmp_path)
    else:
        assert (
            viewers.saved_acquisition_revision(params, tmp_path)[-1]
            != expected["sha256"]
        )
