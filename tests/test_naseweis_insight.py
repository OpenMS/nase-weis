"""How NASE-weis hands its plots and tables to OpenMS-Insight."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl
import pytest

import src.naseweis_insight as insight
from src.naseweis_plots import Plot


class _Streamlit:
    def __init__(self, workspace=None):
        self.session_state = {} if workspace is None else {"workspace": str(workspace)}
        self.downloads = []

    def download_button(self, label, data, file_name, mime, key):
        self.downloads.append({"label": label, "data": data, "file": file_name, "mime": mime, "key": key})


class _Calls(list):
    st: _Streamlit


@pytest.fixture
def drawn(monkeypatch, tmp_path):
    calls = _Calls()
    fake = _Streamlit(tmp_path)

    def draw(cls, path, cache_id, build, call):
        calls.append({"cls": cls, "path": path, "cache_id": cache_id, "built": build(), "call": call})

    monkeypatch.setattr(insight, "_draw", draw)
    monkeypatch.setattr(insight, "st", fake)
    calls.st = fake
    return calls


def test_cache_folder_follows_workspace(monkeypatch, tmp_path):
    monkeypatch.setattr(insight, "st", _Streamlit(tmp_path))
    assert insight.insight_cache() == tmp_path / "insight-cache"
    assert insight.insight_cache(tmp_path / "x") == tmp_path / "x"
    monkeypatch.setattr(insight, "st", _Streamlit())
    assert insight.insight_cache() == Path(".insight-cache")


def test_draw_loads_or_builds_and_calls_the_component(monkeypatch, tmp_path):
    seen = {}

    def component(cls, path, cache_id, build):
        seen.update(cls=cls, path=path, cache_id=cache_id, build=build)
        return lambda **call: seen.update(call=call)

    monkeypatch.setattr(insight, "component", component)
    insight._draw("Cls", tmp_path, "id", "build", {"key": "k"})
    assert seen == {"cls": "Cls", "path": tmp_path, "cache_id": "id", "build": "build", "call": {"key": "k"}}


def _plot(values=(1.0, 2.0), table=None):
    data = pl.DataFrame({"x": list(values), "y": [3.0, 4.0][: len(values)]})
    settings = {"x_column": "x", "y_column": "y", "title": "T", "no_such_option": 1}
    return Plot("LinePlot", data, settings, table if table is not None else data.to_pandas(), height=300)


def test_plot_is_drawn_with_its_rows_and_a_tsv_download(drawn, tmp_path):
    insight.show_plot(_plot(), "p", "scores", state_manager="state")
    (call,) = drawn
    assert call["cls"].__name__ == "LinePlot"
    assert call["path"] == tmp_path / "insight-cache" and call["path"].is_dir()
    assert call["call"] == {"key": "p", "height": 300, "state_manager": "state"}
    assert call["built"]["data"].collect()["x"].to_list() == [1.0, 2.0]
    assert call["built"]["x_column"] == "x"
    # Options the installed openms-insight does not take are left out.
    assert "no_such_option" not in call["built"]
    (download,) = drawn.st.downloads
    assert download["file"] == "scores.tsv" and download["key"] == "p-tsv"
    assert download["data"].splitlines()[0] == "x\ty"


def test_changed_rows_get_another_cache_and_empty_table_no_download(drawn):
    insight.show_plot(_plot(table=pd.DataFrame()), "p", "scores", cache_dir=None)
    insight.show_plot(_plot(values=(5.0, 6.0), table=pd.DataFrame()), "p", "scores")
    first, second = drawn
    assert first["cache_id"] != second["cache_id"]
    assert "state_manager" not in first["call"]
    assert drawn.st.downloads == []


def test_table_columns_index_and_csv(drawn, tmp_path):
    frame = pd.DataFrame(
        {"_row": ["a", "b"], "score": [0.5, 1.25], "count": [1, 2], "ok": [True, False],
         "nested": [[1, 2], [3]], 7: ["s", None]},
        index=[10, 20],
    )
    insight.show_table(frame, "t", "matches", title="Matches", titles={"score": "Score"},
                       cache_dir=tmp_path / "c", height=200, state_manager="state", csv="matches")
    (call,) = drawn
    assert call["cls"].__name__ == "Table"
    assert call["path"] == tmp_path / "c"
    assert call["call"] == {"key": "t", "height": 200, "state_manager": "state"}
    built = call["built"]
    # The row index never collides with a column of the table.
    assert built["index_field"] == "_row_"
    assert built["title"] == "Matches"
    columns = {c["field"]: c for c in built["column_definitions"]}
    assert list(columns) == ["_row", "score", "count", "ok", "nested", "7"]
    assert columns["score"]["title"] == "Score" and columns["score"]["formatter"] == "money"
    assert columns["count"]["sorter"] == "number" and "formatter" not in columns["count"]
    assert "sorter" not in columns["ok"] and "sorter" not in columns["_row"]
    rows = built["data"].collect()
    assert rows["_row_"].to_list() == [0, 1]
    assert rows.schema["nested"] == pl.Utf8
    (download,) = drawn.st.downloads
    assert download["label"] == "Download matches (CSV)" and download["file"] == "matches.csv"
    assert download["data"].decode().splitlines()[0] == "_row,score,count,ok,nested,7"


def test_polars_table_without_csv(drawn):
    insight.show_table(pl.DataFrame({"a": [1]}), "t", "plain")
    (call,) = drawn
    assert call["built"]["index_field"] == "_row"
    assert call["call"] == {"key": "t", "height": 376}
    assert drawn.st.downloads == []
