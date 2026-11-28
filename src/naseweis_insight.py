"""Draw a ``src.naseweis_plots.Plot`` as its OpenMS-Insight component.

The component's cache lives in the workspace (``<workspace>/insight-cache``)
unless a folder is given, and is named by a hash of the rows and settings, so
a changed plot never reuses an old cache. Every plot offers its rows as TSV.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from src.naseweis_plots import Plot
from src.viz.plumbing import cache_name, component, only_accepted


def insight_cache(cache_dir=None) -> Path:
    """The folder the Insight caches of NASE-weis plots live in."""
    if cache_dir is not None:
        return Path(cache_dir)
    workspace = st.session_state.get("workspace")
    return Path(workspace, "insight-cache") if workspace else Path(".insight-cache")


def _draw(cls, path: Path, cache_id: str, build, call: dict) -> None:
    """Load or build the component and draw it; tests replace this to see what is drawn."""
    component(cls, path, cache_id, build)(**call)


def show_plot(plot: Plot, key: str, filename: str, cache_dir=None, state_manager=None, *,
              help: str | None = None, tsv: bool = True) -> None:  # pylint: disable=redefined-builtin
    """Render ``plot`` with its Insight component and, with ``tsv``, a TSV download of its rows."""
    import openms_insight

    cls = getattr(openms_insight, plot.component, None)
    if cls is None:
        st.info(f"This plot needs a newer openms-insight (it has no {plot.component} component).")
        return
    settings = only_accepted(cls, **plot.settings)
    cache_id = cache_name(filename, plot.data, sorted(settings.items(), key=lambda kv: kv[0]))
    path = insight_cache(cache_dir)
    path.mkdir(parents=True, exist_ok=True)
    call = {"key": key, "height": plot.height}
    if state_manager is not None:
        call["state_manager"] = state_manager
    if help:
        call["help"] = help
    _draw(cls, path, cache_id, lambda: {"data": plot.data.lazy(), **settings}, call)
    if tsv and not plot.table.empty:
        st.download_button(
            "Download plotted values (TSV)",
            plot.table.to_csv(sep="\t", index=False),
            f"{filename}.tsv",
            "text/tab-separated-values",
            key=f"{key}-tsv",
        )


def _column(name: str, dtype, titles: dict) -> dict:
    import polars as pl

    column = {"field": name, "title": titles.get(name, name)}
    if dtype.is_numeric() and dtype != pl.Boolean:
        column["sorter"] = "number"
        if dtype.is_float():
            column.update(formatter="money", formatterParams={"thousand": False, "precision": 4})
    return column


def show_table(frame, key: str, filename: str, *, title: str | None = None, titles: dict | None = None,
               cache_dir=None, height: int = 376, state_manager=None, csv: str | None = None) -> None:
    """Draw ``frame`` (pandas or polars) as an Insight ``Table``, with an optional CSV download."""
    import openms_insight
    import pandas as pd
    import polars as pl

    shown = frame
    if isinstance(frame, pd.DataFrame):
        shown = frame.reset_index(drop=True)
        shown.columns = [str(c) for c in shown.columns]
        shown = shown.astype({c: "string" for c in shown.columns if shown[c].dtype == object})
    data = pl.from_pandas(shown) if isinstance(shown, pd.DataFrame) else shown
    data = data.with_columns(
        [pl.col(c).cast(pl.Utf8) for c, t in data.schema.items() if t == pl.Object or t.is_nested()]
    )
    index = "_row"
    while index in data.columns:
        index += "_"
    data = data.with_row_index(index).with_columns(pl.col(index).cast(pl.Int64))
    columns = [_column(c, t, titles or {}) for c, t in data.schema.items() if c != index]
    settings = {"column_definitions": columns, "index_field": index, "title": title}
    cache_id = cache_name(filename, data, sorted(((k, str(v)) for k, v in settings.items()), key=lambda kv: kv[0]))
    path = insight_cache(cache_dir)
    path.mkdir(parents=True, exist_ok=True)
    call = {"key": key, "height": height}
    if state_manager is not None:
        call["state_manager"] = state_manager
    _draw(openms_insight.Table, path, cache_id, lambda: {"data": data.lazy(), **settings}, call)
    if csv:
        st.download_button(f"Download {csv} (CSV)", data.drop(index).write_csv().encode(),
                           f"{filename}.csv", "text/csv", key=f"{key}-csv")
