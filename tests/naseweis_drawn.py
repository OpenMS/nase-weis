"""What NASE-weis draws as Insight components, for AppTest, which cannot see them.

``record(monkeypatch)`` replaces the one place every NASE-weis plot and table
is drawn and returns the list it fills: one ``Drawn`` per component, with its
class name, key, rows (pandas, without the table's own row index) and
settings.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Drawn:
    component: str
    key: str
    data: object
    settings: dict
    toolbar: list | None = None

    @property
    def title(self):
        return self.settings.get("title")


def record(monkeypatch) -> list[Drawn]:
    import src.naseweis_insight as insight

    drawn: list[Drawn] = []
    run = [None]

    def draw(cls, _path, _cache_id, build, call):
        # A rerun inside one AppTest run (st.rerun) starts the list again.
        from streamlit.runtime.scriptrunner import get_script_run_ctx

        ctx = get_script_run_ctx()
        this_run = getattr(ctx, "widget_user_keys_this_run", None)
        if this_run is not run[0]:
            run[0] = this_run
            drawn.clear()
        built = build()
        data = built.pop("data").collect().to_pandas()
        index = built.get("index_field")
        if index in data:
            data = data.drop(columns=index)
        drawn.append(Drawn(cls.__name__, call["key"], data, built, call.get("toolbar_settings")))

    monkeypatch.setattr(insight, "_draw", draw)
    return drawn


def tables(drawn: list[Drawn]):
    """The rows of every Insight table drawn, in order."""
    return [item.data for item in drawn if item.component == "Table"]


def plots(drawn: list[Drawn]):
    return [item for item in drawn if item.component != "Table"]


def attach(app, monkeypatch) -> list[Drawn]:
    """Record what ``app`` draws; each run of its script starts a new list."""
    drawn = record(monkeypatch)
    run = app._run

    def fresh_run(*args, **kwargs):
        drawn.clear()
        return run(*args, **kwargs)

    app._run = fresh_run
    app.drawn = drawn
    return drawn


def table_frames(app):
    """The rows of every Insight table the last run of ``app`` drew."""
    return tables(app.drawn)


def by_key(drawn: list[Drawn], part: str) -> list[Drawn]:
    return [item for item in drawn if part in item.key]


def coverage_shown(app) -> bool:
    """Whether the last run of ``app`` showed a coverage map: drawn by Insight
    (whose own menu downloads it), or the static map with its SVG download."""
    if any(item.component == "SequenceCoverage" for item in app.drawn):
        return True
    return any(item.proto.label == "Download coverage (SVG)" for item in app.get("download_button"))
