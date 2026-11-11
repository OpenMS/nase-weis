"""NASEWEIS workflow: optional RNA decoys, NASE search and result inspection."""

import json
import shutil
from copy import deepcopy
from pathlib import Path

import streamlit as st
from streamlit.components.v1 import html

from src.naseweis_lfq import (
    FEATURE_FINDER_DEFAULTS,
    LFQ_FEATURES,
    LFQ_TABLE,
    LFQ_TOOLS,
    prepare_lfq,
    read_lfq_table,
    run_lfq,
)
from src.naseweis_results import coverage_html, match_table, read_mztab
from src.naseweis_tools import (
    DECOY_INPUTS,
    DECOY_OUTPUTS,
    SEARCH_EXPORTS,
    SEARCH_INPUTS,
    enabled,
    tool_flags,
    tool_schema,
    validate_search_options,
)
from src.naseweis_viewers import render_identification_viewer, render_lfq_viewer
from src.workflow.WorkflowManager import WorkflowManager

SEARCH_DEFAULTS = {
    "oligo:enzyme": "RNase_T1",
    "fdr:decoy_pattern": "DECOY_",
    "fdr:remove_decoys": True,
}
BASIC_SEARCH_PARAMETERS = [
    "oligo:enzyme",
    "oligo:missed_cleavages",
    "oligo:min_size",
    "precursor:mass_tolerance",
    "precursor:mass_tolerance_unit",
    "fragment:mass_tolerance",
    "fragment:mass_tolerance_unit",
    "modifications:fixed",
    "modifications:variable",
    "modifications:variable_max_per_oligo",
]
DECOY_DEFAULTS = {"type": "RNA", "method": "shuffle"}
# Older TOPP versions store CLI flags as strings in ParamXML, so keep the
# known flags as well as those detected from current type="bool" entries.
SEARCH_FLAGS = [
    "precursor:include_unknown_charge",
    "precursor:use_avg_mass",
    "precursor:use_adducts",
    "decharge_ms2",
    "modifications:resolve_ambiguities",
    "fdr:remove_decoys",
]
DECOY_FLAGS = ["only_decoy"]


def validate_lfq_spectra(path: str) -> str | None:
    """Check mzML scan metadata without decoding peak arrays into memory.

    Vendor formats use the native reader in FeatureFinderMetaboIdent, which
    reports reader availability and missing MS1 scans in the run log.
    """
    if Path(path).suffix.lower() != ".mzml":
        return None
    import pyopenms as oms

    try:
        options = oms.PeakFileOptions()
        options.setFillData(False)
        reader = oms.MzMLFile()
        reader.setOptions(options)
        experiment = oms.MSExperiment()
        reader.load(path, experiment)
        if 1 not in experiment.getMSLevels():
            return "LFQ requires MS1 scans in the original acquisition. This mzML has no MS1 scans; disable quantification or select an acquisition containing MS1 and MS2."
    except (RuntimeError, OSError) as error:
        return f"Cannot inspect the LFQ acquisition: {error}"
    return None


@st.cache_data(show_spinner=False)
def load_results(path: str, mtime_ns: int, size: int) -> dict:
    """Cache by full result path and file revision, including the workspace/run."""
    return read_mztab(path)


class Workflow(WorkflowManager):
    """Keep each search and its parameter snapshot as an independent run."""

    multi_run = True

    def __init__(self) -> None:
        super().__init__("NASEWEIS", st.session_state["workspace"])

    def load_example_data(self) -> None:
        """Copy and select the paired NME1 dataset without replacing user files."""
        params = self.parameter_manager.read_parameters_strict()
        for key, folder, filename in [
            ("mzML-files", "mzML", "NME1.mzML"),
            ("fasta-files", "fasta", "NME1_Target.fasta"),
        ]:
            destination = self.workflow_dir / "input-files" / key / filename
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                shutil.copyfile(Path("example-data", folder, filename), destination)
            params[key] = str(destination.resolve())
        self.parameter_manager.write_parameters(params)
        self.parameter_manager.clear_parameter_session_state()

    def upload(self) -> None:
        if st.button("Load example dataset", type="primary", icon=":material/science:"):
            self.load_example_data()
            st.success("NME1 spectra and RNA database loaded. Continue to Configure.")
        for label, key, extensions in [
            ("MS/MS spectra", "mzML-files", ["mzML", "raw", "d"]),
            ("RNA sequences", "fasta-files", ["fasta", "fa", "faa"]),
            ("Reusable database digest (optional)", "digest-files", ["oms"]),
            (
                "bedRMod modification mapping (optional)",
                "bedrmod-mapping-files",
                ["csv"],
            ),
        ]:
            with st.expander(label, expanded=key in {"mzML-files", "fasta-files"}):
                self.ui.upload_widget(key=key, name=label, file_types=extensions)

    def configure(self) -> None:
        advanced = st.toggle(
            "Advanced",
            key="advanced",
            help="Show search, database and LFQ extraction options. Hiding controls keeps their saved values.",
        )
        identification, quantification = st.tabs(["Identification", "Quantification"])
        with identification:
            self._configure_identification(advanced)
        with quantification:
            self._configure_quantification(advanced)

    def _configure_identification(self, advanced: bool) -> None:
        self.ui.select_input_file(
            "mzML-files",
            name="MS/MS spectra",
            multiple=False,
            upload={"file_types": ["mzML", "raw", "d"]},
        )
        if advanced:
            self.ui.input_widget(
                "database-source",
                "FASTA",
                "Search database source",
                "A cached digest keeps its original digestion and decoy provenance.",
                "selectbox",
                options=["FASTA", "Cached digest"],
                reactive=True,
            )
        source = st.session_state.get(
            f"{self.parameter_manager.param_prefix}database-source",
            self.params.get("database-source", "FASTA"),
        )
        if source == "Cached digest":
            self.ui.select_input_file(
                "digest-files",
                name="Cached RNA digest",
                multiple=False,
                upload={"file_types": ["oms"]},
            )
            st.info(
                "The cached database controls digestion. Enable FDR only if its original database includes decoys."
            )
        else:
            self.ui.select_input_file(
                "fasta-files",
                name="RNA sequence database",
                multiple=False,
                upload={"file_types": ["fasta", "fa", "faa"]},
            )
        st.subheader("False discovery rate")
        self.ui.input_widget(
            "add-decoys",
            True,
            "Generate a target-decoy RNA database",
            "Disable for a FASTA that already includes decoys, or when searching a cached digest.",
        )
        if advanced:
            self.ui.input_widget(
                "estimate-fdr",
                self.params.get("add-decoys", True),
                "Estimate false discovery rate",
                "Requires a target-decoy FASTA or a cached digest containing decoys.",
            )
        self.ui.input_widget(
            "FDR_cutoff",
            0.05,
            "Maximum q-value",
            "0.05 corresponds to 5% estimated FDR.",
            "number",
            min_value=0.0,
            max_value=1.0,
            step_size=0.01,
        )
        st.subheader("Search settings")
        st.caption(
            "Mass tolerances are explicit instrument settings. Enable Advanced for additional search controls."
        )
        # Filter the schema as well as inline widgets so the All parameters
        # dialog cannot bypass Advanced. Saved overrides remain untouched.
        excluded = ["fdr:cutoff"]
        if not advanced and self.parameter_manager.create_ini(
            "NucleicAcidSearchEngine"
        ):
            schema = tool_schema(
                self.parameter_manager.ini_dir / "NucleicAcidSearchEngine.ini"
            )
            excluded += sorted(schema.keys() - set(BASIC_SEARCH_PARAMETERS))
        self.ui.input_TOPP(
            "NucleicAcidSearchEngine",
            num_cols=3,
            custom_defaults=SEARCH_DEFAULTS,
            include_parameters=BASIC_SEARCH_PARAMETERS,
            exclude_parameters=excluded,
            flag_parameters=tool_flags(
                self.parameter_manager.ini_dir / "NucleicAcidSearchEngine.ini",
                SEARCH_FLAGS,
            ),
        )
        if advanced:
            with st.expander("Decoy generation settings"):
                self.ui.input_TOPP(
                    "DecoyDatabase",
                    num_cols=3,
                    custom_defaults=DECOY_DEFAULTS,
                    exclude_parameters=[
                        "type",
                        "enzyme",
                        "rnase",
                        "NeighborSearch:",
                        "entrapment_ratio",
                        "entrapment_string",
                        "Decoy:",
                    ],
                    flag_parameters=tool_flags(
                        self.parameter_manager.ini_dir / "DecoyDatabase.ini",
                        DECOY_FLAGS,
                    ),
                )

            st.subheader("Additional result files")
            export_key = f"{self.parameter_manager.param_prefix}search-exports"
            if export_key not in st.session_state:
                st.session_state[export_key] = self.params.get(
                    "search-exports", ["id_out"]
                )
            st.multiselect(
                "Result files",
                list(SEARCH_EXPORTS),
                key=export_key,
                format_func=lambda option: SEARCH_EXPORTS[option][0],
                help="Selected files are retained with this run and can be downloaded from Results.",
                on_change=self.parameter_manager.save_parameters,
            )
            self.parameter_manager.save_parameters()
            self.ui.select_input_file(
                "bedrmod-mapping-files",
                name="Custom modification mapping (optional)",
                multiple=False,
                upload={"file_types": ["csv"]},
            )
            st.caption(
                "RNA mode excludes protein-only entrapment and neighbor searches. The standalone DecoyDatabase page exposes the complete tool."
            )

    def _configure_quantification(self, advanced: bool) -> None:
        st.subheader("Label-free quantification")
        st.write(
            "Extract MS1 features for the oligonucleotides identified in this acquisition, following the RNA TOPPAS workflow."
        )
        self.ui.input_widget(
            "enable-lfq",
            False,
            "Enable label-free quantification",
            "Runs FeatureFinderMetaboIdent and TextExporter after a successful NASE search. Requires MS1 scans in the original acquisition.",
            "checkbox",
            reactive=True,
        )
        lfq_enabled = enabled(st.session_state.get(
            f"{self.parameter_manager.param_prefix}enable-lfq",
            self.params.get("enable-lfq", False),
        ))
        with st.container(border=True):
            st.markdown("**NASE targets → MS1 feature extraction → feature table**")
            st.caption(
                "Uses the original acquisition and the identified oligo formulas, charges and retention times. Saves the target list, featureXML and a TSV with all feature metadata."
            )
        st.info(
            "Each run produces raw feature intensities for one acquisition. This TOPPAS workflow does not align, link or normalize separate runs."
        )
        if not lfq_enabled:
            st.caption("Enable quantification to add these steps to the next run.")
            return
        st.caption(
            "Defaults follow pipeline2_fixed.toppas: full 50 ppm m/z and 240 s RT windows, isotope probability ≥ 0.05, and no elution-model fitting. Enable Advanced to adjust extraction and detection."
        )
        if advanced:
            self.ui.input_TOPP(
                "FeatureFinderMetaboIdent",
                num_cols=3,
                custom_defaults=FEATURE_FINDER_DEFAULTS,
                include_parameters=[
                    "extract:mz_window", "extract:rt_window", "extract:n_isotopes",
                    "extract:isotope_pmin", "detect:peak_width",
                    "detect:min_peak_width", "detect:signal_to_noise", "model:type",
                ],
                flag_parameters=tool_flags(
                    self.parameter_manager.ini_dir / "FeatureFinderMetaboIdent.ini"
                ),
            )
            st.caption(
                "m/z windows ≥ 1 use ppm; values below 1 use Da. A positive isotope probability takes precedence over the isotope count. Automatic RT extraction needs a positive peak width. Model fitting options apply only when a model is selected."
            )

    def execution(self) -> bool:
        """Run the tools and propagate failure; never save Streamlit state here."""
        from src.naseweis_provenance import (
            acquisition_fingerprint,
            assert_acquisition_unchanged,
        )

        cached_digest = self.params.get("database-source", "FASTA") == "Cached digest"
        for key, label in [
            ("mzML-files", "MS spectra"),
            ("digest-files", "cached RNA digest")
            if cached_digest
            else ("fasta-files", "RNA FASTA database"),
        ]:
            selected = self.params.get(key)
            valid_vendor = (
                key == "mzML-files"
                and isinstance(selected, str)
                and Path(selected).suffix.lower() == ".d"
                and Path(selected).is_dir()
            )
            if not isinstance(selected, str) or not (
                Path(selected).is_file() or valid_vendor
            ):
                self.logger.log(
                    f"ERROR: Select one existing {label} file on Configure."
                )
                return False
        needed = ["NucleicAcidSearchEngine"]
        lfq_enabled = enabled(self.params.get("enable-lfq", False))
        if lfq_enabled:
            needed.extend(LFQ_TOOLS)
        if enabled(self.params.get("add-decoys", True)) and not cached_digest:
            needed.append("DecoyDatabase")
        if not all(shutil.which(tool) for tool in needed):
            self.logger.log(
                "ERROR: Install the required OpenMS TOPP tools and add them to PATH: " + ", ".join(needed)
            )
            return False

        params = deepcopy(self.params)
        search = (
            SEARCH_DEFAULTS
            | params.get("_defaults", {}).get("NucleicAcidSearchEngine", {})
            | params.get("NucleicAcidSearchEngine", {})
        )
        decoy = (
            DECOY_DEFAULTS
            | params.get("_defaults", {}).get("DecoyDatabase", {})
            | params.get("DecoyDatabase", {})
        )
        # The executor merges stored defaults again at dispatch. Clear file
        # options from both layers so only the workspace selectors own paths.
        # Native log paths are strings in ParamXML; the executor captures logs.
        for layer in (
            search,
            params.get("_defaults", {}).get("NucleicAcidSearchEngine", {}),
        ):
            for key in SEARCH_INPUTS | set(SEARCH_EXPORTS) | {"out", "log"}:
                layer.pop(key, None)
        for layer in (decoy, params.get("_defaults", {}).get("DecoyDatabase", {})):
            for key in DECOY_INPUTS | DECOY_OUTPUTS | {"log"}:
                layer.pop(key, None)
        use_decoys = enabled(params.get("add-decoys", True)) and not cached_digest
        calculate_fdr = enabled(
            params.get("estimate-fdr", params.get("add-decoys", True))
        )
        exports = list(params.get("search-exports", ["id_out"]))
        if lfq_enabled and "lfq_out" not in exports:
            exports.append("lfq_out")
        error = validate_search_options(search, exports, cached_digest=cached_digest)
        if error:
            self.logger.log(f"ERROR: {error}")
            return False
        if lfq_enabled:
            error = prepare_lfq(params, self.parameter_manager, params["mzML-files"])
            if error is None:
                error = validate_lfq_spectra(params["mzML-files"])
            if error:
                self.logger.log(f"ERROR: {error}")
                return False
        cutoff = float(params.get("FDR_cutoff", 0.05))
        if not 0 <= cutoff <= 1:
            self.logger.log("ERROR: Maximum q-value must be between 0 and 1.")
            return False
        # Preserve native settings; derive only the flags owned by the workflow.
        search["fdr:cutoff"] = cutoff if calculate_fdr else 1.0
        search["fdr:decoy_pattern"] = (
            (
                decoy.get("decoy_string", "DECOY_")
                if use_decoys
                else search.get("fdr:decoy_pattern", "DECOY_")
            )
            if calculate_fdr
            else ""
        )
        search["fdr:remove_decoys"] = (
            enabled(search.get("fdr:remove_decoys", True)) if calculate_fdr else False
        )
        # NASE consumes RNA sequences, while the standalone builder also supports proteins.
        decoy["type"] = "RNA"
        decoy.setdefault("decoy_string", "DECOY_")
        decoy.setdefault("decoy_string_position", "prefix")
        if use_decoys and (not str(decoy["decoy_string"]).strip()):
            self.logger.log("ERROR: The decoy accession marker must not be empty.")
            return False
        if calculate_fdr and not str(search["fdr:decoy_pattern"]).strip():
            self.logger.log(
                "ERROR: FDR estimation requires a nonempty decoy accession marker."
            )
            return False
        if (
            use_decoys
            and calculate_fdr
            and int(decoy.get("shuffle_decoy_ratio", 1)) != 1
        ):
            self.logger.log(
                "ERROR: NASE's FDR estimation requires a 1:1 target-decoy ratio. Use ratio 1 or disable FDR."
            )
            return False
        params["NucleicAcidSearchEngine"], params["DecoyDatabase"] = search, decoy
        params.setdefault("_flag_params", {})
        tools = [("NucleicAcidSearchEngine", SEARCH_FLAGS)]
        if use_decoys:
            tools.append(("DecoyDatabase", DECOY_FLAGS))
        for tool, flags in tools:
            if not self.parameter_manager.create_ini(tool):
                self.logger.log(f"ERROR: Cannot read {tool} parameters.")
                return False
            params["_flag_params"][tool] = tool_flags(
                self.parameter_manager.ini_dir / f"{tool}.ini", flags
            )
        schema = tool_schema(
            self.parameter_manager.ini_dir / "NucleicAcidSearchEngine.ini"
        )
        unknown = sorted(set(search) - schema.keys())
        if use_decoys:
            decoy_schema = tool_schema(
                self.parameter_manager.ini_dir / "DecoyDatabase.ini"
            )
            unknown += [
                f"DecoyDatabase:{key}"
                for key in sorted(set(decoy) - decoy_schema.keys())
            ]
        if unknown:
            self.logger.log(
                "ERROR: Unsupported parameters for these OpenMS tools: "
                + ", ".join(unknown)
            )
            return False
        for option in exports:
            if option not in schema:
                self.logger.log(
                    f"ERROR: This NASE installation cannot generate {SEARCH_EXPORTS[option][0]}; use the pinned nightly runtime."
                )
                return False
        # The manager points at this run's snapshot in both local and queue mode.
        # Imported snapshots cannot claim provenance for the selected input.
        params.pop("_acquisition_fingerprint", None)
        acquisition_record = None
        if Path(params["mzML-files"]).suffix.lower() == ".mzml":
            try:
                acquisition_record = acquisition_fingerprint(params["mzML-files"])
            except (OSError, ValueError) as error:
                self.logger.log(f"ERROR: {error}")
                return False
            params["_acquisition_fingerprint"] = acquisition_record
        self.parameter_manager.write_parameters(params)
        self.params = params
        database = [params["digest-files"] if cached_digest else params["fasta-files"]]
        if use_decoys:
            output = self.file_manager.get_files(
                "decoys.fasta"
                if enabled(decoy.get("only_decoy", False))
                else "target-decoy.fasta",
                set_results_dir="database",
            )
            self.logger.log("Generating RNA decoys...")
            if not self.executor.run_topp(
                "DecoyDatabase", {"in": [database], "out": output}
            ):
                return False
            if enabled(decoy.get("only_decoy", False)):
                combined = self.results_dir / "database" / "target-and-decoy.fasta"
                combined.write_bytes(
                    Path(params["fasta-files"]).read_bytes().rstrip()
                    + b"\n"
                    + Path(output[0]).read_bytes()
                )
                database = [str(combined)]
            else:
                database = output
        tab_out = self.file_manager.get_files(
            "search.mzTab", set_results_dir="mztab_results"
        )
        io = {
            "in": [params["mzML-files"]],
            "digest" if cached_digest else "database": database,
            "out": tab_out,
        }
        mapping = params.get("bedrmod-mapping-files")
        if mapping:
            if not Path(mapping).is_file():
                self.logger.log(
                    "ERROR: The selected custom modification mapping no longer exists."
                )
                return False
            io["bedrmod_chebi_mapping"] = [mapping]
        for option in exports:
            relative = SEARCH_EXPORTS[option][1]
            destination = self.results_dir / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            io[option] = [str(destination)]
        self.logger.log("Searching RNA oligonucleotides...")
        if not self.executor.run_topp("NucleicAcidSearchEngine", io):
            return False
        if acquisition_record is not None:
            try:
                assert_acquisition_unchanged(params["mzML-files"], acquisition_record)
            except ValueError as error:
                self.logger.log(f"ERROR: {error}")
                return False
        if lfq_enabled:
            if not run_lfq(self, io["lfq_out"][0]):
                return False
            if acquisition_record is not None:
                try:
                    assert_acquisition_unchanged(params["mzML-files"], acquisition_record)
                except ValueError as error:
                    self.logger.log(f"ERROR: {error}")
                    return False
        return True

    def results(self) -> None:
        # Legacy runs used extensionless names; retain access after migration.
        result_file = self.results_dir / "mztab_results" / "search.mzTab"
        if not result_file.exists():
            result_file = self.results_dir / "mztab_results" / "tab_out"
        if not result_file.exists():
            st.info("No search results yet. Configure the inputs and start a run.")
            return
        try:
            stat = result_file.stat()
            tables = load_results(
                str(result_file.resolve()), stat.st_mtime_ns, stat.st_size
            )
            matches = match_table(tables)
        except (OSError, ValueError, KeyError) as error:
            st.error(f"Could not read this run's results: {error}")
            return
        cols = st.columns(3)
        cols[0].metric("Spectrum matches", len(matches))
        cols[1].metric(
            "Distinct oligonucleotides",
            matches["Sequence"].nunique() if not matches.empty else 0,
        )
        cols[2].metric("RNA accessions", len(tables.get("NUC", [])))
        snapshot = (
            self.runs.params_file(self.run_id)
            if self.run_id
            else self.parameter_manager.params_file
        )
        online = st.session_state.get("settings", {}).get("online_deployment", False)
        table_tab, spectra_tab, coverage_tab, lfq_tab, download_tab = st.tabs(
            ["Matches", "Spectra", "Sequence coverage", "Quantification", "Export"]
        )
        with table_tab:
            if matches.empty:
                st.info(
                    "No matches passed this search's filters. Review the database, digestion, mass tolerances and FDR cutoff."
                )
            else:
                st.dataframe(matches, hide_index=True, use_container_width=True)
                st.download_button(
                    "Download matches (CSV)",
                    matches.to_csv(index=False).encode(),
                    "naseweis-matches.csv",
                    "text/csv",
                )
        with spectra_tab:
            render_identification_viewer(
                self.results_dir, snapshot, tables.get("OSM", matches),
                self.workflow_dir, online,
            )
        with coverage_tab:
            references = tables.get("NUC")
            if references is not None and not references.empty:
                accession = st.selectbox(
                    "RNA accession", references["accession"].dropna().unique()
                )
                try:
                    plot = coverage_html(tables, accession)
                    html(plot, height=550, scrolling=True)
                    st.download_button(
                        "Download coverage (HTML)",
                        plot,
                        "naseweis-coverage.html",
                        "text/html",
                    )
                except (ValueError, KeyError, IndexError) as error:
                    st.warning(f"Coverage is unavailable: {error}")
            else:
                st.info("No RNA sequence coverage to display for this run.")
        with lfq_tab:
            self.quantification_results()
        with download_tab:
            st.caption(
                "mzTab contains all result tables; idXML opens in TOPPView. Parameters record the settings used for this run."
            )
            for label, path, mime in [
                ("Download mzTab", result_file, "text/tab-separated-values"),
                (
                    "Download idXML",
                    self.results_dir / "idxml_results" / "search.idXML",
                    "application/xml",
                ),
                (
                    "Download legacy idXML",
                    self.results_dir / "idxml_results" / "id_out",
                    "application/xml",
                ),
            ]:
                if path.is_file():
                    st.download_button(label, path.read_bytes(), path.name, mime)
            for option, (label, relative, mime) in SEARCH_EXPORTS.items():
                if option == "id_out":
                    continue
                artifact = self.results_dir / relative
                if artifact.is_file():
                    st.download_button(
                        f"Download {label}", artifact.read_bytes(), artifact.name, mime
                    )
            database_dir = self.results_dir / "database"
            if database_dir.is_dir():
                for artifact in sorted(database_dir.glob("*.fasta")):
                    st.download_button(
                        f"Download {artifact.name}",
                        artifact.read_bytes(),
                        artifact.name,
                        "text/plain",
                    )
            snapshot = (
                self.runs.params_file(self.run_id)
                if self.run_id
                else self.parameter_manager.params_file
            )
            if snapshot.is_file():
                st.download_button(
                    "Download search parameters",
                    snapshot.read_bytes(),
                    "naseweis-parameters.json",
                    "application/json",
                )
            log_path = (
                self.runs.run_dir(self.run_id) if self.run_id else self.workflow_dir
            ) / "logs" / "all.log"
            if log_path.is_file():
                st.download_button("Download run log", log_path.read_bytes(), "naseweis-run.log", "text/plain")

    def select_run(self, label: str = "Run", key: str = "") -> str | None:
        """Keep identification artifacts accessible if downstream LFQ fails."""
        from src.workflow.run_widgets import run_label, selected_run_key

        available = [
            run for run in self.runs.list_runs()
            if run["status"] == "finished" or (
                run["status"] == "failed" and any(
                    (self.runs.results_dir(run["id"]) / "mztab_results" / name).is_file()
                    for name in ("search.mzTab", "tab_out")
                )
            )
        ]
        if not available:
            return None
        choices = {run["id"]: run for run in available}
        shared = selected_run_key(self.runs)
        widget = f"runs-{self.workflow_dir.name}-select-widget{key}"
        selected = st.session_state.get(shared)
        st.session_state[widget] = selected if selected in choices else available[0]["id"]

        def store_selection():
            st.session_state[shared] = st.session_state[widget]

        run_id = st.selectbox(
            label, list(choices), key=widget,
            format_func=lambda value: run_label(choices[value]) + (
                " · Failed / partial results" if choices[value]["status"] == "failed" else ""
            ),
            on_change=store_selection,
            help="Choose a saved run. Failed runs with search artifacts remain available for inspecting partial outputs.",
        )
        store_selection()
        self._use_run(run_id)
        if choices[run_id]["status"] == "failed":
            st.warning("This run failed. Available outputs may be incomplete; inspect this run's saved log on Export before interpreting them.")
        return run_id

    def quantification_results(self) -> None:
        """Render the selected run's artifacts independently of live settings."""
        st.subheader("MS1 feature intensities")
        st.caption(
            "Identification-guided extraction from the original acquisition. Intensities are native feature areas; separate runs have not been aligned, linked or normalized. Multiple charge states may produce separate features."
        )
        table_path = self.results_dir / LFQ_TABLE
        features_path = self.results_dir / LFQ_FEATURES
        if table_path.is_file():
            try:
                table = read_lfq_table(table_path)
                if table.empty:
                    st.info("No MS1 features were extracted for this run's identified targets. Review the target list, extraction windows and MS1 signal.")
                else:
                    cols = st.columns(2)
                    cols[0].metric("Extracted features", len(table))
                    if "label" in table:
                        cols[1].metric("Oligonucleotide labels", table["label"].nunique())
                    leading = [column for column in ("label", "charge", "intensity", "rt", "mz", "expected_rt", "sum_formula") if column in table]
                    st.dataframe(
                        table[leading + [column for column in table if column not in leading]],
                        hide_index=True, use_container_width=True,
                        column_config={
                            "label": "Oligonucleotide", "charge": "Charge",
                            "intensity": st.column_config.NumberColumn("Raw peak area", format="%.3g"),
                            "rt": st.column_config.NumberColumn("RT (s)", format="%.2f"),
                            "mz": st.column_config.NumberColumn("m/z", format="%.5f"),
                            "expected_rt": st.column_config.NumberColumn("Expected RT (s)", format="%.2f"),
                        },
                    )
                    snapshot = self.runs.params_file(self.run_id) if self.run_id else self.parameter_manager.params_file
                    render_lfq_viewer(
                        self.results_dir, snapshot, table, self.workflow_dir,
                        st.session_state.get("settings", {}).get("online_deployment", False),
                    )
            except (OSError, ValueError) as error:
                st.error(f"Could not read this run's feature table: {error}")
        else:
            snapshot = self.runs.params_file(self.run_id) if self.run_id else self.parameter_manager.params_file
            try:
                saved = json.loads(snapshot.read_text())
                if not isinstance(saved, dict):
                    saved = {}
            except (OSError, ValueError):
                saved = {}
            if enabled(saved.get("enable-lfq", False)):
                st.warning("Quantification has no completed feature table for this run. Check the run status and saved log on Export; identification outputs remain available.")
            else:
                st.info("Quantification was not enabled for this run. Enable it on Configure → Quantification and start a new run.")
        for label, path, mime in [
            ("Download LFQ features (featureXML)", features_path, "application/xml"),
            ("Download LFQ feature table (TSV)", table_path, "text/tab-separated-values"),
        ]:
            if path.is_file():
                st.download_button(label, path.read_bytes(), path.name, mime)
