"""NASEWEIS workflow: optional RNA decoys, NASE search and result inspection."""

import shutil
from copy import deepcopy
from pathlib import Path

import streamlit as st
from streamlit.components.v1 import html

from src.naseweis_results import coverage_html, match_table, read_mztab
from src.workflow.ParameterManager import bool_param_paths_from_param_xml_ini
from src.workflow.WorkflowManager import WorkflowManager

SEARCH_DEFAULTS = {"oligo:enzyme": "RNase_T1"}
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
DECOY_FLAGS = ["only_decoy", "keep_terminal_nucleotides"]


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
        for label, key, extension in [
            ("MS/MS spectra", "mzML-files", "mzML"),
            ("RNA sequences", "fasta-files", "fasta"),
        ]:
            with st.expander(label, expanded=True):
                self.ui.upload_widget(key=key, name=label, file_types=[extension])

    def configure(self) -> None:
        self.ui.select_input_file(
            "mzML-files",
            name="MS/MS spectra",
            multiple=False,
            upload={"file_types": ["mzML"]},
        )
        self.ui.select_input_file(
            "fasta-files",
            name="RNA sequence database",
            multiple=False,
            upload={"file_types": ["fasta"]},
        )
        st.subheader("False discovery rate")
        self.ui.input_widget(
            "add-decoys",
            True,
            "Estimate FDR with shuffled RNA decoys",
            "Generate a target-decoy database and filter by q-value.",
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
            "Mass tolerances are explicit instrument settings. Search for an advanced parameter to reveal it."
        )
        self.ui.input_TOPP(
            "NucleicAcidSearchEngine",
            num_cols=3,
            custom_defaults=SEARCH_DEFAULTS,
            include_parameters=[
                "oligo:enzyme",
                "oligo:missed_cleavages",
                "oligo:min_size",
                "precursor:mass_tolerance",
                "precursor:mass_tolerance_unit",
                "fragment:mass_tolerance",
                "fragment:mass_tolerance_unit",
                "modifications:fixed",
                "modifications:variable",
            ],
            exclude_parameters=["fdr:"],
            flag_parameters=SEARCH_FLAGS,
        )
        with st.expander("Decoy generation settings"):
            self.ui.input_TOPP(
                "DecoyDatabase",
                num_cols=3,
                custom_defaults=DECOY_DEFAULTS,
                exclude_parameters=[
                    "type",
                    "decoy_string",
                    "decoy_string_position",
                    "only_decoy",
                    "enzyme",
                    "rnase",
                    "NeighborSearch:",
                    "Decoy:",
                ],
                flag_parameters=DECOY_FLAGS,
            )

    def execution(self) -> bool:
        """Run the tools and propagate failure; never save Streamlit state here."""
        for key, label in [
            ("mzML-files", "mzML spectra"),
            ("fasta-files", "RNA FASTA database"),
        ]:
            selected = self.params.get(key)
            # Reject ambiguous legacy selections instead of combining acquisitions.
            if not isinstance(selected, str) or not Path(selected).is_file():
                self.logger.log(
                    f"ERROR: Select one existing {label} file on Configure."
                )
                return False
        if not all(
            shutil.which(tool) for tool in ["DecoyDatabase", "NucleicAcidSearchEngine"]
        ):
            self.logger.log(
                "ERROR: Install OpenMS TOPP tools (DecoyDatabase and NucleicAcidSearchEngine) and add them to PATH."
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
        use_decoys = params.get("add-decoys", True)
        cutoff = float(params.get("FDR_cutoff", 0.05))
        if not 0 <= cutoff <= 1:
            self.logger.log("ERROR: Maximum q-value must be between 0 and 1.")
            return False
        # Override both enabled and disabled states on every run. In particular,
        # do not leave the previous run's decoy pattern or average-mass mode behind.
        search.update(
            {
                "fdr:cutoff": cutoff if use_decoys else 1.0,
                "fdr:decoy_pattern": "DECOY_" if use_decoys else "",
                "fdr:remove_decoys": bool(use_decoys),
            }
        )
        decoy.update(
            {
                "type": "RNA",
                "decoy_string": "DECOY_",
                "decoy_string_position": "prefix",
                "only_decoy": False,
            }
        )
        params["NucleicAcidSearchEngine"], params["DecoyDatabase"] = search, decoy
        params.setdefault("_flag_params", {})
        for tool, flags in [
            ("NucleicAcidSearchEngine", SEARCH_FLAGS),
            ("DecoyDatabase", DECOY_FLAGS),
        ]:
            if not self.parameter_manager.create_ini(tool):
                self.logger.log(f"ERROR: Cannot read {tool} parameters.")
                return False
            detected = bool_param_paths_from_param_xml_ini(
                self.parameter_manager.ini_dir / f"{tool}.ini", tool
            )
            params["_flag_params"][tool] = sorted(set(flags) | detected)
        # The manager points at this run's snapshot in both local and queue mode.
        self.parameter_manager.write_parameters(params)
        self.params = params
        database = [params["fasta-files"]]
        if use_decoys:
            output = self.file_manager.get_files(
                "target-decoy.fasta", set_results_dir="database"
            )
            self.logger.log("Generating shuffled RNA decoys...")
            if not self.executor.run_topp(
                "DecoyDatabase", {"in": [database], "out": output}
            ):
                return False
            database = output
        tab_out = self.file_manager.get_files(
            "search.mzTab", set_results_dir="mztab_results"
        )
        id_out = self.file_manager.get_files(
            "search.idXML", set_results_dir="idxml_results"
        )
        self.logger.log("Searching RNA oligonucleotides...")
        return self.executor.run_topp(
            "NucleicAcidSearchEngine",
            {
                "in": [params["mzML-files"]],
                "database": database,
                "out": tab_out,
                "id_out": id_out,
            },
        )

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
        table_tab, coverage_tab, download_tab = st.tabs(
            ["Matches", "Sequence coverage", "Export"]
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
