"""Standalone DecoyDatabase access, including its protein-specific capabilities."""

from copy import deepcopy
from pathlib import Path

import streamlit as st

from src.naseweis_tools import DECOY_INPUTS, DECOY_OUTPUTS, tool_flags, tool_schema
from src.workflow.WorkflowManager import WorkflowManager


class DecoyWorkflow(WorkflowManager):
    """Build databases independently of a NASE search and retain every result."""

    multi_run = True

    def __init__(self):
        super().__init__("Decoy Database", st.session_state["workspace"])

    def upload(self):
        self.ui.upload_widget(
            key="database-files",
            name="Sequence databases",
            file_types=["fasta", "fa", "faa"],
        )
        with st.expander("Optional relevant-protein database for NeighborSearch"):
            self.ui.upload_widget(
                key="neighbor-files",
                name="Relevant proteins",
                file_types=["fasta", "fa", "faa"],
            )

    def configure(self):
        advanced = st.toggle(
            "Advanced",
            key="advanced",
            help="Show decoy markers, shuffle controls, entrapment and NeighborSearch options.",
        )
        basic = ["type", "method"]
        excluded = []
        if not advanced and self.parameter_manager.create_ini("DecoyDatabase"):
            schema = tool_schema(self.parameter_manager.ini_dir / "DecoyDatabase.ini")
            excluded = sorted(schema.keys() - set(basic))
        self.ui.select_input_file(
            "database-files",
            name="Sequence databases",
            multiple=True,
            upload={"file_types": ["fasta", "fa", "faa"]},
        )
        self.ui.input_TOPP(
            "DecoyDatabase",
            num_cols=3,
            custom_defaults={"type": "RNA", "method": "shuffle"},
            include_parameters=basic,
            exclude_parameters=excluded,
            flag_parameters=tool_flags(
                self.parameter_manager.ini_dir / "DecoyDatabase.ini", ["only_decoy"]
            ),
        )
        if advanced:
            with st.expander("NeighborSearch input"):
                st.caption(
                    "Entrapment and NeighborSearch operate on proteins. Choose protein sequence type to use them."
                )
                self.ui.select_input_file(
                    "neighbor-files",
                    name="Relevant proteins (optional)",
                    multiple=False,
                    upload={"file_types": ["fasta", "fa", "faa"]},
                )

    def execution(self) -> bool:
        files = self.params.get("database-files", [])
        if (
            not isinstance(files, list)
            or not files
            or not all(isinstance(p, str) and Path(p).is_file() for p in files)
        ):
            self.logger.log("ERROR: Select at least one existing FASTA database.")
            return False
        if not self.parameter_manager.create_ini("DecoyDatabase"):
            self.logger.log(
                "ERROR: DecoyDatabase is not installed or its parameters cannot be read."
            )
            return False
        params = deepcopy(self.params)
        tool = (
            {"type": "RNA", "method": "shuffle"}
            | params.get("_defaults", {}).get("DecoyDatabase", {})
            | params.get("DecoyDatabase", {})
        )
        for key in DECOY_INPUTS | DECOY_OUTPUTS:
            tool.pop(key, None)
        neighbor = params.get("neighbor-files")
        if tool["type"] == "RNA" and (
            int(tool.get("entrapment_ratio", 0)) > 0 or neighbor
        ):
            self.logger.log(
                "ERROR: Entrapment and NeighborSearch require protein sequence type."
            )
            return False
        if neighbor and (not isinstance(neighbor, str) or not Path(neighbor).is_file()):
            self.logger.log("ERROR: The relevant-protein database no longer exists.")
            return False
        if (
            tool.get("method") != "shuffle"
            and int(tool.get("shuffle_decoy_ratio", 1)) != 1
        ):
            self.logger.log(
                "ERROR: Multiple decoys per target require the shuffle method."
            )
            return False
        schema = tool_schema(self.parameter_manager.ini_dir / "DecoyDatabase.ini")
        unknown = sorted(set(tool) - schema.keys())
        if unknown:
            self.logger.log(
                "ERROR: Unsupported DecoyDatabase parameters: " + ", ".join(unknown)
            )
            return False
        params["DecoyDatabase"] = tool
        params.setdefault("_flag_params", {})["DecoyDatabase"] = tool_flags(
            self.parameter_manager.ini_dir / "DecoyDatabase.ini", ["only_decoy"]
        )
        self.parameter_manager.write_parameters(params)
        io = {
            "in": [files],
            "out": self.file_manager.get_files(
                "database.fasta", set_results_dir="database"
            ),
        }
        if neighbor:
            io["NeighborSearch:in_relevant_proteins"] = [neighbor]
            for key, filename in [
                ("NeighborSearch:out_neighbor", "neighbors.fasta"),
                ("NeighborSearch:out_relevant", "relevant-targets-and-decoys.fasta"),
            ]:
                io[key] = self.file_manager.get_files(
                    filename, set_results_dir="database"
                )
        self.logger.log("Building sequence database...")
        return self.executor.run_topp("DecoyDatabase", io)

    def results(self):
        directory = self.results_dir / "database"
        files = sorted(directory.glob("*.fasta")) if directory.is_dir() else []
        if not files:
            st.info(
                "No generated databases yet. Upload sequences, configure the builder and start a run."
            )
            return
        for path in files:
            with path.open(encoding="utf-8") as handle:
                entries = sum(line.startswith(">") for line in handle)
            st.metric(path.name, f"{entries:,} sequences")
            st.download_button(
                f"Download {path.name}", path.read_bytes(), path.name, "text/plain"
            )
        snapshot = (
            self.runs.params_file(self.run_id)
            if self.run_id
            else self.parameter_manager.params_file
        )
        if snapshot.is_file():
            st.download_button(
                "Download database parameters",
                snapshot.read_bytes(),
                "decoy-database-parameters.json",
                "application/json",
            )
