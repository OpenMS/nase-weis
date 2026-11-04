"""Record reviewable NME1 reference results in a matching OpenMS 3.6 runtime.

Run explicitly in the pinned test container; never run automatically in CI.
Use a candidate directory, inspect the diff, then replace accepted fixtures.
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COLUMNS = {
    "NUC": ["accession", "sequence_coverage", "opt_sequence"],
    "OLI": ["sequence", "accession", "start", "end"],
    "OSM": [
        "sequence",
        "retention_time",
        "charge",
        "exp_mass_to_charge",
        "calc_mass_to_charge",
        "search_engine_score[1]",
        "search_engine_score[2]",
        "spectra_ref",
    ],
}


def main(output: Path) -> None:
    import pyopenms

    from src.workflow.RunManager import RunManager
    from src.workflow.tasks import execute_workflow

    help_text = subprocess.run(
        ["NucleicAcidSearchEngine", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    version = re.search(r"Version: (\S+)", help_text.stdout + help_text.stderr).group(1)
    if version != "3.6.0" or pyopenms.__version__ != "3.6.0":
        raise SystemExit(
            "Use the pinned OpenMS/pyOpenMS 3.6.0 container to record these references."
        )
    inputs = ["example-data/mzML/NME1.mzML", "example-data/fasta/NME1_Target.fasta"]
    manifest = {
        "openms_version": version,
        "pyopenms_version": pyopenms.__version__,
        "inputs": {
            path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in inputs
        },
        "settings": {
            "enzyme": "RNase_T1",
            "precursor_tolerance_ppm": 10.0,
            "fragment_tolerance_ppm": 10.0,
            "fdr_cutoff": 0.05,
            "threads": 1,
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="naseweis-reference-") as scratch:
        for use_decoys in [True, False]:
            name = "nme1-fdr" if use_decoys else "nme1-unfiltered"
            workflow = Path(scratch, name, "naseweis")
            workflow.mkdir(parents=True)
            source = workflow / "params.json"
            source.write_text(
                json.dumps(
                    {
                        "mzML-files": str(ROOT / inputs[0]),
                        "fasta-files": str(ROOT / inputs[1]),
                        "add-decoys": use_decoys,
                        "FDR_cutoff": 0.05,
                        "max_threads": 1,
                    }
                )
            )
            runs = RunManager(workflow)
            run_id = runs.new_run(name, source)
            execute_workflow(str(workflow), "Workflow", "src.Workflow", run_id)
            # Independent standard-library reader: a bug in the app's mzTab
            # parser must not silently become the accepted reference.
            headers, records = {}, {key: [] for key in COLUMNS}
            with (
                runs.results_dir(run_id) / "mztab_results/search.mzTab"
            ).open() as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if not row:
                        continue
                    if row[0] in {"NUH", "OLH", "OSH"}:
                        headers[{"NUH": "NUC", "OLH": "OLI", "OSH": "OSM"}[row[0]]] = (
                            row[1:]
                        )
                    elif row[0] in records:
                        fields = dict(zip(headers[row[0]], row[1:], strict=True))
                        records[row[0]].append(
                            [
                                fields[column]
                                for column in COLUMNS[row[0]]
                                if column in headers[row[0]]
                            ]
                        )
            for section, rows in records.items():
                with (output / f"{name}.{section.lower()}.tsv").open(
                    "w", newline=""
                ) as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(
                        [
                            column
                            for column in COLUMNS[section]
                            if column in headers[section]
                        ]
                    )
                    writer.writerows(sorted(rows))
            print(
                f"{name}: {len(records['OSM'])} spectrum matches, {len(records['OLI'])} mappings"
            )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
