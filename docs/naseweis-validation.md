# NASEWEIS modernization validation

Validated locally on 4 November 2026. No GitHub publication or remote deployment
was performed.

## Source revisions

- App: `poshul/nase-weis`, original `main` at `c1cc6f66fc0928295d40a7003937cb8876bc2f27`.
  `poshul/naseweis` was not accessible; the linked hyphenated repository was used.
- Template: `OpenMS/streamlit-template`, `main` at
  `e6e1333` (fetched during this update).
- Previous common template ancestor: `57fccde390ef0ce95a625c7a427840a8d69860c3`.
- Local development branch: `codex/modernize-naseweis`.

The template is merged with its history, rather than copied or squashed.

## Conflict resolutions

- Framework-owned `Dockerfile`, `src/common/common.py`, `src/fileupload.py`,
  `src/workflow/CommandExecutor.py` and `src/workflow/StreamlitUI.py`: template
  versions. NASEWEIS now uses the framework's upload/selection interfaces,
  automatic parameter persistence and worker process implementation.
- `.streamlit/config.toml`: template theme and server configuration, preserving
  NASEWEIS's 500 MB upload limit.
- `app.py`, `content/quickstart.py`, `src/Workflow.py`: preserve the RNA search
  purpose and rewrite the app-specific implementation against the new framework.
- `docker-compose.yml`: app-specific local service, loopback port, health check
  and named workspace volume; uses the pinned template runtime build described below.
- `test_gui.py`: tests the registered NASEWEIS pages and paired example loader.
- Modify/delete conflicts in `content/documentation.py`,
  `content/download_section.py`, `content/file_upload.py`,
  `content/raw_data_viewer.py`, `docs/build_app.md`, `docs/installation.md`,
  `docs/toppframework.py` and `docs/user_guide.md`: retain the current template
  sources as framework reference material. The extra pages are not registered
  in the app navigation.
- `settings.json` merged automatically, then was reconciled explicitly: add the
  new template keys; set the app/repository identity and branding; disable
  analytics, default demo workspaces, developer password and workspace sync.
  Existing online queue and sync implementations remain available for a hosted
  deployment configured separately.

## App changes

The app has Overview, Upload, Configure, Run and Results pages. The example
loader selects a paired mzML and RNA FASTA. Search settings expose the actual
mass tolerance and units, replacing undocumented resolution heuristics.
Every run keeps its own results and parameter snapshot. Failed input checks
and failed TOPP commands report failure to the worker.

Results are read from mzTab, avoiding the obsolete pyOpenMS list-to-dataframe
call. The match table derives q-values from the saved output, not the current
configuration. Coverage supports modified bases and terminal phosphates,
escapes HTML and wraps long references into 50-base blocks without fetching
MODOMICS data. CSV, mzTab, idXML, coverage HTML and run parameters can be exported.
The unused bespoke FASTA upload module was removed, and the generic spectrum
viewer follows the current template again.

## Automated validation

Run from the repository root with the isolated Python 3.12 environment:

```bash
.venv/bin/python -m utils.app_spec
.venv/bin/python -m pytest test.py test_gui.py tests/ -q -rs
.venv/bin/python -m pip check
.venv/bin/pylint $(git ls-files '*.py') --disable=C0103,C0114,C0301,C0411,W0212,W0631,W0602,W1514,W2402,E0401,E1101,F0001,R1732 --errors-only
git diff --check
```

The full suite covers app navigation/example setup, real NASE worker searches
with FDR on/off, missing inputs, TOPP failure propagation, run snapshots,
cache isolation, mzTab parsing, terminal/modified residues, coverage wrapping,
the container registry repair and the inherited framework tests. The
local-process test uses the bundled NME1 acquisition, exercising completion
after the parent process has already loaded an mzML file.

All 530 tests passed with no skips. RQ emits nine deprecation warnings for `job.exc_info` in
inherited queue tests. No broken Python requirements, app-spec errors, pylint
errors or whitespace errors were found.

## Container deployment

Built `naseweis:local` using `docker/Dockerfile.app`, then launched the Compose
service on `127.0.0.1:8501`. This host uses rootless Podman; the native Podman
Compose provider was selected because its Docker-compatible API socket is not
running. The application itself also ran natively on `127.0.0.1:8505`.

The local container build pins the template runtime's manifest:

```
ghcr.io/openms/streamlit-template@sha256:26640815c707bb685d57d197c591956cbc264a4736ab3ccee36c379aa486064e
```

That upstream image contains pyOpenMS/OpenMS 3.6 but omits the generated TOPP
registry: NASE aborted with “The TOPP tool registry ... is empty or unreadable”.
The app-layer build restores installed-tool registrations from the matching
`v3.6.0` source declarations and requires successful `-write_ini` calls for both
NASE and DecoyDatabase. The full source Dockerfiles now copy the generated
registry before deleting the OpenMS build tree. A full OpenMS source rebuild
and ARM/macOS/Windows/Kubernetes deployments were not run.

Real searches inside the deployed container produced:

| Settings | Spectrum matches | Distinct oligonucleotides |
| --- | ---: | ---: |
| RNase T1, 10 ppm precursor/fragment tolerance, FDR q ≤ 0.05 | 163 | 22 |
| Same search without FDR filtering | 299 | 31 |

The FDR run's q-values were all at most 0.05. Both modes produced nonempty
idXML and mzTab outputs and renderable coverage; the worker preserved the
workspace's current parameters while updating its own run snapshot.

Browser checks exercised the example loader, configuration, starting the
workflow, successful completion, results table, coverage and export controls.
The FDR example covers 65.0% of the 340-base NME1 reference. The healthy
container's `/_stcore/health` endpoint returned `ok`. A run created in the
browser was retained in the named workspace volume across container recreation.

## Scope and follow-up

The local deployment is for one operator on their machine. Hosted deployment
configuration remains separate. The Kubernetes manifests were specialized to
the app identity but Kubernetes is disabled in `app.yaml`. Desktop packaging
uses the NASEWEIS name and required tools; installers were not built.

For an existing installation, select its old workspace root explicitly with
`LOCAL_WORKSPACES_DIR` and review the new explicit mass tolerances. The example
counts demonstrate working software; they are not a statistical validation of
FDR estimation for small RNA databases.
