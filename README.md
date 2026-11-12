# NASEWEIS

Search RNA oligonucleotides in MS/MS data with the OpenMS
**NucleicAcidSearchEngine**. Select MS/MS acquisitions and an RNA FASTA database or cached digest,
configure digestion, modifications and mass tolerances, then inspect spectrum
matches and reference sequence coverage. Optional label-free quantification follows
the RNA TOPPAS chain: NASE targets → FeatureFinderMetaboIdent on the original
acquisition → TextExporter. LFQ is off by default. When enabled for multiple
selected samples, it produces a combined modified-RNA abundance matrix,
optional total-area normalization and relative condition comparisons. Native
feature areas, featureXML, metadata and viewer evidence are retained per sample.

Built on [OpenMS/streamlit-template](https://github.com/OpenMS/streamlit-template),
with Streamlit 1.54, pinned **OpenMS nightly 3.7** command-line tools and pyOpenMS **3.7.0.dev20261110** from the official nightly wheel index. The framework is merged from `main`, so
future template updates retain their Git history.

## Run locally with containers

Docker and Docker Compose are sufficient; OpenMS and Python are included.

```bash
git clone https://github.com/poshul/nase-weis.git
cd nase-weis
docker compose up -d --build
```

Open [http://localhost:8501](http://localhost:8501). The service binds to
loopback and stores workspaces in the `naseweis-workspaces` named volume.
`docker compose down` stops the app without deleting its data.

Put your input files in the **`data` folder beside `docker-compose.yml`**. You can
move or copy mzML, FASTA, cached digest and other supported inputs there, including
subfolders. In the app's file browser, choose **Data folder**, select
**Refresh folder** after adding files, then select the files for your workflow.
This folder is mounted read-only at `/data`; results stay in the workspace volume.
Keep the selected inputs in place so saved runs can reopen their original data in
the spectrum and LFQ viewers. Files in `data` are excluded from Git and container
builds.

To use an existing host folder instead, create a local `.env` file beside
`docker-compose.yml` containing its absolute path:

```dotenv
NASEWEIS_DATA_DIR=/absolute/path/to/ms-data
```

Create the folder if needed, then run `docker compose up -d` to apply the mount.
The app shows this host folder in its file browser. The local `.env` file is ignored
by Git and container builds.

Absolute host paths in the app are clickable `file://` folder links. Browsers
may block these links or show a directory in the browser instead of opening
the OS file manager. Use the copy button beside the displayed path and paste
it into your file manager's address or Go to Folder field when needed. The
default `./data` path stays as copyable text relative to `docker-compose.yml`.

```bash
docker compose ps
docker compose logs --tail=100 naseweis
docker compose down
```

Compose uses `docker/Dockerfile.app`: the complete native runtime comes from
an immutable official OpenMS executables image; Python/UI dependencies come
from the pinned template image. Native tools, libraries, data and registry stay
together. A self-contained official pyOpenMS nightly wheel provides Python parsing;
searches and database generation use the **3.7 nightly executables**.
`openms-runtime.json` records the image digest and source commit
`07075086398bc0c2d8d18e47111ed43fbecfc9e1`. New upstream nightlies do not
silently change existing deployments.

The source-build Dockerfiles and desktop default select the same source
commit. Simple images now include the same official nightly native runtime.
See [nightly features and validation](docs/naseweis-nightly.md) for details.

On a host using Podman without a Docker-compatible API socket:

```bash
podman build -f docker/Dockerfile.app -t naseweis:local .
PODMAN_COMPOSE_PROVIDER=/usr/bin/podman-compose docker compose up -d --no-build
```

## Run with Python

Use Python 3.11 or 3.12 and install the pinned OpenMS nightly TOPP tools separately. Both
`NucleicAcidSearchEngine` and `DecoyDatabase` must be on `PATH`; optional LFQ also
requires `FeatureFinderMetaboIdent` and `TextExporter`. Installing the
pyOpenMS wheel alone does not install command-line tools.

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py --server.address 127.0.0.1
```

Native workspaces default to `../workspaces-nase-weis`. Set
`LOCAL_WORKSPACES_DIR=/absolute/path` to choose another location, including
an existing workspace root from an older installation. Input data and results
are kept outside the checkout. The local app has analytics disabled and does
not require Redis or a CAPTCHA.

## Example search

1. Select **Load example dataset** on Overview. Both NME1 input files are
   copied into the workspace and selected automatically.
2. On **Configure**, review the RNA database, RNase, modifications and mass
   tolerances. Enable **Advanced** for MVH scoring, FLASH MS1 precursor correction,
   MS2 deisotoping, filters, mass-accuracy weighting, reporting and export options.
   Settings are saved automatically; hiding Advanced preserves their saved values.
   Select the acquisitions to analyze and assign their conditions and replicates.
   On the **Quantification** tab, optionally enable LFQ for at least two acquisitions containing
   MS1 scans. Its extraction settings are under Advanced and default to the saved
   Herbert TOPPAS workflow. See [LFQ settings and provenance](docs/lfq.md).
3. On **Run**, optionally enter a run name, then select **Start Workflow**.
4. On **Results**, choose a completed run. Inspect matches and sequence
   coverage, plus the **Quantification** tab for LFQ results; download the selected
   native artifacts, CSV, coverage HTML and search parameters.

The example contains 512 MS1 and 2,116 MS2 spectra. With the tested OpenMS
3.7 nightly container and default RNase T1 / 10 ppm settings, the example yields 163
matches at q ≤ 0.05 and 301 without FDR filtering (3.6 yielded 299). Exact scores and counts can
change with OpenMS versions. The example is a software demonstration, not a
calibration of FDR estimation for a small RNA database.

FASTA sequences use A, U, C and G, with bracketed MODOMICS short names for
modified residues, for example `AU[Am]C[m5C]G`. The nightly runtime also accepts Thermo `.raw` and Bruker `.d` inputs.
Upload `.d` directories as ZIP archives or import them from a local folder.
Vendor-reader availability depends on the platform; the bundled mzML example
is the tested scientific reference. Each run has separate results, logs and a parameter
snapshot; a new run does not overwrite earlier results. Legacy NASEWEIS
result directories are adopted by the framework's run manager.

Mass tolerances are now explicit, including ppm/Da units. The previous
resolution-based heuristics are removed: review these settings when reusing
an old workspace. A q-value column appears only when the result file actually
contains q-values, independently of the current configuration.

## Interactive viewers

The **Raw data viewer** uses the existing pyOpenMS-viz/Plotly viewers for MS1
peak maps, individual MS1/MS2 spectra and overlaid TIC, BPC and extracted-ion
chromatograms. Choose an uploaded or locally registered mzML file, then load it.
Overview intensities are summed into display bins; spectra and XICs retain the
original peak values. Indexed mzML files are read one spectrum at a time, and
unindexed files use streaming reads.

On **Results → Spectra**, choose a saved spectrum match to inspect its original
scan and the RNA fragment assignments in NASE's idXML. Complete native spectrum
IDs identify scans; retention time and array position are not substituted.
RNA `a-B`, `d` and `w` ions are supported. Processed fragment masses absent from
the acquired scan remain available in the assignment table.

On **Results → Quantification**, open the LFQ feature viewer to inspect the
feature map and an exact raw MS1 XIC for a selected feature, with its detected
elution boundaries. View settings do not change saved quantification. All plots
support zoom, hover, image export and self-contained HTML downloads.

Result viewers use the selected run's input snapshot, independently of current
search settings. Keep its original mzML available. Vendor acquisitions remain
supported for searches; convert them to mzML for these viewers. See
[viewer behavior and validation](docs/viewers.md).

## Development and verification

Run the complete regression gate before merging future improvements:

```bash
python3 scripts/check_regressions.py --container
```

It checks the current checkout in the pinned OpenMS runtime, compares all NME1
scientific results with reviewed references, rejects skipped checks and enforces
95% application coverage. It leaves a running local app and its workspaces alone.
See [the regression maintenance guide](docs/naseweis-regressions.md) for faster
local checks, reports, CI and intentional baseline updates.

```bash
python -m pip install -r requirements-dev.txt
python -m utils.app_spec
python -m pytest test.py test_gui.py tests/
python -m pytest tests/test_naseweis_integration.py -v
```

The integration tests exercise every native export, FLASH, MVH, cached digests,
custom decoy markers, provided target-decoy inputs, protein entrapment and
NeighborSearch, plus searches with FDR enabled and disabled through
the same worker entry point as the app. They require TOPP tools and skip with
an explicit reason when these are unavailable. The remaining tests cover page
loading, example setup, tool failures, modified sequences, coverage, result
cache isolation, framework parameters, run management and workspace sync.
Some framework tests start localhost servers; run them with socket access.

App-specific code lives in `src/Workflow.py`, `src/DecoyWorkflow.py`,
`src/naseweis_tools.py`, `src/naseweis_results.py`, `src/naseweis_lfq.py` and
the `src/naseweis_*` viewer modules.
Navigation and user-facing pages live in `app.py` and `content/`.
`app.yaml` records the app's identity and file ownership for template syncs.
The template's extra toolbox pages remain available in source but are not
shown in the NASEWEIS navigation. Kubernetes manifests are inherited and
specialized for this app; Kubernetes deployment is disabled in `app.yaml`
and was not deployed as part of local validation.

See [the modernization validation record](docs/naseweis-validation.md) for
the template revision, conflict resolutions and tests performed.

## Citation

- [NucleicAcidSearchEngine: RNA modification mapping by mass spectrometry](https://doi.org/10.1038/s41467-020-14665-7).
- Müller, T. D., Siraj, A., et al. *OpenMS WebApps: Building User-Friendly
  Solutions for MS Analysis.* Journal of Proteome Research (2025).
  [DOI](https://doi.org/10.1021/acs.jproteome.4c00872).
- Pfeuffer, J., Bielow, C., Wein, S. et al. *OpenMS 3 enables reproducible
  analysis of large-scale mass spectrometry data.* Nature Methods 21,
  365–367 (2024). [DOI](https://doi.org/10.1038/s41592-024-02197-7).
