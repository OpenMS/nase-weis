# NASEWEIS

Search RNA oligonucleotides in MS/MS data with the OpenMS
**NucleicAcidSearchEngine**. Upload an MS/MS run and an RNA FASTA database or cached digest,
configure digestion, modifications and mass tolerances, then inspect spectrum
matches and reference sequence coverage.

Built on [OpenMS/streamlit-template](https://github.com/OpenMS/streamlit-template),
with Streamlit 1.54, pinned **OpenMS nightly 3.7** command-line tools and pyOpenMS 3.6 for the UI ParamXML API. The framework is merged from `main`, so
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

```bash
docker compose ps
docker compose logs --tail=100 naseweis
docker compose down
```

Compose uses `docker/Dockerfile.app`: the complete native runtime comes from
an immutable official OpenMS executables image; Python/UI dependencies come
from the pinned template image. Native tools, libraries, data and registry stay
together. A self-contained PyPI pyOpenMS 3.6 wheel provides Python parsing;
searches and database generation use the **3.7 nightly executables**.
`openms-runtime.json` records the image digest and source commit
`553de474c71b4f9926e5fcb264530e3d30b72b88`. New upstream nightlies do not
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
`NucleicAcidSearchEngine` and `DecoyDatabase` must be on `PATH`; installing the
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
   tolerances. Settings are saved automatically.
3. On **Run**, optionally enter a run name, then select **Start Workflow**.
4. On **Results**, choose a completed run. Inspect matches and sequence
   coverage; download the selected native artifacts, CSV, coverage HTML and search parameters.

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
`src/naseweis_tools.py` and `src/naseweis_results.py`.
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
