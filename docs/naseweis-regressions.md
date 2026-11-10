# NASEWEIS regression gate

Run this before merging a NASEWEIS change:

```bash
python3 scripts/check_regressions.py --container
```

This builds the **current checkout** on the immutable OpenMS runtime in
`docker/Dockerfile.app`, then builds a separate test image. The test container
has networking disabled and mounts only its report directory. It never stops
or recreates the running app, publishes ports or mounts users' workspaces.
The production image's real entrypoint is also started in a temporary container
and its health endpoint must return `ok`; that container is always removed.
Docker is the default; use `--engine podman` for native Podman.

The check requires every test to pass, rejects skips and expected failures,
and enforces at least **95% combined statement/branch coverage** of the
app entry point, registered pages and NASEWEIS analysis/result modules.
A test checks that new pages and app-owned analysis modules are not left out
of the coverage configuration. Framework tests run too, but thousands of
framework lines do not dilute the application coverage measurement.

Reports are written to `.local/regressions/junit.xml` and `coverage.xml`.
The main `continuous-integration` workflow calls the reusable `NASEWEIS
regressions` workflow for pull requests and pushes to `main`. The reusable
workflow can also be started manually and uploads reports even when tests fail.
It has no separate push/PR trigger, so each CI run starts the pinned gate once.
Once published, its emitted **NASEWEIS regressions / NASEWEIS pinned OpenMS
regression gate** status check can be made required in branch protection. This local change does not
modify repository settings or publish anything to GitHub.

## What the gate protects

- Real searches through the worker entry point, with FDR enabled and disabled;
  every spectrum assignment, RNA/oligo sequence, charge, retention time,
  experimental/calculated m/z, score, q-value and mapping coordinate is compared
  with the pinned scientific reference.
- All native NASE exports; cached digests producing identical assignments;
  FLASH precursor correction and deisotoping, MVH scoring and report filters;
  independent FDR/decoy generation, custom prefixes/suffixes, decoy-only export;
  RNA/protein shuffling, ratios, protein entrapment and NeighborSearch outputs.
- Invalid/missing/ambiguous input selections, invalid FDR cutoffs, missing
  TOPP tools, unreadable tool parameters and failing pipeline stages.
- Explicit tolerance units, modifications, digestion, parameter precedence,
  RNA-only decoys, FDR toggling and isolation of each run's parameter snapshot.
- Cold start, every registered page, example loading, persistence of changed
  configuration, selecting old/new runs and reading legacy results.
- Empty, missing and malformed results; optional export files; the distinction
  between saved FDR results and the current configuration; coverage and export
  controls staying usable when one artifact is missing or damaged.
- Modified/terminal residues, HTML escaping, coverage coordinates/wrapping,
  result cache refresh and isolation between workspaces, runs and readers.
- The container's native runtime and inherited framework behavior,
  including queued/local execution, cancellation, storage and workspace sync.

Reference comparator tests deliberately remove a match or alter its charge,
score, q-value, reference or coordinates and require each change to be rejected.
These checks demonstrate that a still-successful search cannot silently lose
or reassign results and pass the reference gate.

## Faster local checks

Use the project virtual environment with TOPP tools on `PATH`:

```bash
.venv/bin/python -m pip install -r requirements-dev.txt
.venv/bin/python scripts/check_regressions.py
```

This runs all tests and the coverage gate on your installed tools. It checks
workflow invariants but does not claim scientific equivalence to a different
OpenMS build. To run only the new behavioral tests while editing:

```bash
.venv/bin/python -m pytest test_gui.py tests/test_naseweis*.py -q
```

The ordinary pytest invocation permits missing-tool skips for convenience;
the strict regression gate fails instead. For a full reproducible check,
use `--container`. `--naseweis-baseline` requires the exact nightly version
and source revision recorded in the manifest, plus the pinned pyOpenMS version.

## Scientific references and intentional changes

`tests/data/naseweis/nightly/manifest.json` records the engine/pyOpenMS versions,
source revision and native image digest,
input checksums and baseline settings. The six TSV files record the nucleotide,
oligonucleotide and spectrum-match sections of two NME1 searches. They are
readable review artifacts, not opaque hashes. Numeric comparisons allow
`1e-9` relative/absolute rounding; assignments, row counts and sequence/position
changes must match. Row ordering and runtime paths/timestamps are ignored.

For pinned OpenMS nightly `0707508`, default RNase T1 searches at 10 ppm and one thread produce:

| Mode | Spectrum matches | Oligonucleotide mappings |
| --- | ---: | ---: |
| FDR q ≤ 0.05 | 163 | 22 |
| No FDR | 301 | 32 |

The unfiltered mappings include 31 distinct oligonucleotide sequences. The
original OpenMS 3.6 references remain in the parent directory for comparison.
See [the nightly comparison](naseweis-nightly.md) for changed scores and the
two additional unfiltered matches.
The reference generator parses the raw engine output using Python's CSV
reader, independently of the app's parser. Reference generation is **never**
part of CI or a normal regression run.

To propose an intentional scientific/runtime change:

1. Explain the expected change and update the runtime pin/version guard if
   upgrading OpenMS. Retain an explanation of any changed input data.
2. Build the test image with `scripts/check_regressions.py --container`.
   A changed runtime or result is expected to fail against the old baseline.
3. Generate candidates into a separate directory, using the matching test
   image; for example:

   ```bash
   mkdir -p .local/reference-candidate
   docker run --rm --network none \
     -v "$PWD/.local/reference-candidate:/candidate" \
     --entrypoint /root/miniforge3/envs/streamlit-env/bin/python \
     naseweis-regression-test:local \
     -m scripts.generate_naseweis_reference --output /candidate
   diff -ru tests/data/naseweis/nightly .local/reference-candidate
   ```

4. Review the changed matches, scores, q-values and coverage, rather than just
   accepting a new row count. Replace the fixtures only after accepting the
   scientific change, and describe it in the change being reviewed.
5. Run the full pinned gate again. Do not lower the coverage threshold or
   regenerate references simply to hide an unexpected failure.

These tests enforce known contracts and reviewed scientific outputs. They
cannot establish that every possible regression or scientific error is absent;
new features need new behavior tests and suitable scientific examples.

## Validation of this gate

On 4 November 2026, both the native Python 3.12 run and pinned Python 3.11 /
OpenMS 3.6 container run passed **586 tests with no skips**, with **99.47%**
combined application statement/branch coverage. The production-entrypoint
smoke test passed. An intentionally skipped probe exited with failure, and
the baseline guard rejected the host's development OpenMS version.
Existing RQ tests emit nine deprecation warnings; they do not affect the gate.


The nightly upgrade's container gate passed **650 tests with no skips** and
**99.57%** application coverage. Both search and database workflows and every
registered page had 100% coverage. See the [nightly validation record](naseweis-nightly.md).
