# Nightly NASE and DecoyDatabase support

The app targets the **OpenMS `nightly` branch**, pinned to
[`553de474c71b4f9926e5fcb264530e3d30b72b88`](https://github.com/OpenMS/OpenMS/commit/553de474c71b4f9926e5fcb264530e3d30b72b88).
The executable reports `3.7.0-pre-HEAD-2026-11-03`, revision `553de47`.
`nightly` and `develop` are separate branches; parameters found only on a
newer `develop` are not added here.

## Search capabilities

Configure shows common parameters and provides searchable **All parameters**
dialogs generated from the installed executables' ParamXML. This includes:

- FLASH precursor correction, its mass/cosine/SNR controls and audit output;
  legacy/none/FLASH RNA deisotoping, precursor removal, window mowing and
  largest-peak filtering.
- Hyperscore or MVH scoring, intensity classes, TIC fraction, single-peak
  credit, mass-accuracy weighting and m/z blacklisting.
- Native charge/isotope/adduct options, fragment ions, fixed/variable
  modifications and ambiguity resolution, report ties/top hits, minimum
  fragment coverage and precursor-intensity reporting.
- Independent decoy generation and FDR estimation, configurable markers,
  target-only or target/decoy result retention, and supplied target-decoy inputs.
- FASTA or cached OMS digest input. A cached digest retains its original
  digestion and decoy provenance; its scientific settings must match the search.
- bedRMod reference metadata and an optional modification mapping CSV with
  `mod`/`name`, `chebi_id` and/or `modomics_short_name` columns.
- Native mzML, Thermo RAW and Bruker directory input. The upload/import
  framework accepts zipped `.d` directories. Vendor data were not included
  in the bundled regression dataset; native reader support is platform dependent.

mzTab is always retained. **Additional result files** selects idXML, an OMS
identification database, a reusable OMS digest, LFQ targets, bedRMod,
theoretical/experimental MS2 spectra, preprocessed spectra and FLASH audit TSV.
Every selected file is stored inside its run and offered on Results. Native
optional mzTab diagnostics remain visible in the match table and CSV.

The app rejects incompatible FLASH mass settings, missing a-B ions for
ambiguity resolution, digest re-export, unknown native settings and unavailable
exports. It discovers CLI flags from the runtime schema; `single_peak_credit`
remains a true/false **value**, as required by NASE.

## Standalone database builder

The **Decoy Database** page exposes RNA and protein database preparation
independently of a search, with multiple FASTA inputs, reverse/shuffle, seed,
identity threshold, retry limit, decoy marker position, decoy-only output and
multiple shuffle decoys per target. Protein mode includes enzyme-aware
shuffling, entrapment proteins and NeighborSearch, including its relevant and
neighbor FASTA downloads. Each build retains its parameters and its own results.

NASE's integrated builder always uses RNA. Protein-only modes are rejected
there. Generated decoy-only files are retained separately, and the search uses
a combined target/decoy FASTA so targets are still present. NASE's FDR code has
no decoy-ratio correction parameter in this revision: FDR with generated decoys
requires ratio 1; other ratios remain available in the standalone builder or
when FDR is disabled. Supplied databases/digests must have suitable provenance.
DecoyDatabase enforces its native recognized prefix/suffix rules for markers.

## Runtime and packaging

The locally validated deployment is `docker/Dockerfile.app`, used by Compose
and the strict CI regression gate. It takes native tools **and all their
libraries/data** from the official OpenMS image digest in `openms-runtime.json`.
It copies the template Python environment and installs the self-contained
PyPI pyOpenMS 3.6 wheel for ParamXML/mzML parsing. This avoids loading the
template's locally linked 3.6 wheel against the nightly native library.

The source-build Dockerfiles select the same OpenMS source commit and align
libxml2/Arrow and native Thermo-reader dependencies with upstream. Simple
images use the official nightly native runtime. Desktop builds select the
same source commit. These alternative images and desktop installers require
their platform CI builds; local validation covers the Compose image on Linux
amd64. Updating the pin is an explicit change, accompanied by fresh schemas,
scientific comparisons and the complete regression gate.

## Scientific comparison

Both versions use the unchanged NME1 mzML/FASTA, RNase T1, 10 ppm precursor
and fragment tolerances, one thread and default hyperscore.

| Runtime | Matches at q ≤ 0.05 | Unfiltered matches | Filtered/unfiltered mappings |
| --- | ---: | ---: | ---: |
| OpenMS 3.6.0 | 163 | 299 | 22 / 32 |
| Nightly `553de47` | 163 | 301 | 22 / 32 |

All 163 filtered assignments and all 299 original unfiltered assignments
remain. Among shared assignments, 33 filtered scores and 23 q-values change;
55 unfiltered scores change. Calculated m/z values remain unchanged. The two
additional unfiltered matches are `UAAAAUUUUUUUUACUCAGp` at scan 4263, charge
−8 (score 2.468479), and scan 5166, charge −9 (score 1.713535). Neither passes
the default FDR filter. Their cause has not been isolated to a particular
upstream change; the versioned reference records the observed behavior.

The 3.6 references are preserved under `tests/data/naseweis/`; nightly schemas
and references live under `tests/data/naseweis/nightly/`. CI compares every
assignment, score, q-value and mapping with the nightly reference and rejects
skips. Real-tool feature tests cover all outputs, FLASH/MVH, cached digest
equivalence, marker/decoy-only modes, supplied target-decoy input, protein
entrapment and NeighborSearch. See [regression maintenance](naseweis-regressions.md).


Local validation on 2026-11-04: **650 tests passed with no skips**, **99.57%**
combined application statement/branch coverage; both workflow classes and all
registered pages reached **100%**. The strict gate also passed the real
container entrypoint health check and `pip check`. Browser validation completed
an example search with six native export types and verified their download
controls. A browser-discovered stale file selector in the standalone builder
was replaced by rerendering preparation steps and given a regression test.
