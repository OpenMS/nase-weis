# Optional RNA label-free quantification

On **Configure → Quantification**, enable optional LFQ before starting a run.
Use **Advanced** to adjust extraction settings. LFQ is off by default. It uses the same input run and the
identifications retained by the NASE search, including its selected FDR and
reporting settings. Merely exporting LFQ targets does not enable quantification.

The workflow follows the existing RNA TOPPAS graph:

```text
RNA FASTA → optional DecoyDatabase → NucleicAcidSearchEngine
original mzML ────────────────────→ NucleicAcidSearchEngine
                                   │ lfq_out (target TSV)
original mzML ──────────────────────┼─→ FeatureFinderMetaboIdent
                                   │       │ featureXML
                                   │       └─→ TextExporter → feature TSV
```

FeatureFinderMetaboIdent extracts isotope chromatograms from the original
**MS1** data. NASE MS2 deisotoping and MS1 precursor correction settings affect
identification; their exported processed MS2 spectra are not the LFQ input.
An mzML without MS1 scans is rejected when LFQ is enabled. Vendor inputs follow
the native OpenMS readers and must also contain MS1 data.

## Results and interpretation

Each run retains these files:

| File | Contents |
| --- | --- |
| `exports/lfq-targets.tsv` | Native NASE targets: modified RNA sequence/adduct label, formula, observed charges, and expected retention time |
| `quantification/features.featureXML` | Native detected features, intensities, metadata and quality annotations; viewable in TOPPView |
| `quantification/features.tsv` | TextExporter feature table, including all feature metadata |

The TSV includes `label`, `sum_formula`, `charge`, `expected_rt`, measured `mz`,
retention time and `intensity`. Its first column is normally named `#rt`.
The default intensities are unnormalized feature peak areas in instrument
units. They are not absolute concentrations, and separate charge-state
features are not automatically summed into an oligonucleotide abundance.
Changing the elution model can change how feature intensity is estimated.

A target can have no detected feature, and an identification search can produce
no LFQ targets. An empty target list completes with an empty featureXML and TSV;
it does not synthesize zero-valued abundance observations. Identification
results and targets remain available for inspection. Native errors in
extraction or export cause the run to fail and are recorded in its logs.

The present workflow quantifies one input run at a time. It performs no
cross-run retention-time alignment, feature linking, normalization, missing-value
imputation or statistical comparisons. The native model option
`model:no_imputation` concerns failed elution-model fits within a run; it does
not provide cross-run missing-value imputation.

## TOPPAS reference settings

The starting preset preserves the analytical settings stored in the Herbert
`pipeline2_fixed.toppas` workflow, saved by OpenMS
`3.0.0-pre-HEAD-2022-09-16`. Its SHA256 is
`1cecf741a05e22123d3471c841910097329552b3f662064c1b8a5077f2141000`.
No private experimental data or local filesystem paths are packaged with the
application. The current executable and pyOpenMS versions are pinned in
[`openms-runtime.json`](../openms-runtime.json).

| Setting | Preset | Meaning |
| --- | --- | --- |
| `extract:mz_window` | `50.0` | Full extraction window, ppm when ≥1; Da/Th when <1 |
| `extract:rt_window` | `240.0` | Full extraction window in seconds; zero uses four times the peak width |
| `extract:n_isotopes` | `2` | Fixed isotope count when probability filtering is disabled |
| `extract:isotope_pmin` | `0.05` | Minimum isotope probability; a positive value takes precedence over the count |
| `detect:peak_width` | `0.0` | Expected elution peak width in seconds for smoothing |
| `detect:min_peak_width` | `0.2` | Relative to peak width below 1; seconds otherwise |
| `detect:signal_to_noise` | `0.8` | Feature-detection threshold |
| `model:type` | `none` | Do not fit an elution model |
| `EMGScoring:max_iteration` | `100` | Maximum EMG scoring iterations |
| `EMGScoring:init_mom` | `false` | Preserve the reference's scoring initialization |

Other stored model defaults are `add_zeros=0.2`, `unweighted_fit=false`,
`no_imputation=false`, `each_trace=false`, and checks `min_area=1`,
`boundaries=0.5`, `width=10`, `asymmetry=10`. They are relevant when an elution
model is selected. Parameters and their resolved values are saved with each
run. The table export uses `no_ids=true`, `feature:minimal=false`, and
`feature:add_metavalues=0`, preserving RNA identity in feature metadata.

The related `pipeline2.toppas` has the same extraction preset.
`pipeline2_fixed_2.toppas` instead selects a symmetric elution model.
The earlier `forTrent.toppas` uses 10 ppm, automatic RT windows, 5-second peak
width, four isotopes with probability filtering disabled, and a symmetric
model. All use the same per-run identification → extraction → export graph;
none adds alignment or normalization.

The original Herbert example mzML contains MS2 scans only and cannot validate
successful MS1 quantification. Native integration tests instead use the
existing NME1 example, which contains both negative-mode MS1 and MS2 data, and
a compact synthetic RNA isotope-envelope experiment with a known peak area.
The synthetic data are generated during tests under the repository's
BSD-3-Clause license; see [`tests/data/README.md`](../tests/data/README.md).

## Validation

Run the complete pinned-runtime gate with:

```bash
python scripts/check_regressions.py --container
```

The native LFQ tests run real FeatureFinderMetaboIdent and TextExporter tools,
check RNA identity, negative charge, retention time and a known Gaussian area,
exercise empty targets, and execute NASE → extraction → export on NME1.
Ordinary development runs skip native checks when the pinned tools are absent;
the strict regression gate treats every skip as a failure.
