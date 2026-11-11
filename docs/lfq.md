# Optional RNA label-free quantification

On **Configure → Identification**, select the acquisitions to analyze and assign
their conditions and biological replicates. Imported files are available inputs;
the selection determines which are used. **Browse → Your files** also provides
per-file Use checkboxes, Select all and Clear selection.

LFQ is **off by default**. Enable it on **Configure → Quantification** for an
experiment with at least two selected acquisitions containing MS1 scans. The
analysis uses one saved parameter set and a frozen sample design. Native NASE
identification and feature extraction run independently for every selected file,
then their measurements are combined into one modified-RNA-by-sample matrix.
Identification-only analyses may contain one or more acquisitions.

Choose **Total area** normalization or retain **None**. Advanced exposes the
area measure and extraction settings. Merely exporting LFQ targets does not
enable quantification. Changing settings after a run starts does not change its
saved design, parameters, native outputs or comparison.

Each acquisition follows the existing RNA TOPPAS graph:

```text
RNA FASTA → optional DecoyDatabase → NucleicAcidSearchEngine
original mzML ────────────────────→ NucleicAcidSearchEngine
                                   │ lfq_out (target TSV)
original mzML ──────────────────────┼─→ FeatureFinderMetaboIdent
                                   │       │ featureXML
                                   │       └─→ TextExporter → sample feature TSV

sample TSVs → sum areas by modified RNA label → RNA × sample matrix
            → optional total-area normalization → relative condition comparison
```

FeatureFinderMetaboIdent extracts isotope chromatograms from the original
**MS1** data. NASE MS2 deisotoping and MS1 precursor correction settings affect
identification; their exported processed MS2 spectra are not the LFQ input.
An mzML without MS1 scans is rejected when LFQ is enabled. Vendor inputs follow
the native OpenMS readers and must also contain MS1 data.

## Results and interpretation

Each sample retains these files under `samples/sample_NNN/results/`:

| File | Contents |
| --- | --- |
| `exports/lfq-targets.tsv` | Native NASE targets: modified RNA sequence/adduct label, formula, observed charges, and expected retention time |
| `quantification/features.featureXML` | Native detected features, intensities, metadata and quality annotations; viewable in TOPPView |
| `quantification/features.tsv` | TextExporter feature table, including all feature metadata |

The sample TSV includes `label`, `sum_formula`, `charge`, `expected_rt`, measured `mz`,
retention time and `intensity`. Its first column is normally named `#rt`.
The default intensities are unnormalized feature peak areas in instrument
units. They are not absolute concentrations. This table retains the separate
charge/adduct features; the experiment matrix sums their selected areas by
exact modified RNA label after removing the native trailing adduct suffix.
Sequence modifications and terminal chemistry remain part of that identity.
Changing the elution model can change how feature intensity is estimated.

`experiment.json` records sample identities and statuses. Every sample has its
own effective `params.json`, acquisition fingerprint, native outputs and logs.
The experiment's comparison directory contains the raw abundance matrix, the
analysis matrix, sample design and normalization settings/factors. Result
downloads include both matrices; sample drill-down retains the existing
spectra and feature/XIC viewers. Older single-acquisition saved runs remain
readable.

The comparison follows the recovered `NASE_quant` KNIME workflow: sum areas
across charge/adduct states per modified RNA, outer-join samples, and optionally
divide each sample column by its total measured RNA area. An absent or zero-area
measurement remains missing rather than becoming a fabricated abundance zero.
An empty sample's normalization is reported explicitly. The matrix preserves
floating-point areas.

For a condition comparison, technical acquisitions sharing a condition and
biological replicate are averaged first. The app reports biological-replicate
means, standard deviations, observation counts, ratios and log2 fold changes.
An interactive scatter plots denominator means on x and numerator means on y,
with an equality diagonal and optional RNA labels. It uses the saved area
measure and normalization, and offers logarithmic or linear axes. Unpaired or
nonpositive measurements are excluded explicitly, never placed at zero. SVG
and offline HTML exports retain the plotted values. A single-condition
experiment can compare two individual acquisitions instead of condition means.
Ratios require positive measurements in both conditions. Fractionated experiments
are not supported by this RNA aggregation yet; use fraction 1. Duplicate
acquisition paths or ambiguous duplicate basenames are rejected.

The default measurement remains the TOPPAS **native feature intensity** with
the existing extraction preset. The **fitted model area** option requires an
elution model; missing fitted areas cause an error rather than falling back to
another measure. The Advanced button **Use fitted-area RNA comparison preset**
explicitly applies the recovered KNIME settings: symmetric model, 10 ppm m/z
window, automatic RT window, 60-second expected peak width, five isotopes and
isotope-probability filtering disabled. It selects fitted model areas and
total-area normalization. This is a distinct analytical preset, not a silent
replacement of previously saved TOPPAS settings.

A target can have no detected feature, and an identification search can produce
no LFQ targets. An empty target list completes with an empty featureXML and TSV;
it does not synthesize zero-valued abundance observations. Identification
results and targets remain available for inspection. Native errors in
extraction or export cause the run to fail and are recorded in its logs.

If any sample fails, its partial identification outputs remain inspectable and
the experiment is marked failed. A complete-cohort comparison is not published
from a silently reduced sample set.

This workflow does not perform cross-run retention-time alignment, unidentified
feature matching, identification transfer, missing-value imputation or statistical
significance testing. Targets come from each acquisition's own NASE search.
Condition summaries and fold changes describe the observed measurements. The native model option
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
