# Acquisition and result viewers

NASEWEIS reuses the template's installed **pyopenms-viz 1.0.0** and Plotly
backend. No additional visualization dependency or external visualization
service is required. Plots support zoom, pan, hover and image export; their
HTML downloads include the plotting code and work offline.

## Raw data

Open **Raw data viewer**, choose an uploaded or locally registered mzML, and
select **Load acquisition**. Viewing does not change the input selected for
the next search.

- **MS1 peak map:** retention time versus m/z, colored by summed signal in
  display bins. Optional marginals show the corresponding aggregated signal.
- **Spectra:** filter by MS level or native ID and inspect an individual scan.
  Empty scans remain visible in the index. The displayed spectrum is unbinned;
  relative-intensity display normalization does not modify the native peaks.
- **Chromatograms:** display TIC and BPC together or individually and add a
  full-resolution XIC with explicit ±ppm or ±Da tolerance. Every MS1 scan
  contributes, including zero signal where no peaks fall in the mass window.

The peak-map overview is limited to 40,000 occupied RT/m/z cells. Its bin-center
coordinates are an overview, not individual peak coordinates. Zooming this
overview does not reconstruct peaks within a bin; use the spectrum or XIC
viewer for full-resolution values. TIC/BPC and XIC calculations use all
original MS1 peaks, independently of overview binning.

Indexed mzML files provide metadata and individual spectrum access without
loading the experiment into memory. Unindexed files use native streaming
consumers. Only one scan and the bounded overview are retained while reading;
the app does not construct whole-experiment long MS1/MS2 peak tables.

These viewers accept mzML. Thermo RAW and Bruker inputs remain supported by the
NASE search workflow; convert those acquisitions to mzML for visualization.
MS2-only acquisitions can still be browsed in the spectrum viewer, with an
explanation where an MS1 view is unavailable.

## Identification evidence

Open **Results → Sequence coverage** for an experiment overview, with or without
LFQ. The wrapped reference sequence has packed oligonucleotide bars grouped by
condition or acquisition, a shared spectrum-count color scale and markers for
native RNA modifications, including bracketed bases and D/I/Y codes. Hover shows exact sequences, coordinates
and counts. Terminal phosphate groups are not modification sites. Counts pool
identification evidence across acquisitions; they are not normalized abundances.
Ambiguous mappings retain shared evidence at each native position. Conflicting
reference sequences cannot share a coverage axis. Coverage is also available
for each individual sample and older single-acquisition runs. Download the
self-contained HTML or vector SVG; neither requires the original acquisition.

Open **Results → Spectra** for the selected saved run, or
**Results → Sample details → Spectra** in a multi-sample experiment. Select
**Open spectrum viewer** to load the RNA identification viewer. Search by RNA
sequence or complete native spectrum ID and choose a spectrum match; its
selector includes the signed precursor charge, retention time and native ID.
The acquisition path
comes only from that run's parameter snapshot. Missing source data or an
unreadable snapshot produces a clear message; the current search input is
never silently substituted.

New mzML searches record a SHA256 fingerprint of the acquisition in their
parameter snapshot. Saved-run viewers verify those bytes before showing scan
or XIC evidence, so replacing a same-named upload does not silently change an
old run's source. An identical restored file is accepted. Older runs without
a recorded fingerprint are explicitly marked as unverified; confirm that the
available acquisition is the original one before interpreting their evidence.

The mzTab `ms_run[n]:` prefix is removed from the spectrum reference, and the
complete remaining native ID must identify exactly one MS2 scan. Duplicate or
missing IDs are reported. NASE's internal MS2 array index is not a raw scan
index, and retention time alone does not establish scan identity.

Native RNA annotations come from the selected hit's idXML metadata `label`,
signed precursor charge and peak assignments. The peptide sequence field is
empty for RNA. The viewer preserves multiple assignments to one observed peak
and supports RNA `a-B`, `d` and `w` ion series. It labels only unique m/z matches
within an absolute 0.00001 Da serialization allowance and keeps the acquired
peak heights. NASE's stored annotation intensities describe its scored
spectrum and are not used to replace raw intensities.

The visible **Assigned fragment ions** table lists ion labels, signed fragment
charges, native assigned m/z, matched raw m/z, raw intensity and relative
intensity against the full scan's base peak. Use **Ion series** to filter the
table and spectrum labels together. Select one or more table rows with the
checkboxes to highlight their exact acquired peaks in gold; multiple ions can
refer to one peak. Unplaced or ambiguous assignments remain listed with their
status and cannot highlight another peak. The original native metadata and TSV
download retain all assignments, independently of display filters. Display
settings do not change identification results or saved scientific parameters.

Deisotoping or decharging can produce scored masses absent from the raw scan.
Those assignments remain in the downloadable native fragment table and are
not attached to unrelated raw peaks. Conflicting duplicate native hits are
reported rather than chosen arbitrarily. The viewer does not reconstruct
theoretical spectra or pair experimental/theoretical debug exports by position;
those files can contain pre-filter candidates and parallel output ordering.

## LFQ evidence

Multi-sample **Results → Quantification** includes a condition-versus-condition
RNA abundance scatter, with the denominator on x, numerator on y and an equality
diagonal. It uses the same biological-replicate means, area measure and
normalization as the saved comparison table. Both logarithmic and linear views
omit targets without positive measurements in both groups, with the omitted
count stated explicitly. Open markers are unmodified RNAs; filled markers have
native modification or ambiguity annotations. Select identities to
label, hover for means and replicate counts, and export SVG from the camera
button or download the offline HTML. Experiments with one condition can instead
compare two individual acquisitions; those values are not condition means.

Open **Results → Quantification → Open LFQ feature viewer**. Each feature
retains its native charge, label, RT, m/z and raw area; coincident features are
not merged. Choose a feature to extract its original MS1 signal. The default
XIC half-width is half the saved full extraction window, respecting OpenMS's
ppm/Da convention. The detected feature boundaries and apex are marked.

The XIC contains signal around the selected m/z. Native LFQ areas can include
several isotope traces and model fitting, so they need not equal the integral
of that single-m/z XIC. Adjusting viewer tolerances or the displayed RT range
does not modify the feature map or rerun quantification. A feature map remains
viewable even if its original acquisition is no longer available.

## Isolation and validation

Viewer state is scoped to the workspace, file and saved run. Cached tables
include the resolved source path, size, modification/change times and file
identity. Saved-run reads also bind to the verified digest; platforms without
reliable change timestamps recompute it when evidence is opened. Timestamp-
preserving replacements invalidate reads, and returned tables are copies. Online
viewers reject acquisitions outside the current workspace, including symlink
escapes. Local registered external inputs remain available.

Tests exercise actual pyopenms-viz figures, indexed and unindexed native mzML,
exact native-ID/charge/RNA hit matching, empty and malformed data, multiple
traces, file/run switches and safe offline HTML exports. A pinned-runtime
NME1 integration test compares spectra, TIC/BPC and XIC values against an
independently loaded original experiment and checks native RNA assignments.
Run the complete gate with `python scripts/check_regressions.py --container`.
