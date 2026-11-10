"""Introduce the RNA search workflow and its bundled example."""

import streamlit as st

from src.common.common import page_setup
from src.Workflow import Workflow

page_setup(page="main")
st.caption("OPENMS WEBAPPS · OLIGONUCLEOTIDE IDENTIFICATION")
st.title("NASEWEIS")
st.markdown(
    "Identify RNA oligonucleotides in mass spectrometry data with the **OpenMS NucleicAcidSearchEngine**."
)

with st.container(border=True):
    st.subheader("From spectra to sequence coverage")
    cols = st.columns(3)
    for col, title, text in zip(
        cols,
        ["1 · Add your data", "2 · Set up a search", "3 · Explore the evidence"],
        [
            "Upload an MS/MS run and RNA FASTA or cached digest, or try the NME1 example.",
            "Choose digestion, modifications and mass tolerances. Add decoys for FDR estimation.",
            "Review matches and coverage. Export identification databases, spectra and modification tracks.",
        ],
    ):
        col.markdown(f"**{title}**")
        col.write(text)
    st.page_link(
        "content/topp_workflow_file_upload.py",
        label="Start with your data",
        icon=":material/upload_file:",
    )

st.subheader("Try the NME1 example")
st.write(
    "A small RNA dataset and matching FASTA are included. Load both, review the settings and start a search."
)
if st.button("Load example dataset", type="primary", icon=":material/science:"):
    Workflow().load_example_data()
    st.switch_page("content/naseweis.py")

with st.expander("Input formats and search guidance"):
    st.markdown("""
- **Spectra:** one `.mzML` or Thermo `.raw` file, or a Bruker `.d` directory. Upload `.d` directories as ZIP archives, or import them from a local folder. Vendor readers depend on the OpenMS build and platform.
- **Sequences:** one RNA FASTA database (`.fasta`, `.fa`, `.faa`), or a reusable NASE `.oms` digest. Use A, U, C and G; put MODOMICS short names for modified residues in brackets, such as `AU[Am]C[m5C]G`.
- **Digestion:** select the enzyme used in your experiment, or **no cleavage** for intact oligonucleotides.
- **Mass accuracy:** enter precursor and fragment tolerances and their units directly. The defaults are 10 ppm; adjust them for your instrument and acquisition method.
- **Advanced search:** use **All parameters** for FLASH precursor correction/deisotoping, MVH scoring, fragment coverage, peak filters, modification ambiguity resolution and bedRMod reference metadata. Choose optional output files on Configure.
- **Database preparation:** the **Decoy Database** page runs the complete tool independently, including multi-FASTA input, decoy ratios, protein entrapment and neighbor searches.
- **FDR:** generation and estimation are independent: turn generation off for a database that already includes decoys. The default target-decoy search retains matches at q ≤ 0.05. A small example database illustrates the workflow; it does not establish a reliable FDR estimate for an experiment.

Parameters, input files and completed runs are kept in your workspace. Each new run keeps its own results and parameter snapshot.
""")
st.caption(
    "Powered by OpenMS and pyOpenMS · [NucleicAcidSearchEngine publication](https://doi.org/10.1038/s41467-020-14665-7) · [Source code](https://github.com/poshul/nase-weis)"
)
