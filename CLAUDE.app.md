# NASEWEIS

RNA oligonucleotide searches using DecoyDatabase and NucleicAcidSearchEngine.
App-specific workflow: `src/Workflow.py`; mzTab results and coverage: `src/naseweis_results.py`.
Run tests with `.venv/bin/python -m pytest test_gui.py tests/`.
Real TOPP integration: `.venv/bin/python -m pytest tests/test_naseweis_integration.py -v`.
The required TOPP tools must be on PATH. Start locally with `streamlit run app.py`.
