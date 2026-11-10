# Python/UI dependencies and nightly TOPP tools have independent pinned roots.
# Python bindings and native tools use the verified 2026-11-10 nightly.
# Native libraries/data stay with their matching executables.
ARG TEMPLATE_IMAGE=ghcr.io/openms/streamlit-template@sha256:26640815c707bb685d57d197c591956cbc264a4736ab3ccee36c379aa486064e
ARG OPENMS_NIGHTLY_IMAGE=ghcr.io/openms/openms-executables@sha256:8e3c559cc474a57e12eff10569e2aef033eb867783ab36d45df08e8279880c3d
FROM ${TEMPLATE_IMAGE} AS python-runtime
FROM ${OPENMS_NIGHTLY_IMAGE}
USER root
COPY --from=python-runtime /root/miniforge3/envs/streamlit-env /root/miniforge3/envs/streamlit-env
ENV PATH="/opt/OpenMS/bin:/root/miniforge3/envs/streamlit-env/bin:${PATH}"
WORKDIR /app
RUN NucleicAcidSearchEngine -write_ini /tmp/nase-check.ini && \
    DecoyDatabase -write_ini /tmp/decoy-check.ini && \
    FeatureFinderMetaboIdent -write_ini /tmp/lfq-check.ini && \
    TextExporter -write_ini /tmp/export-check.ini && \
    rm -f /tmp/nase-check.ini /tmp/decoy-check.ini /tmp/lfq-check.ini /tmp/export-check.ini
COPY requirements.txt /app/requirements.txt
RUN /root/miniforge3/envs/streamlit-env/bin/python -m pip install --no-cache-dir -r requirements.txt && \
    /root/miniforge3/envs/streamlit-env/bin/python -c "import pyopenms; print(pyopenms.__version__)"
COPY assets/ /app/assets/
COPY content/ /app/content/
COPY src/ /app/src/
COPY utils/ /app/utils/
COPY docs/ /app/docs/
COPY example-data/ /app/example-data/
COPY gdpr_consent/ /app/gdpr_consent/
COPY .streamlit/config.toml /app/.streamlit/config.toml
COPY app.py settings.json default-parameters.json presets.json app.yaml openms-runtime.json /app/
ENV LOCAL_WORKSPACES_DIR=/workspaces-streamlit-template
EXPOSE 8501
ENTRYPOINT ["/root/miniforge3/envs/streamlit-env/bin/python", "-m", "streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.headless", "true"]
