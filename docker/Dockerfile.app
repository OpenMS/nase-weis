# Python/UI dependencies and nightly TOPP tools have independent pinned roots.
# The Python ParamXML API is kept at pyOpenMS 3.6; NASE/DecoyDatabase and their
# native libraries/data all come from the same verified nightly installation.
ARG TEMPLATE_IMAGE=ghcr.io/openms/streamlit-template@sha256:26640815c707bb685d57d197c591956cbc264a4736ab3ccee36c379aa486064e
ARG OPENMS_NIGHTLY_IMAGE=ghcr.io/openms/openms-executables@sha256:2328f9e3ffe968b92455009d7a327f02fe973d150f8bb5784c3fe37c6288f97d
FROM ${TEMPLATE_IMAGE} AS python-runtime
FROM ${OPENMS_NIGHTLY_IMAGE}
USER root
COPY --from=python-runtime /root/miniforge3/envs/streamlit-env /root/miniforge3/envs/streamlit-env
ENV PATH="/opt/OpenMS/bin:/root/miniforge3/envs/streamlit-env/bin:${PATH}"
WORKDIR /app
RUN NucleicAcidSearchEngine -write_ini /tmp/nase-check.ini && \
    DecoyDatabase -write_ini /tmp/decoy-check.ini && \
    rm -f /tmp/nase-check.ini /tmp/decoy-check.ini
COPY requirements.txt /app/requirements.txt
RUN /root/miniforge3/envs/streamlit-env/bin/python -m pip install --no-cache-dir -r requirements.txt && \
    /root/miniforge3/envs/streamlit-env/bin/python -m pip install --no-cache-dir --force-reinstall --no-deps pyopenms==3.6.0 && \
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
