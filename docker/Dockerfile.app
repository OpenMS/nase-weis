# Reuse the full runtime built by OpenMS/streamlit-template's main branch.
# Pin its multi-platform manifest for reproducible local builds. Override
# TEMPLATE_IMAGE to test a newer template runtime after verifying its version.
ARG TEMPLATE_IMAGE=ghcr.io/openms/streamlit-template@sha256:26640815c707bb685d57d197c591956cbc264a4736ab3ccee36c379aa486064e
FROM ${TEMPLATE_IMAGE}
ARG OPENMS_RELEASE=v3.6.0
WORKDIR /app
COPY docker/restore_tool_registry.py /tmp/restore_tool_registry.py
# Work around the upstream image's omitted build-generated registry. Once
# upstream ships it, this branch becomes a no-op. Tool declarations come
# from the same release the pinned runtime builds, not an unrelated checkout.
RUN if ! ls /openms/share/TOOLS/*.tsv >/dev/null 2>&1; then \
      curl --fail --location --retry 3 \
        "https://raw.githubusercontent.com/OpenMS/OpenMS/${OPENMS_RELEASE}/src/topp/executables.cmake" \
        -o /tmp/openms-tools.cmake && \
      /root/miniforge3/envs/streamlit-env/bin/python /tmp/restore_tool_registry.py \
        /tmp/openms-tools.cmake /openms/bin /openms/share/TOOLS/OpenMS.tsv; \
    fi && \
    NucleicAcidSearchEngine -write_ini /tmp/nase-check.ini && \
    DecoyDatabase -write_ini /tmp/decoy-check.ini && \
    rm -f /tmp/restore_tool_registry.py /tmp/openms-tools.cmake /tmp/nase-check.ini /tmp/decoy-check.ini
COPY requirements.txt /app/requirements.txt
RUN /root/miniforge3/envs/streamlit-env/bin/python -m pip install --no-cache-dir -r requirements.txt
COPY assets/ /app/assets/
COPY content/ /app/content/
COPY src/ /app/src/
COPY utils/ /app/utils/
COPY docs/ /app/docs/
COPY example-data/ /app/example-data/
COPY gdpr_consent/ /app/gdpr_consent/
COPY .streamlit/config.toml /app/.streamlit/config.toml
COPY app.py settings.json default-parameters.json presets.json app.yaml /app/
# Local mode keeps named workspaces without CAPTCHA/Redis. Both native and
# Compose launches use the same application path. Hosted deployments can use
# the full source Dockerfile and their own settings overrides.
ENV LOCAL_WORKSPACES_DIR=/workspaces-streamlit-template
EXPOSE 8501
ENTRYPOINT ["/root/miniforge3/envs/streamlit-env/bin/python", "-m", "streamlit", "run", "app.py", "--server.address", "0.0.0.0", "--server.headless", "true"]
