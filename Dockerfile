FROM python:3.12-slim@sha256:78387bc3881b8273120a12ebe6c1ab22b018ccc2c9adf565ae1ac9b536e184ea

LABEL org.opencontainers.image.title="ORFBounder" \
      org.opencontainers.image.description="Bacterial TIS and TTS peak calling" \
      org.opencontainers.image.version="2.0.0"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /opt/orfbounder
COPY pyproject.toml uv.lock README.md LICENSE ./
COPY orfbounder.py orfbounder_batch.py ./
COPY lib ./lib
COPY helpers ./helpers
RUN python -m pip install --no-cache-dir uv==0.9.7 \
    && uv sync --frozen --no-dev --no-editable --no-cache

ENV PATH="/opt/orfbounder/.venv/bin:${PATH}"

WORKDIR /work
# A command instead of ENTRYPOINT also supports Snakemake/Apptainer shell jobs.
CMD ["orfbounder", "--help"]
