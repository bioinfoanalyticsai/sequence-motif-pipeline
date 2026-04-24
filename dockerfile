# ── Multi-arch base (arm64 + amd64) ──────────────────────────────────────────
FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="your-lab@example.com"
LABEL description="SeqMatcher pipeline: sequence fetch, motif search, PWM scan, BED/BigWig output"
LABEL version="2.2.0"

USER root

# ── System packages ───────────────────────────────────────────────────────────
# procps  → `ps` binary required by Nextflow task-metrics collection
# ca-certificates → needed for conda HTTPS channel downloads
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ── Conda environment (all bioinformatics tools incl. bigWigInfo) ─────────────
# micromamba resolves the correct arch (arm64 / amd64) automatically.
# ucsc-bigwiginfo is available for both architectures on bioconda.
COPY envs/seqmatcher.yml /tmp/seqmatcher.yml
RUN micromamba install -y -n base -f /tmp/seqmatcher.yml \
    && micromamba clean --all --yes

# Make conda-installed binaries the default on PATH
ENV PATH="/opt/conda/bin:${PATH}"

# ── Pipeline scripts ──────────────────────────────────────────────────────────
COPY scripts/ /usr/local/bin/scripts/
RUN chmod +x /usr/local/bin/scripts/*.py

WORKDIR /data
