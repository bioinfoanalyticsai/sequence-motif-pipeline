# Multi-platform support
FROM --platform=linux/amd64 python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git wget build-essential dos2unix \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    pandas \
    biopython \
    logomaker \
    matplotlib \
    numpy \
    pyfaidx

# Working directory inside container
WORKDIR /data

# Copy scripts into container
COPY scripts/ /usr/local/bin/scripts/

# Fix line endings and permissions
RUN dos2unix /usr/local/bin/scripts/*.py && chmod +x /usr/local/bin/scripts/*.py

# No default ENTRYPOINT
