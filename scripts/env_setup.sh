#!/usr/bin/env bash
# Create a clean conda environment 4 this project & install all required command-line tools/Python libraries.

set -euo pipefail

ENV_NAME="filterbt"
PYVER="3.10"

# Ensure conda is available in this shell
if ! command -v conda >/dev/null 2>&1; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
  else
    echo "conda is not available in the current shell. Plz install/initialize conda before running this script."
    exit 1fi
  fi
fi

# 1. Create env
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    conda create -y -n "$ENV_NAME" "python=${PYVER}"
fi

# 2. Activate env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# 3. Install
conda install -y -c conda-forge -c bioconda \
    ncbi-datasets-cli \
    blast \
    parallel \
    pigz \
    biopython \
    pandas \
    requests \
    pyyaml \
    pip

# 4. Alternate downloader
python -m pip install -U pip
python -m pip install ncbi-genome-download

echo "The environment is ready: conda activate $ENV_NAME"