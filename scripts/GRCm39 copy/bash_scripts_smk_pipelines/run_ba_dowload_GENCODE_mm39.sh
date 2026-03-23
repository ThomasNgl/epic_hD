#!/bin/bash
set -e  # Exit if any command fails

# Always run from the directory where this script is located
cd "$(dirname "$0")"

# Get script name, remove 'run_' prefix and '.sh' extension
SCRIPT_BASENAME=$(basename "$0")
SNAKEMAKE_NAME=${SCRIPT_BASENAME#run_}           # remove 'run_' prefix
SNAKEMAKE_NAME=${SNAKEMAKE_NAME%.sh}             # remove '.sh' extension

# Use the derived name for the Snakefile
snakemake --cores 6 \
    -s "${SNAKEMAKE_NAME}.smk" \
    --use-conda \
    --conda-prefix ~/snakemake-conda-envs \
    --rerun-incomplete

echo "Snakemake workflow completed!"
