#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem=64GB
#SBATCH --job-name=ggtyper
#SBATCH --output=logs/slurm_%j_%u_%N.out
#SBATCH --error=logs/slurm_%j_%u_%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com

# ============================================================
# GGTyper multi-sample SV genotyping — SLURM submission script
# ============================================================
# Usage:
#   sbatch submit_job.sh [config.yaml] [sample_sheets/samples.tsv]
#
# Examples:
#   sbatch submit_job.sh                                          # uses defaults
#   sbatch submit_job.sh config.yaml sample_sheets/samples.tsv
#   sbatch submit_job.sh config.yaml sample_sheets/cohort1_samples.tsv
#
# Override experiment name (default: derived from sheet filename):
#   sbatch submit_job.sh config.yaml sample_sheets/cohort1_samples.tsv cohort1_v2
#
# Override resources at submission time:
#   sbatch --mem=128GB --time=5-00:00:00 submit_job.sh
# ============================================================

CONFIGFILE="${1:-config.yaml}"
SAMPLESFILE="${2:-sample_sheets/samples.tsv}"
EXPERIMENT="${3:-}"   # optional; overrides the name derived from the sheet filename

# --- Conda setup ---
MINIFORGE_PATH="/work/hassan/hassan/miniforge"
source "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
conda activate snakemake

# --- Setup ---
mkdir -p logs

TIMESTAMP=$(date +%Y%m%dT%H%M%S)
RUN_LOG="logs/run_${TIMESTAMP}_${SLURM_JOB_ID:-interactive}.log"
exec > >(tee -a "$RUN_LOG") 2>&1

echo "=============================================="
echo "  GGTyper — Multi-sample SV Genotyping"
echo "=============================================="
echo "Job ID:            ${SLURM_JOB_ID:-interactive}"
echo "Node:              $(hostname)"
echo "Working directory: $(pwd)"
echo "Snakemake version: $(snakemake --version)"
echo "Date:              $(date)"
echo "CPUs:              ${SLURM_CPUS_PER_TASK}"
echo "Memory:            ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Config:            $CONFIGFILE"
echo "Samples:           $SAMPLESFILE"
echo "Experiment:        ${EXPERIMENT:-<derived from sheet name>}"
echo "Run log:           $RUN_LOG"
echo "=============================================="

# Build optional experiment_name override
EXTRA_CONFIG="samples_sheet=\"$SAMPLESFILE\""
[ -n "$EXPERIMENT" ] && EXTRA_CONFIG="$EXTRA_CONFIG experiment_name=\"$EXPERIMENT\""

# --- Unlock in case a previous job was killed or timed out ---
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIGFILE" \
    --config $EXTRA_CONFIG \
    --unlock 2>/dev/null || true

# --- Run ---
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIGFILE" \
    --config $EXTRA_CONFIG \
    --use-conda \
    --conda-prefix /work/hassan/hassan/snakemake-conda \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds

EXIT_CODE=$?

echo "=============================================="
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "Pipeline finished successfully"
else
    echo "Pipeline FAILED with exit code $EXIT_CODE"
    echo ""
    echo "--- Last 50 lines of Snakemake log ---"
    SMLOG=$(ls -t .snakemake/log/*.snakemake.log 2>/dev/null | head -1)
    if [ -n "$SMLOG" ]; then
        echo "Log: $SMLOG"
        tail -50 "$SMLOG"
    fi
fi
echo "=============================================="
echo "Full run log: $RUN_LOG"
echo "=============================================="

exit "$EXIT_CODE"
