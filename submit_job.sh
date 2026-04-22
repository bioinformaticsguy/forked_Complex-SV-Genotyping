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
#   sbatch submit_job.sh [config.yaml] [samples.tsv]
#
# Examples:
#   sbatch submit_job.sh                         # uses defaults
#   sbatch submit_job.sh config.yaml samples.tsv
#
# Override resources at submission time:
#   sbatch --mem=128GB --time=5-00:00:00 submit_job.sh
# ============================================================

CONFIGFILE="${1:-config.yaml}"
SAMPLESFILE="${2:-samples.tsv}"

# --- Conda setup ---
MINIFORGE_PATH="/root/miniforge3"
source "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
conda activate genotyping

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
echo "Run log:           $RUN_LOG"
echo "=============================================="

# --- Unlock in case a previous job was killed or timed out ---
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIGFILE" \
    --config samples_sheet="$SAMPLESFILE" \
    --unlock 2>/dev/null || true

# --- Run ---
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIGFILE" \
    --config samples_sheet="$SAMPLESFILE" \
    --use-conda \
    --conda-prefix /root/miniforge3/envs \
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
