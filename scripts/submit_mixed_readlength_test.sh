#!/usr/bin/env bash

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"

DATA_ROOT="/data/humangen_sfb1665_seqdata/short_read/processed_data/phase_i"
READ_LENGTH_TABLE="/data/humangen_sfb1665_seqdata/short_read/processed_data/bam_read_lengths.tsv"
SAMPLE_151="A4842_DNA_01"
SAMPLE_150="GS240"
SAMPLE_SHEET="sample_sheets/real_mixed_150_151_samples.tsv"
CONFIG_FILE="configs/config.yaml"
EXPERIMENT="real_mixed_150_151"
MINIFORGE_PATH="${GGTYPER_MINIFORGE_PATH:-/work/hassan/hassan/miniforge}"
BUILD_ENV="${GGTYPER_BUILD_ENV:-genotyping}"
SNAKEMAKE_ENV="${GGTYPER_SNAKEMAKE_ENV:-snakemake}"
BUILD_JOBS="${GGTYPER_BUILD_JOBS:-8}"
BUILD=true
DRY_RUN=true
SUBMIT=true

usage() {
    cat <<'EOF'
Prepare and submit a two-sample GGTyper test with 151-bp and 150-bp reads.

Usage:
  scripts/submit_mixed_readlength_test.sh [options]

Defaults:
  151-bp sample: A4842_DNA_01
  150-bp sample: GS240

Options:
  --sample-151 ID       Sample expected to have 151-bp reads
  --sample-150 ID       Sample expected to have 150-bp reads
  --data-root PATH      Directory containing the phase_i sample directories
  --read-lengths PATH   bam_read_lengths.tsv path
  --sample-sheet PATH   Generated sample sheet, relative to the repository
  --config PATH         Snakemake config, relative to the repository
  --experiment NAME     Output experiment name
  --skip-build          Use the existing ./ggtyper binary
  --skip-dry-run        Submit without a Snakemake dry run
  --dry-run-only        Prepare and dry-run, but do not submit
  -h, --help            Show this help

Environment overrides:
  GGTYPER_MINIFORGE_PATH  Miniforge installation (default /work/hassan/hassan/miniforge)
  GGTYPER_BUILD_ENV       Conda build environment (default genotyping)
  GGTYPER_SNAKEMAKE_ENV   Conda Snakemake environment (default snakemake)
  GGTYPER_BUILD_JOBS      Parallel make jobs (default 8)

Typical cluster command:
  bash scripts/submit_mixed_readlength_test.sh
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

absolute_from_repo() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s/%s\n' "$REPO_ROOT" "$1" ;;
    esac
}

validate_sample_id() {
    local sample_id="$1"
    [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] ||
        die "Unsafe sample ID: $sample_id"
}

show_and_check_read_length() {
    local sample_id="$1"
    local expected="$2"
    local matching_rows

    matching_rows="$(awk -F '\t' -v sample="$sample_id" 'index($0, sample) { print }' "$READ_LENGTH_TABLE")"
    [[ -n "$matching_rows" ]] ||
        die "No read-length row found for $sample_id in $READ_LENGTH_TABLE"

    printf '\nRead-length inventory row(s) for %s:\n%s\n' "$sample_id" "$matching_rows"

    if ! awk -F '\t' -v sample="$sample_id" -v expected="$expected" '
        index($0, sample) {
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^[0-9]+([.]0+)?$/ && ($i + 0) == expected) found = 1
            }
        }
        END { exit(found ? 0 : 1) }
    ' "$READ_LENGTH_TABLE"; then
        die "$sample_id does not have a standalone $expected-bp value in $READ_LENGTH_TABLE"
    fi
}

activate_conda_env() {
    local environment="$1"
    local conda_setup="$MINIFORGE_PATH/etc/profile.d/conda.sh"

    [[ -r "$conda_setup" ]] || die "Cannot read Conda setup script: $conda_setup"
    # shellcheck disable=SC1090
    source "$conda_setup"
    conda activate "$environment"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-151)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SAMPLE_151="$2"
            shift 2
            ;;
        --sample-150)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SAMPLE_150="$2"
            shift 2
            ;;
        --data-root)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            DATA_ROOT="$2"
            shift 2
            ;;
        --read-lengths)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            READ_LENGTH_TABLE="$2"
            shift 2
            ;;
        --sample-sheet)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        --config)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            CONFIG_FILE="$2"
            shift 2
            ;;
        --experiment)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            EXPERIMENT="$2"
            shift 2
            ;;
        --skip-build)
            BUILD=false
            shift
            ;;
        --skip-dry-run)
            DRY_RUN=false
            shift
            ;;
        --dry-run-only)
            DRY_RUN=true
            SUBMIT=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1 (use --help)"
            ;;
    esac
done

validate_sample_id "$SAMPLE_151"
validate_sample_id "$SAMPLE_150"
[[ "$SAMPLE_151" != "$SAMPLE_150" ]] || die "The 151-bp and 150-bp samples must differ"
[[ "$EXPERIMENT" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe experiment name: $EXPERIMENT"

SAMPLE_SHEET="$(absolute_from_repo "$SAMPLE_SHEET")"
CONFIG_FILE="$(absolute_from_repo "$CONFIG_FILE")"

BAM_151="$DATA_ROOT/$SAMPLE_151/alignment/$SAMPLE_151.markdup.bam"
VCF_151="$DATA_ROOT/$SAMPLE_151/sv_calls/$SAMPLE_151.manta.vcf.gz"
BAM_150="$DATA_ROOT/$SAMPLE_150/alignment/$SAMPLE_150.markdup.bam"
VCF_150="$DATA_ROOT/$SAMPLE_150/sv_calls/$SAMPLE_150.manta.vcf.gz"

cd "$REPO_ROOT"

[[ -r "$READ_LENGTH_TABLE" ]] || die "Cannot read $READ_LENGTH_TABLE"
[[ -r "$CONFIG_FILE" ]] || die "Cannot read config file $CONFIG_FILE"
[[ -r "$REPO_ROOT/submit_job.sh" ]] || die "Cannot read submit_job.sh"

for required_file in \
    "$BAM_151" "$BAM_151.bai" "$VCF_151" "$VCF_151.tbi" \
    "$BAM_150" "$BAM_150.bai" "$VCF_150" "$VCF_150.tbi"; do
    [[ -r "$required_file" ]] || die "Required input is missing or unreadable: $required_file"
done

show_and_check_read_length "$SAMPLE_151" 151
show_and_check_read_length "$SAMPLE_150" 150

if ! grep -Eq 'tolerance[[:space:]]*=[[:space:]]*1' src/readLengthPolicy.hpp; then
    die "This checkout does not contain the 1-bp read-length tolerance"
fi

if [[ "$BUILD" == true ]]; then
    printf '\nBuilding GGTyper in Conda environment %s...\n' "$BUILD_ENV"
    activate_conda_env "$BUILD_ENV"
    make -B -j "$BUILD_JOBS" \
        INCLUDE_PATH="$CONDA_PREFIX/include" \
        LIB_PATH="$CONDA_PREFIX/lib"
fi

[[ -x "$REPO_ROOT/ggtyper" ]] ||
    die "The GGTyper executable is missing. Run without --skip-build."

mkdir -p "$(dirname -- "$SAMPLE_SHEET")" logs
printf 'sample_id\tbam_path\tvcf_path\n%s\t%s\t%s\n%s\t%s\t%s\n' \
    "$SAMPLE_151" "$BAM_151" "$VCF_151" \
    "$SAMPLE_150" "$BAM_150" "$VCF_150" > "$SAMPLE_SHEET"

printf '\nGenerated sample sheet:\n'
column -t -s $'\t' "$SAMPLE_SHEET" 2>/dev/null || sed -n '1,3p' "$SAMPLE_SHEET"

if [[ "$DRY_RUN" == true ]]; then
    printf '\nRunning Snakemake dry run in Conda environment %s...\n' "$SNAKEMAKE_ENV"
    activate_conda_env "$SNAKEMAKE_ENV"
    snakemake \
        --snakefile workflow/Snakefile \
        --configfile "$CONFIG_FILE" \
        --config samples_sheet="$SAMPLE_SHEET" experiment_name="$EXPERIMENT" \
        --cores 1 \
        --dry-run \
        --printshellcmds
fi

if [[ "$SUBMIT" == true ]]; then
    command -v sbatch >/dev/null 2>&1 || die "sbatch is not available on this host"
    printf '\nSubmitting Slurm job...\n'
    sbatch \
        --job-name=ggtyper_mix_test \
        submit_job.sh \
        "$CONFIG_FILE" \
        "$SAMPLE_SHEET" \
        "$EXPERIMENT"
    printf 'Results will be written to: %s/output/%s/\n' "$REPO_ROOT" "$EXPERIMENT"
else
    printf '\nDry run completed; no Slurm job was submitted.\n'
fi
