#!/bin/bash

# --- Setup ---
CONFIG_FILE="${PIPELINE_CONFIG:-pipeline.conf}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

OUT_DIR="${OUT_DIR:-dump/admixture}"
mkdir -p "$OUT_DIR"
PREFIX="${PREFIX:-${ADMIXTURE_PREFIX:-admixture}}"
KLIST="${KLIST:-2 3 4 5 6 7 8 9 10}"
THREADS="${THREADS:-1}"
RANDOM_SEED="${RANDOM_SEED:-43}"
ADMIXTURE_BIN="${ADMIXTURE_BIN:-$(command -v admixture || true)}"
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
DATA_PREPROCESS_DIR="${DATA_PREPROCESS_DIR:-dump/common}"
DATA_PREPROCESS_PREFIX="${DATA_PREPROCESS_PREFIX:-common}"
INPUT_BED="${INPUT_BED:-$DATA_PREPROCESS_DIR/${DATA_PREPROCESS_PREFIX}_ALL.pruned.bed}"

if [[ -z "$ADMIXTURE_BIN" && -x "$HOME/envs/bio_tools/bin/admixture" ]]; then
    ADMIXTURE_BIN="$HOME/envs/bio_tools/bin/admixture"
fi

if [[ -z "$ADMIXTURE_BIN" || ! -x "$ADMIXTURE_BIN" ]]; then
    echo "ERROR: admixture not found or not executable."
    echo "Please activate your conda env (e.g., conda activate bio_tools) or set ADMIXTURE_BIN in pipeline.conf."
    exit 1
fi

if [[ ! -f "$INPUT_BED" ]]; then
    echo "ERROR: shared preprocessed input not found: $INPUT_BED"
    echo "Please run preprocess.sh first, for example from analysis.ipynb."
    exit 1
fi
INPUT_BED_ABS="$(cd "$(dirname "$INPUT_BED")" && pwd)/$(basename "$INPUT_BED")"

# shellcheck disable=SC1091
source "$(cd "$(dirname "$0")" && pwd)/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi

# Convert to absolute path so it stays valid after pushd.
METRICS_FILE="$(cd "$(dirname "$METRICS_FILE")" && pwd)/$(basename "$METRICS_FILE")"
pushd "$OUT_DIR" >/dev/null

run_one_k() {
    local K="$1"
    local fit_exit_code=0

    echo "Running ADMIXTURE for K=$K"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "admixture" "fit" "$K" "cv_log_K${K}.out" \
            "$ADMIXTURE_BIN" --cv "-j1" "--seed=${RANDOM_SEED}" "$INPUT_BED_ABS" "$K"
        fit_exit_code=$?
    else
        "$ADMIXTURE_BIN" --cv "-j1" "--seed=${RANDOM_SEED}" "$INPUT_BED_ABS" "$K" | tee "cv_log_K${K}.out"
        fit_exit_code=${PIPESTATUS[0]}
    fi

    if [[ "$fit_exit_code" -ne 0 ]]; then
        echo "ERROR: ADMIXTURE failed at K=$K with exit code $fit_exit_code. See $OUT_DIR/cv_log_K${K}.out"
    fi
    return "$fit_exit_code"
}

# Run one single-threaded job per K in parallel.
declare -a job_pids=()
declare -a job_ks=()

for K in $KLIST; do
    run_one_k "$K" &
    job_pids+=("$!")
    job_ks+=("$K")
done

overall_exit_code=0
for i in "${!job_pids[@]}"; do
    if ! wait "${job_pids[$i]}"; then
        overall_exit_code=1
        echo "ERROR: ADMIXTURE background job failed for K=${job_ks[$i]}"
    fi
done

if [[ "$overall_exit_code" -ne 0 ]]; then
    exit "$overall_exit_code"
fi

# Extract CV results for easy viewing
grep "CV error" cv_log_K*.out | sed 's/cv_log_K\(.*\).out:CV error (K=\(.*\)): \(.*\)/\2 \3/' | sort -n > cv_results.txt

popd >/dev/null

echo "Pipeline complete. Results and logs are in $OUT_DIR/"