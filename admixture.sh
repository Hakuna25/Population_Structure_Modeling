#!/bin/bash

# --- Setup ---
CONFIG_FILE="${PIPELINE_CONFIG:-pipeline.conf}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

OUT_DIR="dump/admixture"
mkdir -p "$OUT_DIR"
sample_info="${SAMPLE_INFO:-1000Genomes/igsr_samples.tsv}"
PREFIX="${PREFIX:-${ADMIXTURE_PREFIX:-admixture}}"
KLIST="${KLIST:-2 3 4 5 6 7 8 9 10}"
THREADS="${THREADS:-10}"
RANDOM_SEED="${RANDOM_SEED:-43}"
ADMIXTURE_BIN="${ADMIXTURE_BIN:-$(command -v admixture || true)}"
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [[ -z "$ADMIXTURE_BIN" && -x "$HOME/envs/bio_tools/bin/admixture" ]]; then
    ADMIXTURE_BIN="$HOME/envs/bio_tools/bin/admixture"
fi

if [[ -z "$ADMIXTURE_BIN" || ! -x "$ADMIXTURE_BIN" ]]; then
    echo "ERROR: admixture not found or not executable."
    echo "Please activate your conda env (e.g., conda activate bio_tools) or set ADMIXTURE_BIN in pipeline.conf."
    exit 1
fi

# shellcheck disable=SC1091
source "$SCRIPT_DIR/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi

#  Preprocess
if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    benchmark_run "$METRICS_FILE" "admixture" "preprocess" "NA" "$OUT_DIR/preprocess.log" \
        bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
    preprocess_exit_code=$?
else
    bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
    preprocess_exit_code=$?
fi

if [[ "$preprocess_exit_code" -ne 0 ]]; then
    echo "ERROR: preprocess failed with exit code $preprocess_exit_code. See $OUT_DIR/preprocess.log"
    exit "$preprocess_exit_code"
fi

# Run ADMIXTURE
# Convert to absolute path so it stays valid after pushd
METRICS_FILE="$(cd "$(dirname "$METRICS_FILE")" && pwd)/$(basename "$METRICS_FILE")"
pushd "$OUT_DIR" >/dev/null

for K in $KLIST; do
    echo "Running ADMIXTURE for K=$K"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "admixture" "fit" "$K" "cv_log_K${K}.out" \
            "$ADMIXTURE_BIN" --cv "-j${THREADS}" "--seed=${RANDOM_SEED}" "${PREFIX}_ALL.pruned.bed" "$K"
        fit_exit_code=$?
    else
        "$ADMIXTURE_BIN" --cv "-j${THREADS}" "--seed=${RANDOM_SEED}" "${PREFIX}_ALL.pruned.bed" "$K" | tee "cv_log_K${K}.out"
        fit_exit_code=${PIPESTATUS[0]}
    fi

    if [[ "$fit_exit_code" -ne 0 ]]; then
        echo "ERROR: ADMIXTURE failed at K=$K with exit code $fit_exit_code. See $OUT_DIR/cv_log_K${K}.out"
        exit "$fit_exit_code"
    fi
done

# Extract CV results for easy viewing
grep "CV error" cv_log_K*.out | sed 's/cv_log_K\(.*\).out:CV error (K=\(.*\)): \(.*\)/\2 \3/' | sort -n > cv_results.txt
popd >/dev/null

echo "Pipeline complete. Results and logs are in dump/admixture/"