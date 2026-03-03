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
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# shellcheck disable=SC1091
source "$SCRIPT_DIR/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi

#  Preprocess
if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    benchmark_run "$METRICS_FILE" "admixture" "preprocess" "NA" "$OUT_DIR/preprocess.log" \
        bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
else
    bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
fi

# 4. Run ADMIXTURE
pushd "$OUT_DIR" >/dev/null

for K in $KLIST; do
    echo "Running ADMIXTURE for K=$K"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "admixture" "fit" "$K" "cv_log_K${K}.out" \
            admixture --cv -j "$THREADS" "${PREFIX}_ALL.pruned.bed" "$K"
    else
        admixture --cv -j "$THREADS" "${PREFIX}_ALL.pruned.bed" "$K" | tee "cv_log_K${K}.out"
    fi
done

# 5. Extract CV results for easy viewing
grep "CV error" cv_log_K*.out | sed 's/cv_log_K\(.*\).out:CV error (K=\(.*\)): \(.*\)/\2 \3/' | sort -n > cv_results.txt
popd >/dev/null

echo "Pipeline complete. Results and logs are in dump/admixture/"