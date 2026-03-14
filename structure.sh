#!/bin/bash

# --- Setup ---
CONFIG_FILE="${PIPELINE_CONFIG:-pipeline.conf}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

# Create the specific output directory
OUT_DIR="${OUT_DIR:-dump/structure}"
mkdir -p "$OUT_DIR"
PREFIX="${PREFIX:-${STRUCTURE_PREFIX:-structure}}"
THREADS="${THREADS:-1}"
REPLICATES="${REPLICATES:-1}"
KLIST="${KLIST:-2 3 4 5 6 7 8 9 10}"
RANDOM_SEED="${RANDOM_SEED:-43}"
# Per-individual population labels (see analysis.ipynb, Cell 3)
# Use --ind instead of --pop because FAM samples are not contiguously grouped by population
IND_FILE_DEFAULT="$OUT_DIR/${PREFIX}_ind_file.tsv"
IND_FILE="${IND_FILE:-$IND_FILE_DEFAULT}"
FASTSTRUCTURE_BIN="${FASTSTRUCTURE_BIN:-$HOME/.local/bin/fastStructure}"
STRUCTURE_THREADER_BIN="${STRUCTURE_THREADER_BIN:-$(command -v structure_threader || true)}"
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
DATA_PREPROCESS_DIR="${DATA_PREPROCESS_DIR:-dump/common}"
DATA_PREPROCESS_PREFIX="${DATA_PREPROCESS_PREFIX:-common}"
INPUT_BED="${INPUT_BED:-$DATA_PREPROCESS_DIR/${DATA_PREPROCESS_PREFIX}_ALL.pruned.bed}"

if [[ -z "$STRUCTURE_THREADER_BIN" && -x "$HOME/.local/bin/structure_threader" ]]; then
    STRUCTURE_THREADER_BIN="$HOME/.local/bin/structure_threader"
fi

# Convert to absolute path
METRICS_FILE="$(cd "$(dirname "$METRICS_FILE")" && pwd)/$(basename "$METRICS_FILE")"

if [[ ! -f "$INPUT_BED" ]]; then
    echo "ERROR: shared preprocessed input not found: $INPUT_BED"
    echo "Please run preprocess.sh first, for example from analysis.ipynb."
    exit 1
fi
INPUT_BED="$(cd "$(dirname "$INPUT_BED")" && pwd)/$(basename "$INPUT_BED")"

# shellcheck disable=SC1091
source "$(cd "$(dirname "$0")" && pwd)/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi

# Validate ind file before running
if [[ ! -f "$IND_FILE" ]]; then
    echo "ERROR: Individual population file not found: $IND_FILE"
    echo "Please generate it first (see analysis.ipynb, Cell 3)."
    exit 1
fi

if [[ -z "$STRUCTURE_THREADER_BIN" || ! -x "$STRUCTURE_THREADER_BIN" ]]; then
    echo "ERROR: structure_threader not found or not executable."
    echo "Please install it (e.g., pip install structure_threader --user) or set STRUCTURE_THREADER_BIN in pipeline.conf."
    exit 1
fi

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer, got: $THREADS"
    exit 2
fi

RESULT_DIR="$OUT_DIR/${PREFIX}_ALL"
if [[ -e "$RESULT_DIR" && ! -d "$RESULT_DIR" ]]; then
    echo "ERROR: fastStructure result path exists and is not a directory: $RESULT_DIR"
    exit 2
fi
mkdir -p "$RESULT_DIR"

run_one_k() {
    local K="$1"
    local fit_exit_code=0
    local meanq_file
    local -a run_cmd

    meanq_file="$RESULT_DIR/fS_run_K.${K}.meanQ"
    if [[ -f "$meanq_file" ]]; then
        if [[ "$FORCE_RERUN" == "1" ]]; then
            echo "FORCE_RERUN=1: rerun existing fastStructure meanQ for K=$K: $meanq_file"
        else
            echo "Skipping fastStructure for K=$K because $meanq_file already exists"
            return 0
        fi
    fi

    # Disable structure_threader's internal bestK tests and plotting.
    # This script launches one K per invocation, so those global post-processing
    # steps can incorrectly fail even when fastStructure itself produced meanQ.
    run_cmd=(
        "$STRUCTURE_THREADER_BIN" run
        -Klist "$K"
        -R "$REPLICATES"
        -i "$INPUT_BED"
        -o "$RESULT_DIR"
        -t 1
        --ind "$IND_FILE"
        --seed "$RANDOM_SEED"
        --no_tests True
        --no_plots True
        -fs "$FASTSTRUCTURE_BIN"
    )

    echo "Running fastStructure for K=$K"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "faststructure" "fit" "$K" "$OUT_DIR/structure_threader_K${K}.log" \
            "${run_cmd[@]}"
        fit_exit_code=$?
    else
        "${run_cmd[@]}"
        fit_exit_code=$?
    fi

    if [[ "$fit_exit_code" -ne 0 ]]; then
        echo "ERROR: fastStructure failed at K=$K with exit code $fit_exit_code. See $OUT_DIR/structure_threader_K${K}.log"
    fi
    return "$fit_exit_code"
}

# Run up to THREADS single-threaded jobs in parallel.
overall_exit_code=0
declare -A active_pids=()

reap_one_job() {
    local finished_pid

    if wait -n -p finished_pid; then
        unset "active_pids[$finished_pid]"
        return 0
    fi

    if [[ -n "$finished_pid" ]]; then
        unset "active_pids[$finished_pid]"
        overall_exit_code=1
        return 0
    fi

    return 1
}

for K in $KLIST; do
    while (( ${#active_pids[@]} >= THREADS )); do
        reap_one_job || break
    done

    run_one_k "$K" &
    active_pids[$!]=1
done

while (( ${#active_pids[@]} > 0 )); do
    reap_one_job || break
done

if [[ "$overall_exit_code" -ne 0 ]]; then
    exit "$overall_exit_code"
fi

#  K_selection summary for fastStructure
if command -v chooseK.py >/dev/null 2>&1; then
    if [[ -f "$OUT_DIR/best_K.txt" ]]; then
        echo "Skipping chooseK.py because $OUT_DIR/best_K.txt already exists"
    else
        echo "Running chooseK.py"
        if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
            benchmark_run "$METRICS_FILE" "faststructure" "K_selection" "all" "$OUT_DIR/best_K.txt" \
                chooseK.py --input="$RESULT_DIR/fS_run_K"
        else
            chooseK.py --input="$RESULT_DIR/fS_run_K" | tee "$OUT_DIR/best_K.txt"
        fi
    fi
else
    echo "chooseK.py not found. Skipping model selection summary."
fi

echo "Pipeline complete. All files are in $OUT_DIR/"