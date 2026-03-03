#!/bin/bash

# --- Setup ---
CONFIG_FILE="${PIPELINE_CONFIG:-pipeline.conf}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

# Create the specific output directory
OUT_DIR="dump/structure"
mkdir -p "$OUT_DIR"
sample_info="${SAMPLE_INFO:-1000Genomes/igsr_samples.tsv}"
PREFIX="${PREFIX:-${STRUCTURE_PREFIX:-structure}}"
THREADS="${THREADS:-10}"
KLIST="${KLIST:-2 3 4 5 6 7 8 9 10}"
# Generated from sample metadata (see analysis.ipynb, first code cell)
POP_FILE_DEFAULT="$OUT_DIR/${PREFIX}_pop_file.tsv"
POP_FILE="${POP_FILE:-$POP_FILE_DEFAULT}"
FASTSTRUCTURE_BIN="${FASTSTRUCTURE_BIN:-$HOME/.local/bin/fastStructure}"
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# shellcheck disable=SC1091
source "$SCRIPT_DIR/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi


# Preprocess
if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    benchmark_run "$METRICS_FILE" "faststructure" "preprocess" "NA" "$OUT_DIR/preprocess.log" \
        bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}" --method structure --klist "$KLIST" --threads "$THREADS" --faststructure-bin "$FASTSTRUCTURE_BIN"
else
    bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}" --method structure --klist "$KLIST" --threads "$THREADS" --faststructure-bin "$FASTSTRUCTURE_BIN"
fi


# 4. Run STRUCTURE / fastStructure via threader
# We run this from the base directory so the paths stay consistent
if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    benchmark_run "$METRICS_FILE" "faststructure" "fit" "all" "$OUT_DIR/structure_threader.log" \
        structure_threader run -Klist $KLIST -R 1 -i "$OUT_DIR/${PREFIX}_ALL.pruned.bed" -o "$OUT_DIR/${PREFIX}_ALL" -t "$THREADS" --pop "$POP_FILE" -fs "$FASTSTRUCTURE_BIN"
else
    structure_threader run -Klist $KLIST -R 1 -i "$OUT_DIR/${PREFIX}_ALL.pruned.bed" -o "$OUT_DIR/${PREFIX}_ALL" -t "$THREADS" --pop "$POP_FILE" -fs "$FASTSTRUCTURE_BIN"
fi

# 5. K_selection summary for fastStructure
if command -v chooseK.py >/dev/null 2>&1; then
    echo "Running chooseK.py"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "faststructure" "K_selection" "all" "$OUT_DIR/best_K.txt" \
            chooseK.py --input="$OUT_DIR/${PREFIX}_ALL/fS_run"
    else
        chooseK.py --input="$OUT_DIR/${PREFIX}_ALL/fS_run" | tee "$OUT_DIR/best_K.txt"
    fi
else
    echo "chooseK.py not found. Skipping model selection summary."
fi

echo "Pipeline complete. All files are in $OUT_DIR/"