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
RANDOM_SEED="${RANDOM_SEED:-43}"
# Per-individual population labels (see analysis.ipynb, Cell 3)
# Use --ind instead of --pop because FAM samples are not contiguously grouped by population
IND_FILE_DEFAULT="$OUT_DIR/${PREFIX}_ind_file.tsv"
IND_FILE="${IND_FILE:-$IND_FILE_DEFAULT}"
FASTSTRUCTURE_BIN="${FASTSTRUCTURE_BIN:-$HOME/.local/bin/fastStructure}"
STRUCTURE_THREADER_BIN="${STRUCTURE_THREADER_BIN:-$(command -v structure_threader || true)}"
BENCHMARK_ENABLED="${BENCHMARK_ENABLED:-1}"
METRICS_FILE="${METRICS_FILE:-$OUT_DIR/metrics.tsv}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [[ -z "$STRUCTURE_THREADER_BIN" && -x "$HOME/.local/bin/structure_threader" ]]; then
    STRUCTURE_THREADER_BIN="$HOME/.local/bin/structure_threader"
fi

# Convert to absolute path
METRICS_FILE="$(cd "$(dirname "$METRICS_FILE")" && pwd)/$(basename "$METRICS_FILE")"

# shellcheck disable=SC1091
source "$SCRIPT_DIR/benchmark.sh"

if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    init_metrics_file "$METRICS_FILE"
fi


# Preprocess
if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
    benchmark_run "$METRICS_FILE" "faststructure" "preprocess" "NA" "$OUT_DIR/preprocess.log" \
        bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
else
    bash "$SCRIPT_DIR/preprocess.sh" --out-dir "$OUT_DIR" --prefix "$PREFIX" --sample-info "$sample_info" --chr-start "${CHR_START:-1}" --chr-end "${CHR_END:-22}"
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

# Run STRUCTURE / fastStructure via threader
# We run this from the base directory so the paths stay consistent
for K in $KLIST; do
    echo "Running fastStructure for K=$K"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "faststructure" "fit" "$K" "$OUT_DIR/structure_threader_K${K}.log" \
            "$STRUCTURE_THREADER_BIN" run -Klist "$K" -R 1 -i "$OUT_DIR/${PREFIX}_ALL.pruned.bed" -o "$OUT_DIR/${PREFIX}_ALL" -t "$THREADS" --ind "$IND_FILE" --seed "$RANDOM_SEED" -fs "$FASTSTRUCTURE_BIN"
    else
        "$STRUCTURE_THREADER_BIN" run -Klist "$K" -R 1 -i "$OUT_DIR/${PREFIX}_ALL.pruned.bed" -o "$OUT_DIR/${PREFIX}_ALL" -t "$THREADS" --ind "$IND_FILE" --seed "$RANDOM_SEED" -fs "$FASTSTRUCTURE_BIN"
    fi
done

#  K_selection summary for fastStructure
if command -v chooseK.py >/dev/null 2>&1; then
    echo "Running chooseK.py"
    if [[ "$BENCHMARK_ENABLED" == "1" ]]; then
        benchmark_run "$METRICS_FILE" "faststructure" "K_selection" "all" "$OUT_DIR/best_K.txt" \
            chooseK.py --input="$OUT_DIR/${PREFIX}_ALL/fS_run_K"
    else
        chooseK.py --input="$OUT_DIR/${PREFIX}_ALL/fS_run_K" | tee "$OUT_DIR/best_K.txt"
    fi
else
    echo "chooseK.py not found. Skipping model selection summary."
fi

echo "Pipeline complete. All files are in $OUT_DIR/"