#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="$ROOT_DIR/dump/test"
mkdir -p "$OUT_DIR"
LOG_FILE="$OUT_DIR/test.log"
: > "$LOG_FILE"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] Start test" | tee -a "$LOG_FILE"

for script in preprocess.sh admixture.sh structure.sh benchmark.sh test.sh; do
    bash -n "$ROOT_DIR/$script"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] bash -n passed: $script" | tee -a "$LOG_FILE"
done

TMP_DIR="$(mktemp -d "$OUT_DIR/tmp.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

if [[ ! -f "$ROOT_DIR/1000Genomes/igsr_samples.tsv" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] missing real sample metadata: 1000Genomes/igsr_samples.tsv" | tee -a "$LOG_FILE"
    exit 1
fi

PRE_DIR="$TMP_DIR/pre"
mkdir -p "$PRE_DIR"

for chr in 1 2 3; do
    for ext in bed bim fam; do
        src="$ROOT_DIR/dump/admixture/admixture_chr${chr}.pruned.${ext}"
        dst="$PRE_DIR/demo_chr${chr}.pruned.${ext}"
        if [[ ! -f "$src" ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] missing required file: $src" | tee -a "$LOG_FILE"
            exit 1
        fi
        cp "$src" "$dst"
    done
done

bash "$ROOT_DIR/preprocess.sh" \
    --out-dir "$PRE_DIR" \
    --prefix "demo" \
    --sample-info "$ROOT_DIR/1000Genomes/igsr_samples.tsv" \
    --chr-start 1 \
    --chr-end 3 \
    --skip-merge

MERGE_LIST="$PRE_DIR/mergelist.txt"
rm -f "$MERGE_LIST"
for i in $(seq 2 3); do
    echo "$PRE_DIR/demo_chr${i}.pruned" >> "$MERGE_LIST"
done

plink --bfile "$PRE_DIR/demo_chr1.pruned" --merge-list "$MERGE_LIST" --make-bed --out "$PRE_DIR/demo_ALL.pruned"

if [[ ! -f "$PRE_DIR/demo_samples.txt" || ! -f "$PRE_DIR/mergelist.txt" || ! -f "$PRE_DIR/demo_ALL.pruned.bed" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] preprocess example test failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] preview demo_samples.txt" | tee -a "$LOG_FILE"
head -n 5 "$PRE_DIR/demo_samples.txt" | tee -a "$LOG_FILE"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] preview mergelist.txt" | tee -a "$LOG_FILE"
cat "$PRE_DIR/mergelist.txt" | tee -a "$LOG_FILE"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] All tests passed" | tee -a "$LOG_FILE"
