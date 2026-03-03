#!/bin/bash

OUT_DIR=""
PREFIX=""
SAMPLE_INFO="1000Genomes/igsr_samples.tsv"
CHR_START="1"
CHR_END="22"
SKIP_MERGE="0"
RUN_CHR_PROCESS="0"
MAF="0.01"
LD_WINDOW="50"
LD_STEP="10"
LD_R2="0.1"
METHOD=""
KLIST=""
THREADS=""
FASTSTRUCTURE_BIN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --sample-info) SAMPLE_INFO="$2"; shift 2 ;;
        --chr-start) CHR_START="$2"; shift 2 ;;
        --chr-end) CHR_END="$2"; shift 2 ;;
        --skip-merge) SKIP_MERGE="1"; shift ;;
        --run-chr-process) RUN_CHR_PROCESS="1"; shift ;;
        --maf) MAF="$2"; shift 2 ;;
        --ld-window) LD_WINDOW="$2"; shift 2 ;;
        --ld-step) LD_STEP="$2"; shift 2 ;;
        --ld-r2) LD_R2="$2"; shift 2 ;;
        --method) METHOD="$2"; shift 2 ;;
        --klist) KLIST="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --faststructure-bin) FASTSTRUCTURE_BIN="$2"; shift 2 ;;
        --help)
            echo "Usage: bash preprocess.sh --out-dir DIR --prefix NAME [--sample-info FILE] [--chr-start N] [--chr-end N] [--skip-merge] [--run-chr-process] [--maf V] [--ld-window N] [--ld-step N] [--ld-r2 V] [--method structure --klist \"2 3\" --threads N --faststructure-bin PATH]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 2
            ;;
    esac
done

if [[ -z "$OUT_DIR" || -z "$PREFIX" ]]; then
    echo "OUT_DIR and PREFIX are required"
    exit 2
fi

mkdir -p "$OUT_DIR"

grep "1000 Genomes phase 3 release" "$SAMPLE_INFO" | awk '{print $1 "\t" $1}' > "$OUT_DIR/${PREFIX}_samples.txt"

if [[ "$RUN_CHR_PROCESS" == "1" ]]; then
    for chr_number in $(seq "$CHR_START" "$CHR_END"); do
        vcf="1000Genomes/1000G_chr${chr_number}_pruned.vcf.gz"

        echo "Processing Chromosome $chr_number..."

        plink --vcf "$vcf" --maf "$MAF" --keep "$OUT_DIR/${PREFIX}_samples.txt" --double-id --make-bed --out "$OUT_DIR/${PREFIX}_chr${chr_number}"
        plink --bfile "$OUT_DIR/${PREFIX}_chr${chr_number}" --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" --out "$OUT_DIR/tmp_chr${chr_number}"
        plink --bfile "$OUT_DIR/${PREFIX}_chr${chr_number}" --extract "$OUT_DIR/tmp_chr${chr_number}.prune.in" --make-bed --out "$OUT_DIR/${PREFIX}_chr${chr_number}.pruned"
    done
fi

echo "Creating merge list..."
MERGE_LIST="$OUT_DIR/mergelist.txt"
rm -f "$MERGE_LIST"
for i in $(seq "$(( CHR_START + 1 ))" "$CHR_END"); do
    echo "$OUT_DIR/${PREFIX}_chr${i}.pruned" >> "$MERGE_LIST"
done

if [[ "$SKIP_MERGE" != "1" ]]; then
    plink --bfile "$OUT_DIR/${PREFIX}_chr1.pruned" --merge-list "$MERGE_LIST" --make-bed --out "$OUT_DIR/${PREFIX}_ALL.pruned"
fi

if [[ "$METHOD" == "structure" ]]; then
    echo "Running STRUCTURE_THREADER"
    echo "Using fastStructure binary: $FASTSTRUCTURE_BIN"
    echo "K list: $KLIST"
    echo "Threads: $THREADS"
fi
