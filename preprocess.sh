#!/bin/bash

set -euo pipefail

CONFIG_FILE="${PIPELINE_CONFIG:-pipeline.conf}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

OUT_DIR=""
PREFIX=""
SAMPLE_INFO="${SAMPLE_INFO:-1000Genomes/igsr_samples.tsv}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-$(command -v bcftools || true)}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
REF_FA="${REF_FA:-1000Genomes/reference/human_g1k_v37.fasta}"
REF_FAI="${REF_FAI:-${REF_FA}.fai}"
REF_FA_URL="${REF_FA_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz}"
REF_FA_GZ_TRUNCATE_SIZE="${REF_FA_GZ_TRUNCATE_SIZE:-891946027}"
PREPROCESS_THREADS="${PREPROCESS_THREADS:-8}"
CHR_START="${CHR_START:-1}"
CHR_END="${CHR_END:-22}"
RUN_CHR_PROCESS="${RUN_CHR_PROCESS:-1}"
SKIP_MERGE="0"
SKIP_LD_PRUNE="0"
SKIP_NORMALIZE="0"
VCF_TEMPLATE="${VCF_TEMPLATE:-1000Genomes/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz}"
OUTPUT_PRUNED="${OUTPUT_PRUNED:-}"
MAF="${MAF:-0.01}"
LD_WINDOW="${LD_WINDOW:-50}"
LD_STEP="${LD_STEP:-10}"
LD_R2="${LD_R2:-0.1}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --sample-info) SAMPLE_INFO="$2"; shift 2 ;;
        --ref-fa) REF_FA="$2"; shift 2 ;;
        --threads) PREPROCESS_THREADS="$2"; shift 2 ;;
        --chr-start) CHR_START="$2"; shift 2 ;;
        --chr-end) CHR_END="$2"; shift 2 ;;
        --skip-chr-process) RUN_CHR_PROCESS="0"; shift ;;
        --skip-merge) SKIP_MERGE="1"; shift ;;
        --skip-ld-prune) SKIP_LD_PRUNE="1"; shift ;;
        --run-ld-prune) SKIP_LD_PRUNE="0"; shift ;;
        --skip-normalize) SKIP_NORMALIZE="1"; shift ;;
        --run-normalize) SKIP_NORMALIZE="0"; shift ;;
        --run-chr-process) RUN_CHR_PROCESS="1"; shift ;;
        --vcf-template) VCF_TEMPLATE="$2"; shift 2 ;;
        --output-pruned) OUTPUT_PRUNED="1"; shift ;;
        --output-unpruned) OUTPUT_PRUNED="0"; shift ;;
        --maf) MAF="$2"; shift 2 ;;
        --skip-maf) MAF=""; shift ;;
        --ld-window) LD_WINDOW="$2"; shift 2 ;;
        --ld-step) LD_STEP="$2"; shift 2 ;;
        --ld-r2) LD_R2="$2"; shift 2 ;;
        --help)
            echo "Usage: bash preprocess.sh --out-dir DIR --prefix NAME [--sample-info FILE] [--ref-fa FILE] [--threads N] [--chr-start N] [--chr-end N] [--run-chr-process] [--skip-chr-process] [--skip-merge] [--skip-ld-prune] [--skip-normalize] [--skip-maf] [--vcf-template TEMPLATE] [--output-pruned|--output-unpruned] [--maf V] [--ld-window N] [--ld-step N] [--ld-r2 V]"
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

if [[ -z "$OUTPUT_PRUNED" ]]; then
    if [[ "$SKIP_LD_PRUNE" == "1" ]]; then
        OUTPUT_PRUNED="0"
    else
        OUTPUT_PRUNED="1"
    fi
fi

MERGE_INPUT_SUFFIX=""
if [[ "$OUTPUT_PRUNED" == "1" ]]; then
    MERGE_INPUT_SUFFIX=".pruned"
fi

if [[ "$CHR_START" == "$CHR_END" ]]; then
    DATASET_LABEL="chr${CHR_START}"
else
    DATASET_LABEL="ALL"
fi
MERGED_OUTPUT_PREFIX="$OUT_DIR/${PREFIX}_${DATASET_LABEL}${MERGE_INPUT_SUFFIX}"

if ! command -v plink >/dev/null 2>&1; then
    echo "plink not found in PATH"
    exit 127
fi

if [[ "$SKIP_NORMALIZE" != "1" ]]; then
    if [[ -z "$BCFTOOLS_BIN" || ! -x "$BCFTOOLS_BIN" ]]; then
        CONDA_BIN="${CONDA_EXE:-$(command -v conda || true)}"
        if [[ -n "$CONDA_BIN" && -x "$CONDA_BIN" ]]; then
            CONDA_ROOT="$(dirname "$(dirname "$CONDA_BIN")")"
            BCFTOOLS_BIN="$CONDA_ROOT/envs/bcftools_tools/bin/bcftools"
        fi
    fi

    if [[ -z "$BCFTOOLS_BIN" || ! -x "$BCFTOOLS_BIN" ]]; then
        echo "bcftools not found or not executable."
        echo "Please activate an environment with bcftools in PATH, or set BCFTOOLS_BIN in pipeline.conf."
        exit 127
    fi

    if [[ -z "$SAMTOOLS_BIN" || ! -x "$SAMTOOLS_BIN" ]]; then
        CONDA_BIN="${CONDA_EXE:-$(command -v conda || true)}"
        if [[ -n "$CONDA_BIN" && -x "$CONDA_BIN" ]]; then
            CONDA_ROOT="$(dirname "$(dirname "$CONDA_BIN")")"
            SAMTOOLS_BIN="$CONDA_ROOT/envs/bcftools_tools/bin/samtools"
        fi
    fi

    if [[ -z "$SAMTOOLS_BIN" || ! -x "$SAMTOOLS_BIN" ]]; then
        echo "samtools not found or not executable."
        echo "Please activate an environment with samtools in PATH, or set SAMTOOLS_BIN in pipeline.conf."
        exit 127
    fi

    if ! command -v wget >/dev/null 2>&1; then
        echo "wget not found in PATH"
        exit 127
    fi

    if ! command -v gunzip >/dev/null 2>&1; then
        echo "gunzip not found in PATH"
        exit 127
    fi

    if ! command -v gzip >/dev/null 2>&1; then
        echo "gzip not found in PATH"
        exit 127
    fi

    mkdir -p "$(dirname "$REF_FA")"
    REF_FA_GZ="${REF_FA}.gz"
    if [[ ! -f "$REF_FA_GZ" ]]; then
        echo "Downloading $(basename "$REF_FA_GZ")..."
        wget -O "$REF_FA_GZ" "$REF_FA_URL"
    fi

    if [[ ! -s "$REF_FA" ]]; then
        if ! gzip -t "$REF_FA_GZ" 2>/dev/null; then
            REF_FA_GZ_SIZE="$(stat -c%s "$REF_FA_GZ")"
            if [[ "$REF_FA_GZ_SIZE" -gt "$REF_FA_GZ_TRUNCATE_SIZE" ]]; then
                if ! command -v truncate >/dev/null 2>&1; then
                    echo "truncate not found in PATH"
                    exit 127
                fi
                echo "Truncating $(basename "$REF_FA_GZ") to $REF_FA_GZ_TRUNCATE_SIZE bytes..."
                truncate -s "$REF_FA_GZ_TRUNCATE_SIZE" "$REF_FA_GZ"
                gzip -t "$REF_FA_GZ"
            else
                echo "Reference gzip is corrupt or incomplete: $REF_FA_GZ"
                exit 1
            fi
        fi

        echo "Decompressing $(basename "$REF_FA_GZ")..."
        gunzip -c "$REF_FA_GZ" > "${REF_FA}.tmp"
        mv "${REF_FA}.tmp" "$REF_FA"
    fi

    rm -f "$REF_FAI"
    "$SAMTOOLS_BIN" faidx "$REF_FA"
fi

awk '/1000 Genomes phase 3 release/ {print $1 "\t" $1}' "$SAMPLE_INFO" > "$OUT_DIR/${PREFIX}_samples.txt"

if [[ ! -s "$OUT_DIR/${PREFIX}_samples.txt" ]]; then
    echo "No samples extracted from $SAMPLE_INFO"
    exit 2
fi

process_chromosome() {
    local chr_number="$1"
    local vcf="${VCF_TEMPLATE//\{chr\}/$chr_number}"
    local dedup_vcf="$OUT_DIR/${PREFIX}_chr${chr_number}_dedup.vcf.gz"
    local chr_prefix="$OUT_DIR/${PREFIX}_chr${chr_number}"
    local final_chr_prefix="${chr_prefix}${MERGE_INPUT_SUFFIX}"
    local bim_path="${chr_prefix}.bim"
    local bim_tmp_path="${chr_prefix}_bim_tmp"
    local vcf_for_plink="$vcf"

    echo "Processing Chromosome $chr_number..."

    if [[ ! -f "$vcf" ]]; then
        echo "Missing VCF for chromosome $chr_number: $vcf"
        exit 2
    fi

    if [[ "$SKIP_NORMALIZE" == "1" ]]; then
        echo "Skipping VCF normalization for chromosome $chr_number."
    else
        "$BCFTOOLS_BIN" norm --threads 1 -f "$REF_FA" -m -any -d all -Oz -o "$dedup_vcf" "$vcf"
        vcf_for_plink="$dedup_vcf"
    fi

    plink_args=(
        --threads 1
        --vcf "$vcf_for_plink"
        --keep "$OUT_DIR/${PREFIX}_samples.txt"
        --double-id
        --make-bed
        --out "$chr_prefix"
        --keep-allele-order
    )
    if [[ -n "$MAF" ]]; then
        plink_args+=(--maf "$MAF")
        plink_args+=(--snps-only just-acgt)
    else
        plink_args+=(--maf 0.01)
    fi
    plink "${plink_args[@]}"

    awk '
        BEGIN {
            OFS = "\t"
        }
        {
            variant_id = $1 ":" $4 ":" $5 ":" $6
            seen_count[variant_id] += 1
            if (seen_count[variant_id] > 1) {
                print "ERROR: duplicate variant key found after dedup: " variant_id > "/dev/stderr"
                exit 1
            }
            $2 = variant_id
            print
        }
    ' "$bim_path" > "$bim_tmp_path"
    mv "$bim_tmp_path" "$bim_path"

    if [[ "$SKIP_LD_PRUNE" == "1" ]]; then
        echo "Skipping LD pruning for chromosome $chr_number."
        if [[ "$MERGE_INPUT_SUFFIX" != "" ]]; then
            cp "${chr_prefix}.bed" "${final_chr_prefix}.bed"
            cp "${chr_prefix}.bim" "${final_chr_prefix}.bim"
            cp "${chr_prefix}.fam" "${final_chr_prefix}.fam"
        fi
    else
        plink --threads 1 --bfile "$chr_prefix" --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" --out "$OUT_DIR/tmp_chr${chr_number}"
        plink --threads 1 --bfile "$chr_prefix" --extract "$OUT_DIR/tmp_chr${chr_number}.prune.in" --make-bed --out "$final_chr_prefix"
    fi
}

if [[ "$RUN_CHR_PROCESS" == "1" ]]; then
    active_jobs=0
    overall_exit_code=0
    echo "Running chromosome preprocessing with up to $PREPROCESS_THREADS parallel jobs."
    for chr_number in $(seq "$CHR_START" "$CHR_END"); do
        process_chromosome "$chr_number" &
        ((active_jobs += 1))

        if (( active_jobs >= PREPROCESS_THREADS )); then
            if ! wait -n; then
                overall_exit_code=1
            fi
            active_jobs=$((active_jobs - 1))
        fi
    done

    while (( active_jobs > 0 )); do
        if ! wait -n; then
            overall_exit_code=1
        fi
        active_jobs=$((active_jobs - 1))
    done

    if [[ "$overall_exit_code" -ne 0 ]]; then
        echo "ERROR: chromosome preprocessing failed."
        exit "$overall_exit_code"
    fi
else
    echo "Skipping chromosome-level preprocessing because RUN_CHR_PROCESS=0."
fi

MERGE_LIST="$OUT_DIR/mergelist.txt"
rm -f "$MERGE_LIST"
if [[ "$SKIP_MERGE" != "1" ]]; then
    BASE_PREFIX="$OUT_DIR/${PREFIX}_chr${CHR_START}${MERGE_INPUT_SUFFIX}"
    if [[ ! -f "${BASE_PREFIX}.bed" || ! -f "${BASE_PREFIX}.bim" || ! -f "${BASE_PREFIX}.fam" ]]; then
        echo "Missing required chromosome PLINK files for merge: ${BASE_PREFIX}.*"
        exit 2
    fi

    if [[ "$CHR_START" == "$CHR_END" ]]; then
        if [[ "$BASE_PREFIX" != "$MERGED_OUTPUT_PREFIX" ]]; then
            cp "${BASE_PREFIX}.bed" "${MERGED_OUTPUT_PREFIX}.bed"
            cp "${BASE_PREFIX}.bim" "${MERGED_OUTPUT_PREFIX}.bim"
            cp "${BASE_PREFIX}.fam" "${MERGED_OUTPUT_PREFIX}.fam"
        fi
    else
        echo "Creating merge list..."
        for i in $(seq "$(( CHR_START + 1 ))" "$CHR_END"); do
            echo "$OUT_DIR/${PREFIX}_chr${i}${MERGE_INPUT_SUFFIX}" >> "$MERGE_LIST"
        done
        plink --threads "$PREPROCESS_THREADS" --bfile "$BASE_PREFIX" --merge-list "$MERGE_LIST" --make-bed --out "$MERGED_OUTPUT_PREFIX"
    fi
fi
