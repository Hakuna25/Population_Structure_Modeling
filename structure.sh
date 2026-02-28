#!/bin/bash

# --- Setup ---
# Create the specific output directory
mkdir -p dump/structure
vcf_files='phase3_vcf_files.txt'
sample_info='1000Genomes/igsr_samples.tsv'
PREFIX='structure'

readarray -t vcfs < "$vcf_files"

# 1. Extract samples (Saving to the target directory)
cat "$sample_info" | grep "1000 Genomes phase 3 release" | awk '{print $1 "\t" $1}' > dump/structure/${PREFIX}_samples.txt

# 2. Process Chromosomes
# for vcf in "${vcfs[@]}"; do
#     file_raw=$(basename "$vcf" .vcf.gz)
#     chr_number=$(echo "$file_raw" | sed 's/.*chr\([0-9]\+\).*/\1/')
    
#     echo "Processing Chromosome $chr_number..."
    
#     # All intermediate files now go directly into dump/structure/
#     plink --vcf "$vcf" --maf 0.01 --keep dump/structure/${PREFIX}_samples.txt --double-id --make-bed --out dump/structure/${PREFIX}_chr${chr_number}
    
#     # LD Pruning
#     plink --bfile dump/structure/${PREFIX}_chr${chr_number} --indep-pairwise 50 10 0.1 --out dump/structure/tmp_chr${chr_number}
#     plink --bfile dump/structure/${PREFIX}_chr${chr_number} --extract dump/structure/tmp_chr${chr_number}.prune.in --make-bed --out dump/structure/${PREFIX}_chr${chr_number}.pruned
# done

# 3. Merge Files
echo "Creating merge list..."
MERGE_LIST="dump/structure/mergelist.txt"
rm -f "$MERGE_LIST"

# Note: Start from 2 because chr1 is the reference for the merge
for i in $(seq 2 22); do
    echo "dump/structure/${PREFIX}_chr${i}.pruned" >> "$MERGE_LIST"
done

plink --bfile dump/structure/${PREFIX}_chr1.pruned --merge-list "$MERGE_LIST" --make-bed --out dump/structure/${PREFIX}_ALL.pruned

# 4. Run STRUCTURE / fastStructure via threader
# We run this from the base directory so the paths stay consistent
echo "Running STRUCTURE_THREADER"

structure_threader run -Klist 2 3 4 5 6 7 8 9 10 -R 1 -i dump/structure/${PREFIX}_ALL.pruned.bed -o dump/structure/${PREFIX}_ALL -t 10 --pop dump/structure/structure_pop_file.tsv -fs /home/xueqian/.local/bin/fastStructure

# 5. Handle "CV" results
# IMPORTANT: As discussed, fastStructure uses chooseK.py, not --cv. 
# If you used ADMIXTURE, the grep would work. For fastStructure:
# python /path/to/fastStructure/chooseK.py --input=dump/structure/results/${PREFIX}_ALL.pruned > dump/structure/best_K.txt

echo "Pipeline complete. All files are in dump/structure/"