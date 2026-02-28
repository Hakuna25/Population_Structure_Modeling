#!/bin/bash

# --- Setup ---
mkdir -p dump/admixture
vcf_files='phase3_vcf_files.txt'
sample_info='1000Genomes/igsr_samples.tsv'
PREFIX='admixture'

readarray -t vcfs < "$vcf_files"

# 1. Extract samples
cat "$sample_info" | grep "1000 Genomes phase 3 release" | awk '{print $1 "\t" $1}' > ${PREFIX}_samples.txt

# 2. Process Chromosomes (Uncomment the loop below to run preprocessing)
# for vcf in "${vcfs[@]}"; do
#     file_raw=$(basename "$vcf" .vcf.gz)
#     chr_number=$(echo "$file_raw" | sed 's/.*chr\([0-9]\+\).*/\1/')
#     
#     plink --vcf "$vcf" --maf 0.01 --keep ${PREFIX}_samples.txt --double-id --make-bed --out ./dump/${PREFIX}_chr${chr_number}
#     plink --bfile ./dump/${PREFIX}_chr${chr_number} --indep-pairwise 50 10 0.1
#     plink --bfile ./dump/${PREFIX}_chr${chr_number} --extract plink.prune.in --make-bed --out ./dump/${PREFIX}_chr${chr_number}.pruned
# done

# 3. Merge Files
echo "Creating merge list..."
rm -f admixture_filenames.txt
for i in $(seq 2 22); do
    echo "dump/${PREFIX}_chr${i}.pruned" >> admixture_filenames.txt
done

plink --bfile dump/${PREFIX}_chr1.pruned --merge-list admixture_filenames.txt --make-bed --out ./dump/${PREFIX}_ALL.pruned

# 4. Run ADMIXTURE
pushd dump/admixture >/dev/null
# Copy the merged bed/bim/fam into the admixture folder for the tool to work
ln -s ../${PREFIX}_ALL.pruned.* .

for K in {7..10}; do
    echo "Running ADMIXTURE for K=$K"
    admixture --cv ${PREFIX}_ALL.pruned.bed $K | tee log_K${K}.out
done

# 5. Extract CV results for easy viewing
grep "CV error" log_K*.out | sed 's/log_K\(.*\).out:CV error (K=\(.*\)): \(.*\)/\2 \3/' | sort -n > cv_results.txt
popd >/dev/null

echo "Pipeline complete. Results and logs are in dump/admixture/"