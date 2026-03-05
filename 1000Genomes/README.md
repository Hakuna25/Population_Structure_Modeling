## 1000Genomes Phase 3 — LD-Pruned Data

Fetch the LD-pruned 1000 Genomes Phase 3 data from the course DataHub directory:
`~/public/1000Genomes/`

You only need the pruned VCF files for chromosomes 1–22 (one file per chromosome) ```1000G_chr{1..22}_pruned.vcf.gz``` and the sample metadata file `igsr_samples.tsv` for this project.

Place all files directly under this folder:
```
1000Genomes/
    igsr_samples.tsv
    1000G_chr1_pruned.vcf.gz
    1000G_chr2_pruned.vcf.gz
    ...
    1000G_chr22_pruned.vcf.gz
```

You can also use the original VCF files if you want to test the full pipeline. In this situation, the additional required files include:
```1000Genomes/
    ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    ...
    ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```
You also need to change the `RUN_CHR_PROCESS="0"` in `preprocess.sh` to `RUN_CHR_PROCESS="1"` to enable pruning from raw files.