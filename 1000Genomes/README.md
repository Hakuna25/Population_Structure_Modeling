## 1000Genomes Phase 3 — LD-Pruned Data

Fetch the LD-pruned 1000 Genomes Phase 3 data from the course DataHub directory:
`~/public/1000Genomes/`

You will need the pruned VCF files for chromosomes 1–22 (one file per chromosome) ```1000G_chr{1..22}_pruned.vcf.gz```, the sample metadata file `igsr_samples.tsv` and the original VCF file for chromosome 20 `ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` for this project.

### Minimal required files

```text
1000Genomes/
├── igsr_samples.tsv
├── 1000G_chr1_pruned.vcf.gz
├── 1000G_chr2_pruned.vcf.gz
├── ...
├── 1000G_chr22_pruned.vcf.gz
└── ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
```


### Reference genome

The reference genome used by `bcftools norm` is:
- `human_g1k_v37.fasta` (plus index files) under `1000Genomes/reference/`

If it is missing, `preprocess.sh` will automatically:
- Download `human_g1k_v37.fasta.gz` from the 1000 Genomes FTP.
- Apply the htslib-recommended `truncate` fix when needed.
- Decompress to a plain `.fasta` file.
- Rebuild the `.fai` index locally with `samtools faidx`.
The path to the FASTA is controlled by `REF_FA` in `pipeline.conf`.