# Comparing Global Ancestry Analysis Methods on 1000 Genomes Phase 3 Data

<!-- [Xueqian Bai](https://hakuna25.github.io/), [Jialin Wu](https://jlwu.cn/), [Yutong Liang](https://www.lyt0112.com/) -->
Xueqian Bai, Jialin Wu, Yutong Liang

## Introduction
Global ancestry analysis is a fundamental task in population genetics, aiming to infer the ancestral origins of individuals based on their genetic data. This analysis provides insights into human migration patterns, population structure, and evolutionary history. In this project, we focus on comparing two widely used model-based methods for global ancestry analysis: [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655) and [fastSTRUCTURE](https://doi.org/10.1534/genetics.114.164350). ADMIXTURE is a maximum-likelihood
admixture model that estimates individual ancestry proportions efficiently for large SNP datasets using numerical optimization. In contrast, fastSTRUCTURE adopts a Bayesian framework and applies variational inference to approximate posterior distributions of ancestry proportions. Both methods utilize genotype data to estimate the proportion of ancestry from different populations for each individual. By applying these methods to the 1000 Genomes Phase 3 dataset, we aim to evaluate their performance and compare the inferred ancestry proportions across different populations.

## Requirements
- Python 3.11
- Main environment: plink==1.90b7.7, admixture==1.3.1, faststructure==1.0
- Separate environment: bcftools

## Installation
```
conda env create -f env.yml
conda activate bio_tools
pip install structure_threader --user

# create a dedicated bcftools environment
conda env create -f env_bcftools.yml

# install official ADMIXTURE 1.3.1 binary
mkdir -p tools
cd tools
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.1.tar.gz
tar -xzf admixture_linux-1.3.1.tar.gz
cd ..

# configure pipeline binary paths in pipeline.conf
# ADMIXTURE_BIN="$PWD/tools/admixture_linux-1.3.1/admixture"
# BCFTOOLS_BIN="$HOME/anaconda3/envs/bcftools_tools/bin/bcftools"
# SAMTOOLS_BIN="$HOME/anaconda3/envs/bcftools_tools/bin/samtools"
# REF_FA="1000Genomes/reference/human_g1k_v37.fasta"
# PREPROCESS_THREADS="8"
# DATA_PREPROCESS_DIR="dump/common"
# DATA_PREPROCESS_PREFIX="common"
```
## 🚀 Quick Start (View Results)
If you want to skip the computation and dive straight into the findings:

Open `analysis.ipynb`: This notebook contains the pre-rendered experimental results, visualizations, and detailed analysis. It also includes extra experiment blocks for `no LD pruning` and `MAF=0.05`, alongside the original baseline workflow.

## 🧬 Data and assets
[1000 Genomes Phase 3 Data](https://www.nature.com/articles/nature15393) is a comprehensive release of the 1000 Genomes Project dataset, providing whole-genome sequencing–based variant calls for 2,504 individuals from 26 populations across five continental groups. In this project, we focus on autosomal biallelic SNPs and use an LD-pruned version of the autosomal data.

For reproducibility, please place all input 1000 Genomes files in `1000Genomes/`.
For detailed download instructions and expected filenames, see [1000Genomes/README.md](https://github.com/Hakuna25/Population_Structure_Modeling/blob/main/1000Genomes/README.md).

For `bcftools norm`, `preprocess.sh` will automatically download the GRCh37 reference FASTA archive into `1000Genomes/reference/`, apply the htslib-recommended `truncate` fix when needed, decompress it to a plain `.fasta` file, and rebuild the `.fai` index locally with `samtools faidx`.

Official reference files used by the script:
- [human_g1k_v37.fasta.gz](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)

## ⚙️ Reproducible Workflow
Follow these steps to reproduce the environment and the full analysis from scratch:
1. Create environment and install dependencies.
2. Place 1000 Genomes input files in `1000Genomes/`.
3. Edit `pipeline.conf` for config. 
4. Open `analysis.ipynb` and run the shared `preprocess.sh` step once.
5. Run `admixture.sh` or `structure.sh`; both scripts consume the shared `dump/common/common_ALL.pruned.*` outputs directly.
 We’ve already provided the preprocessed, merged 1000 Genomes Phase 3 chromosome dataset in `dump/`, so you can run the tools directly. If you’re interested, you can also reproduce the preprocessing steps using the snippets included in the code.

### Additional Notes
- `preprocess.sh` contains the shared preprocessing steps:
  - sample extraction
  - per-chromosome VCF normalization and deduplication with `bcftools norm -f "$REF_FA" -m -any -d exact`
  - PLINK conversion to chromosome-level binary files
  - variant ID normalization for downstream merge safety
  - LD pruning from the shared chromosome-level files
  - merge-list creation
  - merged PLINK generation
  - optional chromosome-stage reuse via `RUN_CHR_PROCESS=0`

- `bcftools` is resolved through `BCFTOOLS_BIN`. By default, the pipeline first checks `PATH`, then falls back to `<conda-root>/envs/bcftools_tools/bin/bcftools`.
- `samtools` is resolved through `SAMTOOLS_BIN`. By default, the pipeline first checks `PATH`, then falls back to `<conda-root>/envs/bcftools_tools/bin/samtools`.
- The reference FASTA path is controlled by `REF_FA` in `pipeline.conf`. If the FASTA is missing, `preprocess.sh` downloads `human_g1k_v37.fasta.gz`, repairs the archive with `truncate` when needed, decompresses it to `human_g1k_v37.fasta`, and rebuilds the `.fai` index locally with `samtools faidx`.
- Chromosome preprocessing parallelism is controlled by `PREPROCESS_THREADS` in `pipeline.conf` or `preprocess.sh --threads`; this sets how many chromosome jobs run concurrently.

- `bash test.sh` provides a quick test for the environment setup and method implementation. Test outputs are written to `dump/test/`.

## ⚡ Runtime & Memory Benchmarking
Both pipelines support built-in benchmarking (wall-clock runtime + peak memory usage).

- Enable/disable in `pipeline.conf`:
  - `BENCHMARK_ENABLED="1"` (default enabled)
- Output metric files:
  - `dump/admixture/metrics.tsv`
  - `dump/structure/metrics.tsv`

Each row in `metrics.tsv` contains:
- `method`: `admixture` or `faststructure`
- `step`: e.g. `preprocess`, `fit`, `K_selection`
- `K`: corresponding K value
- `elapsed_sec`: time in seconds
- `max_rss_kb`: peak resident memory (KB)

## 🌳 Repository structure
- `1000Genomes/`: contains the original input data, specifically the igsr_samples.tsv and the prunded VCF files.
- `dump/`:
    - `common/`: shared outputs from `preprocess.sh`.
    - `admixture/`: the primary workspace for Admixture. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
    - `structure/`: the primary workspace for fastStructure. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
- Root scripts: `preprocess.sh` builds the shared analysis input once, `admixture.sh` and `structure.sh` run the two methods, and `analysis.ipynb` is the main entry point for orchestration and results analysis.

## 📊 Current Results
The detailed current results are shown in `analysis.ipynb`, here is a brief summary:

- For ADMIXTURE, the cross-validation (CV) error drops quickly for small K=5 as the best clustering resolution for this dataset. For fastSTRUCTURE, the marginal likelihood shows a clear peak and elbow at K=7, indicating a similar estimate of the underlying number of clusters.

- Both methods successfully identify the major continental clusters corresponding to the 1000 Genomes super-populations (AFR, AMR, EAS, EUR, SAS). fastSTRUCTURE tends to produce **cleaner** looking plots with less background noise.

- fastSTRUCTURE is substantially faster than ADMIXTURE for all tested K values. While ADMIXTURE’s runtime increases sharply as K grows—reaching over 14 hours at K=10—fastSTRUCTURE remains consistently below ~1.5 hours even at its peak. There show the Splitter Effect when increasing K.

- Both methods achieve high superpopulation assignment, with EUR, EAS, and SAS near 100% accuracy.

- On Soft-label performance (K=5): the overall Brier score is 0.11971 for ADMIXTURE and 0.12735 for fastSTRUCTURE, indicating similar performance.

- AMR shows complex multi-component ancestry in both methods, consistent with admixed population history.

## Citation
```bibtex
@article{alexander2009fast,
  title={Fast model-based estimation of ancestry in unrelated individuals},
  author={Alexander, David H and Novembre, John and Lange, Kenneth},
  journal={Genome research},
  volume={19},
  number={9},
  pages={1655--1664},
  year={2009},
  publisher={Cold Spring Harbor Lab}
}
```
```bibtex
@article{raj2014faststructure,
  title={fastSTRUCTURE: variational inference of population structure in large SNP data sets},
  author={Raj, Anil and Stephens, Matthew and Pritchard, Jonathan K},
  journal={Genetics},
  volume={197},
  number={2},
  pages={573--589},
  year={2014},
  publisher={Oxford University Press}
}
```
