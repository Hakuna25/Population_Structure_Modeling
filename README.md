# Comparing Global Ancestry Analysis Methods on 1000 Genomes Phase 3 Data

Xueqian Bai, Jialin Wu, Yutong Liang

## Introduction
Global ancestry analysis is a fundamental task in population genetics, aiming to infer the ancestral origins of individuals based on their genetic data. This analysis provides insights into human migration patterns, population structure, and evolutionary history. In this project, we focus on comparing two widely used model-based methods for global ancestry analysis: [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655) and [fastSTRUCTURE](https://doi.org/10.1534/genetics.114.164350). ADMIXTURE is a maximum-likelihood
admixture model that estimates individual ancestry proportions efficiently for large SNP datasets using numerical optimization. In contrast, fastSTRUCTURE adopts a Bayesian framework and applies variational inference to approximate posterior distributions of ancestry proportions. Both methods utilize genotype data to estimate the proportion of ancestry from different populations for each individual. By applying these methods to the 1000 Genomes Phase 3 dataset, we aim to evaluate their performance and compare the inferred ancestry proportions across different populations.

## Requirements
- Python 3.11
- Main environment: plink==1.90b7.7, admixture==1.3.1, faststructure==1.0
- Separate environment: bcftools == 1.23, samtools == 1.23

## Installation
```
# create a dedicated bcftools environment
conda env create -f env_bcftools.yml

# create main environment
conda env create -f env.yml
conda activate bio_tools
pip install structure_threader --user


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

- Open `analysis_setting_1.ipynb` for comparison on the main whole-genome LD-pruned.
- Open `analysis_setting_2.ipynb` for comparison on chromosome 20 with different MAF thresholds and LD pruning strategies.

## 🧬 Data and assets
[1000 Genomes Phase 3 Data](https://www.nature.com/articles/nature15393) is a comprehensive release of the 1000 Genomes Project dataset, providing whole-genome sequencing–based variant calls for 2,504 individuals from 26 populations across five continental groups. 

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
4. Open `analysis_setting_1.ipynb` or `analysis_setting_2.ipynb` for whole pipeline.


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
- `1000Genomes/`: input data, sample metadata, LD-pruned VCFs, and the reference FASTA files under `reference/`.
- `dump/`:
  - `common/`: shared outputs from `preprocess.sh`.
  - `admixture/`: the primary workspace for ADMIXTURE. Contains all intermediate PLINK files (`.bed/.bim/.fam`),  logs, fitted `P/Q` matrices, and benchmarking outputs.
  - `structure/`: the primary workspace for fastSTRUCTURE. Contains all intermediate PLINK files, logs, `meanQ/meanP` matrices, and benchmarking outputs.
  - `special_chr20_*`: chromosome 20 experiments under different MAF / LD settings.
- `plots/`: exported figures generated by the notebooks.
- `tools/`: local binary cache used for downloaded third-party tools such as ADMIXTURE.
- `.gitignore`: local ignore rules for generated files and environments.
- `README.md`: project overview, setup instructions, workflow, and summarized results.
- `admixture.sh`: runs ADMIXTURE on the shared `dump/common/common_ALL.pruned.*` dataset.
- `analysis_setting_1.ipynb`: main whole-autosome comparison between ADMIXTURE and fastSTRUCTURE.
- `analysis_setting_2.ipynb`: chromosome 20 comparison under different MAF / LD configurations.
- `benchmark.sh`: measuring runtime and peak-memory from pipeline runs.
- `env.yml`: main Conda environment for PLINK, ADMIXTURE, plotting, and evaluation.
- `env_bcftools.yml`: dedicated Conda environment for `bcftools` and `samtools`.
- `pipeline.conf`: central configuration for binary paths, threads, and preprocessing directories.
- `preprocess.sh`: shared preprocessing pipeline that normalizes VCFs, builds PLINK files, LD pruning, and merges chromosomes.
- `structure.sh`: runs fastSTRUCTURE on the shared `dump/common/common_ALL.pruned.*` dataset.
- `utils.py`: helper utilities used by the notebooks and plotting workflow.


## 📊 Results and Conclusion
The full analysis is documented in `analysis_setting_1.ipynb` and `analysis_setting_2.ipynb`. Below is a README-sized summary of the notebook outputs.

### Setting 1: All Autosomes

On the whole-autosome dataset, both methods recover the five 1000 Genomes super-populations clearly, with fastSTRUCTURE producing visually cleaner blocks and ADMIXTURE preserving more low-level mixed ancestry signal. K-selection points to a stable range around `K=5-7` rather than one single sharp optimum. For ADMIXTURE, the CV curve has a clear elbow at `K=5`: the error drops from `0.50675` at `K=2` to `0.47627` at `K=5`, and then improves by only `0.00105` from `K=5` to `K=10`. For fastSTRUCTURE, the per-`K` marginal likelihood is best at `K=6` (`-0.86454`) and stays nearly flat for `K=5..10`; separately, `chooseK` reports `5` components needed to explain the structure in the data and `7` as the maximum model complexity. Taken together, these results support using `K=5` as the main super-population baseline while inspecting `K=4-7` for finer structure.

#### Whole-autosome K-value Exploration:

![Whole-autosome CV and marginal-likelihood curves](plots/all_autosomes_cv_ml_vs_k.png)

#### Whole-autosome Population Structure:
![Whole-autosome ancestry panels](plots/dump_admixture_combined_membership_panels.png)

#### Whole-autosome Runtime and Memory:

| Method | Total fit time across `K=2..10` | Mean fit time per K | Peak RSS |
| --- | ---: | ---: | ---: |
| ADMIXTURE | 267,246 s | 29,694 s | 6.19 GB |
| fastSTRUCTURE | 158,028 s | 17,559 s | 3.95 GB |

- fastSTRUCTURE is clearly more memory efficient on the whole-autosome benchmark: its peak memory stays between about `2.57 GB` and `3.95 GB`, versus `4.76 GB` to `6.19 GB` for ADMIXTURE, or roughly `54%` to `64%` of ADMIXTURE's memory footprint across the tested `K` values.
- Runtime is more nuanced than the memory story. fastSTRUCTURE is faster overall in aggregate, reducing total fit time by about `41%`, but it is not uniformly faster at every `K`: in this run it is slower than ADMIXTURE at `K=3-5`, with the largest slowdown at `K=4` (`22,612 s` vs `10,560 s`).
- At larger `K`, the trend reverses strongly in fastSTRUCTURE's favor. The biggest speedup appears at `K=10`, where fastSTRUCTURE needs `11,975 s` compared with `58,395 s` for ADMIXTURE, and ADMIXTURE also shows pronounced runtime jumps at `K=6` and `K=8-10`.

![Whole-autosome runtime and memory vs K](plots/all_autosomes_runtime_memory_vs_k.png)

#### Bonus Evaluation (Not in report)
| Evaluation on whole autosomes | ADMIXTURE | fastSTRUCTURE |
| --- | ---: | ---: |
| Hard-assignment accuracy at `K=5` | 0.90935 | 0.90974 |
| Macro-F1 at `K=5` | 0.86697 | 0.86894 |
| Multi-label Brier score at `K=5` | 0.12257 | 0.12085 |

- EAS and EUR are essentially perfect under the hard-assignment evaluation, AFR is above 99% for both methods, and AMR is the hardest group because many samples split between AMR and EUR components.
- The per-super-population Brier scores show the same pattern: AMR is by far the most difficult class (`0.76408` for ADMIXTURE, `0.80059` for fastSTRUCTURE), while EAS, EUR, and SAS are modeled much more cleanly.
- As `K` increases, both methods show the same high-level continental backbone, but ADMIXTURE shows a more visible "splitter effect" and finer within-region fragmentation.


### Setting 2: Chromosome 20 (Chr20) Quality-Control (QC) Experiments

The setting compares four QC configurations and shows that LD pruning and MAF thresholds strongly affect the computational cost while preserving the main continental structure.

####  Chr20 Population Structure:

`MAF=0.01 + LD prune`

![Chr20 MAF=0.01 + LD prune](plots/special_chr20_maf001_pruned_admixture_combined_membership_panels.png)

`MAF=0.01 without LD prune`

![Chr20 MAF=0.01 without LD prune](plots/special_chr20_maf001_no_ld_admixture_combined_membership_panels.png)

`MAF=0.005 + LD prune`

![Chr20 MAF=0.005 + LD prune](plots/special_chr20_maf0005_pruned_admixture_combined_membership_panels.png)

`MAF=0.05 + LD prune`

![Chr20 MAF=0.05 + LD prune](plots/special_chr20_maf005_pruned_admixture_combined_membership_panels.png)

####  Chr20 K-value exploration Across QC Configurations:

![Chr20 CV and marginal-likelihood curves](plots/cv_ml_vs_k.png)

| Chr20 setting | SNPs | ADMIXTURE mean fit time / peak RSS | fastSTRUCTURE mean fit time / peak RSS | fastSTRUCTURE `chooseK` |
| --- | ---: | ---: | ---: | ---: |
| `MAF=0.005 + LD prune` | 67,555 | 1517.7 s / 480.7 MB | 1439.9 s / 351.5 MB | 6 |
| `MAF=0.01 + LD prune` | 32,319 | 801.8 s / 233.9 MB | 784.4 s / 424.4 MB | 6 |
| `MAF=0.05 + LD prune` | 6,021 | 137.3 s / 50.3 MB | 416.7 s / 425.4 MB | 5 |
| `MAF=0.01 without LD prune` | 267,644 | 6997.8 s / 1885.8 MB | 6079.3 s / 1292.6 MB | 5 |

- Removing LD pruning is the most expensive choice for both methods: compared with `MAF=0.01 + LD prune`, the SNP count jumps from `32,319` to `267,644`, and both runtime and memory increase sharply.
- Raising the MAF threshold to `0.05` makes the problem much smaller and especially benefits ADMIXTURE, but the population structure visualization panels become visibly coarser because only `6,021` SNPs remain.
- Lowering the MAF threshold to `0.005` keeps more rare variation (`67,555` SNPs) and yields panels close to the whole-autosome analysis, with slightly more fine-scale heterogeneity, especially in ADMIXTURE.
- Across all four Chr20 settings, ADMIXTURE CV tends to keep decreasing slightly toward larger `K`, while fastSTRUCTURE peaks earlier and more stably at `K=5` or `K=6`, so the practical resolution remains around the same range as the whole-genome experiment.
- In the fastSTRUCTURE memory plot, `MAF=0.01 + LD prune` and `MAF=0.05 + LD prune` are visually the closest pair: from `K=3` onward the two curves are nearly on top of each other, with differences as small as `0.02 MB` and an average gap of about `10.48 MB` across `K=2..10`.

#### Chr20 Runtime and Memory Across QC Configurations:

![Chr20 runtime and memory vs K](plots/runtime_memory_vs_k.png)



## ✈️ Future Directions
- Explore additional methods such as PCA-based approaches or non-parametric clustering for a broader comparison.
- Extend the analysis to other datasets with different population structures or sequencing technologies.


## References
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
```bibtex
@article{10002015global,
  title={A global reference for human genetic variation},
  author={1000 Genomes Project Consortium and others},
  journal={Nature},
  volume={526},
  number={7571},
  pages={68},
  year={2015}
}
```
