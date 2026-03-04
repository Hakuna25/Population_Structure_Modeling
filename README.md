# CSE284 Final Project

[Xueqian Bai](https://hakuna25.github.io/), [Jialin Wu](https://jlwu.cn/), [Yutong Liang](https://www.lyt0112.com/)

## Introduction
Global ancestry analysis is a fundamental task in population genetics, aiming to infer the ancestral origins of individuals based on their genetic data. This analysis provides insights into human migration patterns, population structure, and evolutionary history. In this project, we focus on comparing two widely used model-based methods for global ancestry analysis: [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655) and [fastSTRUCTURE](https://doi.org/10.1534/genetics.114.164350). Both methods utilize genotype data to estimate the proportion of ancestry from different populations for each individual. By applying these methods to the 1000 Genomes Phase 3 dataset, we aim to evaluate their performance and compare the inferred ancestry proportions across different populations.

## Requirements
- Python 3.11
- plink==1.90b7.7, admixture==1.3.0, faststructure==1.0

## Installation
```
conda env create -f env.yml
pip install structure_threader --user
```

## Data and assets
[1000 Genomes Phase 3 Data](https://www.nature.com/articles/nature15393) is a comprehensive release of the 1000 Genomes Project dataset, providing whole-genome sequencing–based variant calls for 2,504 individuals from 26 populations across five continental groups. In this project, we focus on autosomal biallelic SNPs and use an LD-pruned version of the autosomal data.

For reproducibility, please place all input 1000 Genomes files in `1000Genomes/`.
For detailed download instructions and expected filenames, see [1000Genomes/README.md](https://github.com/Hakuna25/Population_Structure_Modeling/blob/main/1000Genomes/README.md).

## Reproducible Workflow
1. Create environment and install dependencies.
2. Place 1000 Genomes input files in `1000Genomes/`.
3. Edit `pipeline.conf` for config. 
4. Open `analysis.ipynb` and run from top to bottom.
 We’ve already provided the preprocessed, merged 1000 Genomes Phase 3 chromosome dataset in `dump/`, so you can run the tools directly. If you’re interested, you can also reproduce the preprocessing steps using the snippets included in the code.

### Additional Notes
- `preprocess.sh` contains shared preprocessing steps used by both methods:
  - sample extraction
  - merge-list creation
  - merged PLINK generation
  - optional per-chromosome preprocessing loop (`--run-chr-process`)

- `bash test.sh` provides a quick test for the environment setup and method implementation. Test outputs are written to `dump/test/`.

## Runtime & Memory Benchmarking
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

## Repository structure
- `1000Genomes/`: contains the original input data, specifically the igsr_samples.tsv and the prunded VCF files.
- `dump/`:
    - `admixture/`: the primary workspace for Admixture. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
    - `structure/`: the primary workspace for fastStructure. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
- Root scripts: `preprocess.sh` are the shared preprocessing commands; `admixture.sh` and `structure.sh` are method implementations; `analysis.ipynb` is the main entry point for method calling and results analysis.

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
