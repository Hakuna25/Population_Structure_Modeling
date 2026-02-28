# CSE284 Final Project

[Xueqian Bai](https://hakuna25.github.io/), Jialin Wu, Yutong Liang

## Introduction
Briefly describe the project here...

## Requirements
- Python 3.11
- plink==1.90b7.7, admixture==1.3.0, faststructure==1.0
## Installation
```
# create conda environment
conda env create -f env.yml
pip install structure_threader --user
```
## Data and assets
Here, describe the dataset...
We have already included in repo ./1000Genomes... For more info, u can refer to ...(1000genomes link)
## Running the tools
The bash scripts in this repo cover data preprocessing and population structure analysis. We’ve already provided the preprocessed, merged 1000 Genomes Phase 3 chromosome dataset in ./dump, so you can run the tools directly. If you’re interested, you can also reproduce the preprocessing steps using the snippets included in the code.
```
# from the repo root:
bash admixture.sh   # ADMIXTURE
bash structure.sh   # fastStructure
```
Run ```analysis.ipynb``` to view the visualizations and analysis results.
## Repository structure
- `1000Genomes/`: contains the original input data, specifically the igsr_samples.tsv and the prunded VCF files.
- `dump/`:
    - `admixture/`: the primary workspace for Admixture. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
    - `structure/`: the primary workspace for fastStructure. Contains all intermediate PLINK files (.bed, .bim, .fam) for chromosomes and running logs.
- `script/`: contains main bash scripts, and Python scripts used for downstream visualization.

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