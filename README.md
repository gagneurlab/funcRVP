# FuncRVP: Bayesian genotype to phenotype model leveraging gene embeddings
-------------------------
Rare variant association testing is a promising approach to identify effector genes for common traits. However, ensuring sensitive and robust rare variant association testing is challenging due to the scarcity of high-impact rare-allele carriers. Here we introduce FuncRVP, a Bayesian rare variant association framework that addresses this issue by leveraging functional gene embeddings, i.e. multidimensional representations of gene function. FuncRVP models the accumulated effects of rare variants on traits as a weighted sum of rare-variant gene impairment scores. A prior, learnt from data, regularizes the weight of each gene depending on the location of the gene in a functional gene embedding. Want to know more about FuncRVP? Read our preprint: [https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2](https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2)

To try out the FuncRVP model on simulated example data run the [python notebook in the example directory](https://github.com/gagneurlab/funcRVP/tree/orga-configs/example).

This repo contains the scripts to train our model and benchmark the results.  

## Data
--------------------------
To run FuncRVP you need the following data processed:
- **Gene impairment scores matrix (samples x genes)**. This can be either simple count of pLoF variants in a gene or something more informative like a DeepRVAT score.
- **Phenotype (samples x 1)**. The phenotypes for the corresponding samples.
- **Covariates (samples x covariates)**. The covariates, i.e., age, sex, genetic PCs for the corresponding samples.
- **Gene embeddings (genes x embedding_dimension)**. The gene embeddings for the genes you want to use for the phenotype prediction. Only the gene for which an embedding is provided will be used for phenotype prediction. We test multiple embeddings in the manuscript, all of which can be found in the [`scripts/utils/gene_embeddings`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/utils/gene_embeddings) folder. 

## Requirements
--------------------------
- Linux, Python (tested with v3.9, v3.11, and v3.12)
- NVIDIA GPU (tested on A40 and L40S) for training on large datasets

Training FuncRVP for 50 epochs on ~300,000 samples x ~18,000 genes from the UK biobank takes about 10 minutes on an NVIDIA L40S for a single trait.


# Repo structure
--------------------------

## Scripts

The [`scripts`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts) directory contains the codes for the running the model, hyperparameter optimization and consolidating results.

## Bash

The [`bash`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash) directory contains the shell scripts used for the running the python scripts on a SLURM cluster.

## Figures

The [`figures`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/figures) directory contains R scripts used to create figures in our manuscript. The scripts are labeled according to the figure they generate in the manuscript.

## Supplementary Figures

The [`supplementary`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/supplementary) directory contains R scripts used to create figures in our manuscript. The scripts are labeled according to the figure they generate in the manuscript.
