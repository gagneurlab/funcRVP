# Bayesian genotype to phenotype model leveraging gene embeddings
---
Rare variant association testing is a promising approach to identify effector genes for common traits. However, ensuring sensitive and robust rare variant association testing is challenging due to the scarcity of high-impact rare-allele carriers. Here we introduce FuncRVP, a Bayesian rare variant association framework that addresses this issue by leveraging functional gene embeddings, i.e. multidimensional representations of gene function. FuncRVP models the accumulated effects of rare variants on traits as a weighted sum of rare-variant gene impairment scores. A prior, learnt from data, regularizes the weight of each gene depending on the location of the gene in a functional gene embedding. Want to know more about or work? Read our preprint: [https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2](https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2)

To try out the FuncRVP model on simulated example data run the [python notebook in the example directory](https://github.com/gagneurlab/funcRVP/tree/main/example).

This repo contains the scripts to train our model and benchmark the results. The [scripts](https://github.com/gagneurlab/funcRVP/tree/main/scripts) directory contains the codes for the running the model, hyperparameter optimization and consolidating results. The [plotting_codes](https://github.com/gagneurlab/funcRVP/tree/main/plotting_codes) directory contains codes to recreate figures from our manuscript.

# Requirements
---
- Linux, Python (tested with v3.9)
- NVIDIA GPU (tested on A40) for training on biobank scale data

