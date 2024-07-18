# Bayesian genotype to phenotype model leveraging gene embeddings
---
Rare variant association testing is a promising approach to identify effector genes for common traits. However, ensuring sensitive and robust rare variant association testing is challenging due to the scarcity of high-impact rare-allele carriers. Here we introduce FuncRVP, a Bayesian rare variant association framework that addresses this issue by leveraging functional gene embeddings, i.e. multidimensional representations of gene function. FuncRVP models the accumulated effects of rare variants on traits as a weighted sum of rare-variant gene impairment scores. A prior, learnt from data, regularizes the weight of each gene depending on the location of the gene in a functional gene embedding.

This repo contains the model scripts to train and evaluate our model.
