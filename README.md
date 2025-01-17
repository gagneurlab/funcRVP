# Bayesian genotype to phenotype model leveraging gene embeddings
-------------------------
Rare variant association testing is a promising approach to identify effector genes for common traits. However, ensuring sensitive and robust rare variant association testing is challenging due to the scarcity of high-impact rare-allele carriers. Here we introduce FuncRVP, a Bayesian rare variant association framework that addresses this issue by leveraging functional gene embeddings, i.e. multidimensional representations of gene function. FuncRVP models the accumulated effects of rare variants on traits as a weighted sum of rare-variant gene impairment scores. A prior, learnt from data, regularizes the weight of each gene depending on the location of the gene in a functional gene embedding. Want to know more about FuncRVP? Read our preprint: [https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2](https://www.biorxiv.org/content/10.1101/2024.07.22.604535v2)

To try out the FuncRVP model on simulated example data run the [python notebook in the example directory](https://github.com/gagneurlab/funcRVP/tree/main/example).

This repo contains the scripts to train our model and benchmark the results. The [`scripts`](https://github.com/gagneurlab/funcRVP/tree/main/scripts) directory contains the codes for the running the model, hyperparameter optimization and consolidating results. The [`plotting_codes`](https://github.com/gagneurlab/funcRVP/tree/main/plotting_codes) directory contains R scripts used to create figures in our manuscript.

## Requirements
-------------------------- 
- Linux, Python (tested with v3.9)
- NVIDIA GPU (tested on A40) for training on biobank scale data

Training FuncRVP for 50 epochs on ~250,000 samples x ~18,000 genes from the UK biobank takes about 30 minutes on a NVIDIA A40 GPU for a single trait. 


# Scripts

## Training FuncRVP
--------------------

- The `funcrvp_training.py` script is used for training FuncRVP on a single trait. `submit_funcrvp_training.sh` is used to submit the script as a SLURM cluster job. To submit training jobs for multiple traits, `submit_funcrvp_training_jobs.sh` is used.

- The `funcrvp_hyperopt_repeats.py` script is used for training FuncRVP on a single trait and hyperparameter optimization. The same combination of hyperparameters can be submitted to evaluate the robustness of the hyperpatameter combination. `submit_funcrvp_hyperopt.sh` is used to submit the script as a SLURM cluster job. To submit training jobs for multiple traits, `submit_funcrvp_hyperopt_jobs.sh` is used.


## Supplementary scripts
--------------------

- `dataloader_clean.py` is the file containing the dataloader class. This class converts precomputed gene burdens, covariates, and trait measurements into the format required by FuncRVP. This is used by all FuncRVP and burden test training scripts. The script requires a gene burden score table in the format: (samples x genes), and a covariate table in the format: (samples x covariates).
  
- `g2b_bayes_cov_skipcon.py` is the file containing the FuncRVP model class. It is used by all FuncRVP training scripts. The script requires a gene burden score table in the format: (samples x genes), and a covariate table in the format: (samples x covariates).

- `train_ols.py` performs a burden test across all genes for a given trait.
  
- The notebook `read_results.ipynb` is used after running FuncRVP model training. It compiles phenotype prediction and gene effect estimates of FuncRVP (or the burden test) across all traits into two tables.

- The directory [model_prior_mean](https://github.com/gagneurlab/funcRVP/tree/main/scripts/model_prior_mean) contains scripts to run a variant of the FuncRVP model that models the mean of the prior instead of the variance.


# Plotting codes

This folder contains the R scripts used to create the figures in the FuncRVP manuscript. The scripts are labeled according to the figure they generate in the manuscript.
