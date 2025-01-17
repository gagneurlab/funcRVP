This folder contains scripts to train the FuncRVP model.

- dataloader_clean.py: Loads genotype data, phenotype data, and embeddings.
- g2p_bayes_cov_skipcon.py: Defines the model class
- train_model.py: Runs model without hyperparameter optimization on a user defined architecture
- hyperopt_repeats.py: Hyperparameter optimization
- submit_hyperopt_jobs.sh & submit_fixed_arch_jobs.sh: Bash scripts to submit training jobs on a SLURM cluster

- The directory [model_prior_mean](https://github.com/gagneurlab/funcRVP/tree/main/scripts/model_prior_mean) contains scripts to run a variant of the FuncRVP model that models the mean of the prior instead of the variance.
