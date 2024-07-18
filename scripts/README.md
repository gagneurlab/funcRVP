This folder contains scripts to train the funcRVP model.

- dataloader.py: Loads genotype data, phenotype data, and embeddings.
- g2p_bayes_cov_skipcon.py: Model class
- example_fit_model.ipynb: Notebook to run the model
- hyperopt_repeats.py: Hyperparameter optimization
- run_architecture.py: Runs model without hyperparameter optimization on a user defined architecture
- submit_hyperopt_jobs.sh & submit_fixed_arch_jobs.sh: Bash scripts to submit training jobs to the cluster
