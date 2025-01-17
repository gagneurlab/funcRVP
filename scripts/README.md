This folder contains scripts to train the FuncRVP model.

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
