This folder contains scripts to train the FuncRVP model.

# To train FuncRVP
--------------------

- The `funcrvp_training.py` script is used for training FuncRVP on a single trait. `submit_funcrvp_training.sh` is used to submit the script as a SLURM cluster job. To submit training jobs for multiple traits, `submit_funcrvp_training_jobs.sh` is used.

- The `funcrvp_hyperopt_repeats.py` script is used for training FuncRVP on a single trait and hyperparameter optimization. The same combination of hyperparameters can be submitted to evaluate the robustness of the hyperpatameter combination. `submit_funcrvp_hyperopt.sh` is used to submit the script as a SLURM cluster job. To submit training jobs for multiple traits, `submit_funcrvp_hyperopt_jobs.sh` is used.

- The directory [model_prior_mean](https://github.com/gagneurlab/funcRVP/tree/main/scripts/model_prior_mean) contains scripts to run a variant of the FuncRVP model that models the mean of the prior instead of the variance.
