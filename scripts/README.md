# This folder contains the scritps used to train FuncRVP and the baseline models

## Training FuncRVP
--------------------

- The [`training/train_model.py`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/training/train_model.py) script is used for training FuncRVP on a single trait. 

To train FuncRVP run the following command from the root dir (or run the relevant shell script in [`bash`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash)): `python -m scripts.training.train_model trainer --trait=<INSERT TRAIT> --config-path=<INSERT PATH>`

- The [`models/g2p_bayes_cov_skipcon.py`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/models/g2p_bayes_cov_skipcon.py) is the file containing the FuncRVP model class. This class is imported in the training scripts.

## Baseline models
----------------------

- The [`baselines/train_ols.py`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/) contains the functions to perform a burden test for each gene, and also subsequently perform phenotype prediction using the significant (using a user defined threshold) genes. 

To perform the burden test run the following command from the root dir (or run the relevant shell script in [`bash`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash)): `python -m scripts.baselines.train_ols_multi_anc association-test --trait=<INSERT TRAIT> --config-path=<INSERT PATH>`

To perform phenotype prediction run the following command from the root dir (or run the relevant shell script in [`bash`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash)): `python -m scripts.baselines.train_ols lm-phenopred --trait=<INSERT TRAIT> --config-path=<INSERT PATH> --pval-threshold=<INSERT P-VALUE>`

- [`model_prior_mean`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/model_prior_mean) contains scripts to run a version of the FuncRVP model that models the mean of the prior rather than the variance.

## Util scripts
--------------------

- [`utils/dataloader_old.py`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/utils/dataloader_old.py) is the file containing the dataloader class. This class converts precomputed gene burdens, covariates, and trait measurements into the format required by FuncRVP. This is used by all FuncRVP and burden test training scripts. The script requires a gene burden score table in the format: (samples x genes), and a covariate table in the format: (samples x covariates).

- [`utils/utils.py`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/scripts/utils/utils.py) contains other important utility functions used by other scripts.
