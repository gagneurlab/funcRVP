# This folder contains shell scripts to run scripts on a SLURM cluster

Always run these scripts **from the root folder**.

- [`run_trainer.sh`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash/run_trainer.sh) runs the FuncRVP training. `sbatch ./bash/run_trainer.sh <INSERT TRAIT> <INSERT CONFIG_FILE>`

- [`submit_all_funcrvp.sh`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash/submit_all_funcrvp.sh) runs the FuncRVP training across multiple traits (defined in the config). Modify the paths directly in the shell script before running. `./bash/submit_all_funcrvp.sh`

- [`run_lm_phenopred.sh`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash/run_lm_phenopred.sh) runs the phenotype prediction using a linear model on burden test significant genes (user defined significance threshold). If the burden test results are not present in the same output folder (defined in the config file), it will run the burden test again. `sbatch ./bash/run_lm_phenopred.sh <INSERT TRAIT> <INSERT CONFIG_FILE> <INSERT P-VALUE>`

- [`run_ols.sh`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash/run_ols.sh) runs the burden test. `sbatch ./bash/run_ols.sh <INSERT TRAIT> <INSERT CONFIG_FILE>`

- [`submit_all_ols.sh`](https://github.com/gagneurlab/funcRVP/tree/orga-configs/bash/submit_all_ols.sh) runs the phenotype prediction using a linear model on burden test significant genes and the burden test on multiple traits (defined in the config). Modify the paths directly in the shell script before running. `./bash/submit_all_ols.sh`