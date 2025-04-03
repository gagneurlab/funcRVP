#!/bin/bash
#SBATCH --job-name=ols_rvat  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/ols-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/ols-%j.stderr         # Error log file
#SBATCH --cpus-per-task=32         # Number of CPU cores per task
#SBATCH --mem=200G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

python -m scripts.baselines.train_ols lm-phenopred --trait=$1 --config-path=$2

# python -m scripts.baselines.train_ols association-test --trait=LDL_direct --config-path=run_config_test.yaml
# python -m scripts.baselines.train_ols lm-phenopred --trait=LDL_direct --config-path=run_config_test.yaml