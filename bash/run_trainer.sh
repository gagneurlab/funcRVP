#!/bin/bash
#SBATCH --job-name=g2p  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/trainer-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/trainer-%j.stderr         # Error log file
#SBATCH --cpus-per-task=42         # Number of CPU cores per task
#SBATCH --mem=166G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:l40s:1
#SBATCH --exclude=ouga05,ouga06,ouga08

python -m scripts.training.train_model trainer --trait=$1 --config-path=$2

# python -m scripts.training.train_model trainer --trait=LDL_direct --config-path=run_config_test.yaml

# Example job
# sbatch -p urgent ./bash/run_trainer.sh LDL_direct run_config_local.yaml