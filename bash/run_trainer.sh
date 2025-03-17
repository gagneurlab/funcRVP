#!/bin/bash
#SBATCH --job-name=g2p  # Job name
#SBATCH --output=./logs/trainer_out-%j.log       # Output log file
#SBATCH --error=./logs/trainer_err-%j.log         # Error log file
#SBATCH --cpus-per-task=16         # Number of CPU cores per task
#SBATCH --mem=128G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:a40:1

python -m scripts.training.train_model trainer --trait=$1 --config-path=$2

# python -m scripts.training.train_model trainer --trait=LDL_direct --config-path=run_config.yaml