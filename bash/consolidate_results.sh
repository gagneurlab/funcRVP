#!/bin/bash
#SBATCH --job-name=g2p  # Job name
#SBATCH --output=./logs/conres_out-%j.log       # Output log file
#SBATCH --error=./logs/conres_err-%j.log         # Error log file
#SBATCH --cpus-per-task=32         # Number of CPU cores per task
#SBATCH --mem=128G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0

python -m scripts.utils.utils consolidate_results --config-path=$1