#!/bin/bash
#SBATCH --job-name=g2p_ho  # Job name
#SBATCH --output=./logs/ho-%j_output.log       # Output log file
#SBATCH --error=./logs/ho-%j_error.log         # Error log file
#SBATCH --cpus-per-task=16         # Number of CPU cores per task
#SBATCH --mem=128G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:a40:1
# SBATCH --exclude=ouga12

python funcrvp_hyperopt_repeats.py $1 $2 $3 $4 $5
