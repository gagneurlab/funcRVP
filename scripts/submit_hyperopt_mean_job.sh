#!/bin/bash
#SBATCH --job-name=g2p_hyperopt  # Job name
#SBATCH --output=./logs/meanho-%j-%a_output.log       # Output log file
#SBATCH --error=./logs/meanho-%j-%a_error.log         # Error log file
#SBATCH --cpus-per-task=16         # Number of CPU cores per task
#SBATCH --mem=64G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:a40:1
#SBATCH --exclude=ouga12,ouga14

python hyperopt_repeats_mean.py $1 $2 $3 $4

