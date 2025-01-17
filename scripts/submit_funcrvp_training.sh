#!/bin/bash
#SBATCH --job-name=g2pTrait  # Job name
#SBATCH --output=./logs/g2pTrait-%j_output.log       # Output log file
#SBATCH --error=./logs/g2pTrait-%j_error.log         # Error log file
#SBATCH --cpus-per-task=16         # Number of CPU cores per task
#SBATCH --mem=128G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:a40:1
#SBATCH --exclude=ouga12

python funcrvp_training.py --trait=$1 --embedding=$2 --genotype=$3 --test_split_size=$4 --version=$5
# python funcrvp_training.py --trait=$1 --embedding=$2 --genotype=$3 --test_split_size=$4 --version=$5 $6 --random_seed=$7

