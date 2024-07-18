#!/bin/bash
#SBATCH --job-name=g2pArch  # Job name
#SBATCH --output=./logs/g2pArch-%j_output.log       # Output log file
#SBATCH --error=./logs/g2pArch-%j_error.log         # Error log file
#SBATCH --cpus-per-task=16         # Number of CPU cores per task
#SBATCH --mem=64G                 # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:a40:1
#SBATCH --exclude=ouga14

python run_architecture.py --trait=$1 --embedding=$2 --genotype=$3 --version=$4

# sbatch submit_arch_job_shuffledpheno.sh High_light_scatter_reticulocyte_count omics_pops deepRVAT v5e5arch_5e4const 1234