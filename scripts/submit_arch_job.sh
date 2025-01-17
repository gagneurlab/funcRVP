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

# sbatch -p urgent submit_arch_job.sh High_light_scatter_reticulocyte_percentage esm2 deepRVAT 0.25 v1NEWsplit
# sbatch -p urgent submit_arch_job.sh LDL_direct pops deepRVAT 0.25 NEWsplit
# sbatch -p urgent submit_arch_job.sh BodyMassIndex omics_pops deepRVAT 0.25 v1NEWsplit
# sbatch -p urgent submit_arch_job.sh systolic_blood_pressure omics_pops deepRVAT 0.25 v1NEWsplit_noEmb --no_emb
# sbatch -p urgent submit_arch_job.sh LDL_direct omics_pops deepRVAT 0.25 NEWsplit_randEmb --random_emb 0
