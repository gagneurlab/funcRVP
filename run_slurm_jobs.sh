#!/bin/bash
# Simplified Snakemake SLURM Submission Script

set -e  # Exit on error

# ----- Configuration -----
snakemake="snakemake"  # Path to Snakemake (assuming it's in PATH)
snakefile="${SNAKEFILE:-Snakefile}"  # Snakefile path (default to 'Snakefile')
number_of_snakemake_cores="${N_CORES:-256}"  # Snakemake cores (default 256)
number_of_snakemake_jobs="${N_JOBS:-300}" # Number of snakemake jobs (default 128)
amount_of_snakemake_memory="${MEM_MB:-128000}" # Amount of snakemake memory
number_of_snakemake_gpus="${N_GPUS:-1}" # Number of gpus

# Project directory and log folder
project_folder="$(dirname "$snakefile")"
project_folder="$(realpath "$project_folder")"
logs="$project_folder/logs"
job_name="$(basename "$project_folder")-$(date +"%Y-%m-%d_%H-%M-%S")"

# ----- SLURM Configuration -----
sbatch_args="--requeue" #Add requeue to sbatch parameters
#Kerberos and auks configuration (optional)
if [ "${AUKS_ENABLED:-false}" = true ]; then
    kinit="/usr/bin/kinit"
    $kinit -r 7d
    auks -a
    sbatch_args="$sbatch_args --auks=done"
fi
# Create log folder
mkdir -p "$logs"

# ----- Snakemake Command -----
$snakemake --keep-going \
    --default-resources ntasks=1 mem_mb=32000 gpu=0 \
    --cluster "sbatch $sbatch_args \
        --ntasks {resources.ntasks} \
        --cpus-per-task {threads} \
        --mem {resources.mem_mb}M \
        --job-name=$job_name-{rule} \
        --gres=gpu:{resources.gpu} \
        --output=$logs/%x-%j.out \
        --error=$logs/%x-%j.err \
        --parsable \
        " \
    --cluster-status="${0%/*}/slurm-status.py" \
    --cores "$number_of_snakemake_cores" \
    -j "$number_of_snakemake_jobs" \
    --resources mem_mb="$amount_of_snakemake_memory" gpu="$number_of_snakemake_gpus" \
    --snakefile "$snakefile" "$@"