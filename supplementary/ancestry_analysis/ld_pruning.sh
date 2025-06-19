#!/bin/bash
#SBATCH --job-name=LD_prune  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/LD_prune-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/LD_prune-%j.stderr         # Error log file
#SBATCH --cpus-per-task=64         # Number of CPU cores per task
#SBATCH --mem=128G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04


# plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/1kgp_b37_common_dedup"
# plink_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/1kgp_b37_common_dedup_ldpruned"
plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/ukbb_imputation_gagneur_ids_common_dedup"
plink_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/ukbb_imputation_gagneur_ids_common_dedup_ldpruned"

# Generate list of SNPs to keep after pruning 1KGP data - Use the same one for UKBB as well
# Run this first before the netx step
# plink2 \
#     --pfile "$plink_infile" \
#     --indep-pairwise 50 5 0.5 \
#     --out 1kgp_pruning_50_5_0.5

plink2 \
    --pfile "$plink_infile" \
    --extract 1kgp_pruning_50_5_0.5.prune.in \
    --make-pgen \
    --out "$plink_outfile"