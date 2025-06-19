#!/bin/bash
#SBATCH --job-name=var_inter  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/var_inter-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/var_inter-%j.stderr         # Error log file
#SBATCH --cpus-per-task=100         # Number of CPU cores per task
#SBATCH --mem=400G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

# plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/1kgp_b37/all_phase3_unrelated_uniqueID"
plink_infile="/s/project/uk_biobank/processed/Imputation/ukbb_imputation"
snplist="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/ukbb_1kgp_b37_common.snplist"
plink_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/ukbb_imputation_gagneur_ids_common"

plink2 \
    --pfile "$plink_infile" \
    --extract "$snplist" \
    --make-pgen \
    --out "$plink_outfile"

# plink2 --pfile 1kgp_b37_common --rm-dup force-first --max-alleles 2 --make-pgen --out 1kgp_b37_common_dedup
# plink2 --pfile ukbb_imputation_gagneur_ids_common --rm-dup force-first --max-alleles 2 --make-pgen --out ukbb_imputation_gagneur_ids_common_dedup