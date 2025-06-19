#!/bin/bash
#SBATCH --job-name=PCAplink  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/PCA-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/PCA-%j.stderr         # Error log file
#SBATCH --cpus-per-task=64         # Number of CPU cores per task
#SBATCH --mem=128G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04

# plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/1kgp_b37_common_dedup_ldpruned"
# pca_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/pca_1kgp_b37_qc"

# plink2 \
#   --pfile "$plink_infile" \
#   --pca 10 \
#   --out "$pca_outfile"

# Save PCA with allele frequency weights
# plink2 --pfile "$plink_infile" \
#        --freq counts --ac-founders \
#        --pca allele-wts 10 \
#        --out "$pca_outfile"

# Step 2a-Prep 1: Create a variant info file defining A1 as ALT
# Assuming standard .pvar format: #CHROM POS ID REF ALT
# cat 1kgp_b37_common_dedup_ldpruned.pvar | awk 'NR>190 {print $3, $5}' > force_a1_alt.txt

# Step 2a-Prep 2: Create a new PGEN with A1 forced to ALT using the file
# echo "Creating new PGEN with A1 forced to ALT..."
# plink2 --pfile 1kgp_b37_common_dedup_ldpruned --a1-allele force_a1_alt.txt 1 2 --make-pgen --out 1kgp_b37_common_dedup_ldpruned_FORCED_A1_ALT

# # Step 2a-Main: Run PCA on the data with forced A1=ALT
# echo "Running PCA on data with forced A1=ALT..."
# plink2 --pfile plink/1kgp_b37_common_dedup_ldpruned_FORCED_A1_ALT --freq counts --ac-founders --pca allele-wts 10 --out pca_1kgp_loadings_FORCED_A1_ALT

# Filtering loadings file to keep only A1=ALT rows...
# awk 'BEGIN{OFS="\t"} NR==1 || $5 == $4 {print}' pca_1kgp_loadings_FORCED_A1_ALT.eigenvec.allele > pca_1kgp_loadings_FORCED_A1_ALT_filtered.eigenvec.allele
# awk 'BEGIN{OFS="\t"} NR==1 || $5 == $3 {print}' pca_1kgp_loadings_FORCED_A1_ALT.eigenvec.allele > pca_1kgp_loadings_FORCED_A1_ALT_filtered2.eigenvec.allele

# Use approx flag for UKBiobank
plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/plink/ukbb_imputation_gagneur_ids_common_dedup_ldpruned"
pca_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/pca_projected_ukbb_imputation_qc"
acount="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/pca_1kgp_loadings_FORCED_A1_ALT.acount"
eigenvec_allele="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1kgp_b37_common/pca_1kgp_loadings_FORCED_A1_REF_filtered.eigenvec.allele"

plink2 --pfile "$plink_infile" \
       --read-freq "$acount" \
       --score "$eigenvec_allele" 2 5 header-read no-mean-imputation variance-standardize \
       --score-col-nums 6-15 \
       --out "$pca_outfile"