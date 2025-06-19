#!/bin/bash
#SBATCH --job-name=QCplink  # Job name
#SBATCH --output=/s/project/geno2pheno/funcrvp/logs/imp_qc-%j.stdout       # Output log file
#SBATCH --error=/s/project/geno2pheno/funcrvp/logs/imp_qc-%j.stderr         # Error log file
#SBATCH --cpus-per-task=32         # Number of CPU cores per task
#SBATCH --mem=16G                  # Memory allocation per node (adjust as needed)
#SBATCH --gres=gpu:0
#SBATCH --exclude=ouga03,ouga04


# ukbb_imputation_file="/s/project/uk_biobank/processed/Imputation/ukbb_imputation"
plink_infile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/1000genomes_b37/all_phase3_unrelated_uniqueID"
qc_outfile="/s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_1000genomes_b37_common/all_phase3_unrelated_uniqueID"

plink2 \
  --pfile "$plink_infile" \
  --maf 0.1 \
  --geno 0.05 \
  --hwe 1e-5 \
  --write-snplist allow-dups \
  --write-samples --no-id-header \
  --allow-extra-chr \
  --out "$qc_outfile"

  # plink2 --pfile all_phase3 --remove deg2_phase3.king.cutoff.out.id --make-pgen --out all_phase3_unrelated
  # plink2 --pfile all_phase3_unrelated --set-all-var-ids chr@:#:\$r\>\$a --make-pgen --out all_phase3_unrelated_uniqueID --new-id-max-allele-len 662

  # plink2 --pfile all_hg38 --remove deg2_hg38.king.cutoff.out.id --make-pgen --out all_hg38_unrelated --allow-extra-chr
  # plink2 --pfile all_hg38_unrelated --set-all-var-ids @:#:\$r\>\$a --new-id-max-allele-len 487 --allow-extra-chr --make-pgen --out uniqueID/all_hg38_unrelated_uniqueID
  # plink2 --pfile all_hg38_unrelated_uniqueID --rm-dup force-first --make-pgen --out all_hg38_unrelated_uniqueID_rm_dup --allow-extra-chr

  # HapMap3 
  # plink2 --pedmap hapmap3_r1_b36_fwd_consensus.qc.poly.recode --set-all-var-ids @:#:\$r\>\$a --make-pgen --out uniqueID/hapmap3_hg37_uniqueID
  # plink2 --pfile hapmap3_hg37_uniqueID --maf 0.1 --geno 0.05 --hwe 1e-5 --write-snplist allow-dups --write-samples --no-id-header --out /s/project/geno2pheno/funcrvp/paper_revisions/supplementary_results/ancestry_analysis/backmanQC_hapmap3_common/hapmap3_hg37_uniqueID