import pandas as pd
import argparse
import pickle
from tqdm import tqdm
import os
import ast

# Create ArgumentParser object
parser = argparse.ArgumentParser(description='model run info')

# Add arguments
parser.add_argument('--output_genes', type=str, help='Output path for genes')
parser.add_argument('--output_pheno', type=str, help='Output path for pheno')
parser.add_argument('--study_name', type=str, help='Study name')
parser.add_argument('--embedding', type=str, default='omics_pops', help='Embedding')
parser.add_argument('--genotype', type=str, default='deepRVAT', help='Genotype matrix type')
parser.add_argument('--test_split_size', type=float, default=0.25, help='Test split size')
parser.add_argument('--logs_dir', type=str, help='Logs directory')
parser.add_argument('--all_beta_results', type=str, help='All beta results')
parser.add_argument('--all_pred_results', type=str, help='All pred results')
parser.add_argument('--traits', type=str, help='List of traits')

# Parse command-line arguments
args = parser.parse_args()

output_genes = args.output_genes
output_pheno = args.output_pheno
study_name = args.study_name
embedding = args.embedding
genotype = args.genotype
test_size = args.test_split_size
logs_dir = args.logs_dir
all_beta_results = args.all_beta_results
all_pred_results = args.all_pred_results
traits = args.traits

TRAITS = ast.literal_eval(traits)

with open('/s/project/geno2pheno/data/replication_sets/genebass_backman_gene_dict.pkl', 'rb') as f:
    replication_dict = pickle.load(f)
    
with open('/s/project/geno2pheno/data/trait_phenocode_dict.pkl', 'rb') as f:
    genebass_phenocode_dict = pickle.load(f)

genebass_phenocode_dict_flipped = {v: k for k, v in genebass_phenocode_dict.items()}

gb_rep = pd.read_parquet('/s/project/geno2pheno/data/replication_sets/genebass_backman_41traits.pq')

genes_dt_list = []
pheno_dt_list = []
skip_list = []
for trait in tqdm(TRAITS):
    try:
        phenocode = genebass_phenocode_dict[trait]
        # Concatenate gene effect files across all traits
        mean_betas_df = pd.read_parquet(f"{os.path.join(all_beta_results.split(trait+'_mean_betas.pq')[0], trait+'_mean_betas.pq')}")
        mean_betas_df["phenocode"] = str(phenocode)
        mean_betas_df["replicated"] = mean_betas_df.index.isin(replication_dict[phenocode])
        genes_dt_list.append(mean_betas_df)
        
        # Concatenate phenotype prediction files across all traits
        bayes_pred_df = pd.read_parquet(f"{os.path.join(all_pred_results.split(trait+'_bayes_pred.pq')[0], trait+'_bayes_pred.pq')}").reset_index()
        bayes_pred_df["phenocode"] = str(phenocode)
        bayes_pred_df["genotype"] = genotype
        pheno_dt_list.append(bayes_pred_df)
    except:
        skip_list.append(trait)
        continue

genes_dt = pd.concat(genes_dt_list)
pheno_dt = pd.concat(pheno_dt_list)

genes_dt.to_parquet(output_genes, index=False)
pheno_dt.to_parquet(output_pheno, index=False)
print(skip_list)