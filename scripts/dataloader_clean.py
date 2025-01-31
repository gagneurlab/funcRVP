import os
import sys
import pandas as pd
import numpy as np
from sklearn.preprocessing import quantile_transform

EMBEDDINGS = {
    "string": "/s/project/geno2pheno/data/embeddings/STRING_d256.tsv",
    "string_exp": "/s/project/geno2pheno/data/embeddings/STRING_EXP_d256.tsv",
    "omics": "/s/project/geno2pheno/data/embeddings/omics_d256.tsv",
    "pops": "/s/project/geno2pheno/data/embeddings/pops_mat_d256.tsv",
    "pops_exp": "/s/project/geno2pheno/data/embeddings/pops_mat_exp_d256.tsv",
    "omics_pops": "/s/project/geno2pheno/data/embeddings/pops_mat_pca256_omics.tsv",
    "omics_pops_exp": "/s/project/geno2pheno/data/embeddings/pops_mat_exp256_omics.tsv",
    "esm2": "/s/project/geno2pheno/data/embeddings/ESM2_PCA_d512.tsv",
    "gene2vec": "/s/project/geno2pheno/data/embeddings/gene2vec_d200.tsv",
    "enformer": "/s/project/geno2pheno/data/embeddings/enformer_d3072.tsv",
    "enformer_small": "/s/project/geno2pheno/data/embeddings/enformer_d512.tsv",
    "frogs": "/s/project/geno2pheno/data/embeddings/frogs_archs4_d256.tsv"
}

genotype_path_dict = {
    'plof': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_vep_plof.parquet",
    'anti_plof': "/s/project/uk_biobank/processed/derived_datasets/temp_files/ukbb_wes_500k_plof_mask_301124.parquet",
    'deeprvat': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_final_090924_medianshifted.parquet",
    'deeprvat_pub': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_final_090924.parquet",
    'deepRVAT_allmodels': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_notclean_all_30_models_151024.parquet",
    'deeprvat_plof': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_plof_031224_medianshifted.parquet",
    'deeprvat_pub_plof': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_plof_301124.parquet",
    'deeprvat_pub_no_plof': "/s/project/uk_biobank/processed/derived_datasets/temp_files/ukbb_wes_500k_DeepRVAT_no_plof_301124.parquet",
    'alphamissense': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_alphamissense_300924.parquet",
    'am_plof': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_sum_alphamissense_plof_151024.parquet",
    'primateai': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_primateai_011024.parquet",
    'cadd': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_cadd_011024.parquet",
    'max_plof_am_deeprvat': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_max_plof_am_deeprvat_061224.parquet",
    'max_plof_am_deeprvat_medianshifted': "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_max_plof_am_deeprvat_061224_medianshifted.parquet",
}

def trait_INT_zscore(trait):
    pheno_dt = pd.read_parquet(f"/s/project/uk_biobank/processed/decoded_phenotypes/{trait}/data.parquet")
    trait_mean = pheno_dt[f'{trait}_raw'].mean()
    trait_sd = pheno_dt[f'{trait}_raw'].std()
    pheno_dt['zscore'] = (pheno_dt[f'{trait}_raw'] - trait_mean)/trait_sd
    pheno_dt['raw'] = pheno_dt[f'{trait}_raw']
    pheno_dt['INT'] = pheno_dt[f'{trait}']
    pheno_dt['new_INT'] = quantile_transform(pheno_dt[['zscore']],output_distribution='normal', random_state=0, copy=True)
    pheno_dt['trait'] = trait
    
    pheno_dt = pheno_dt.rename(columns={'eid':'individual'})
    pheno_dt['individual'] = pheno_dt['individual'].astype('str')
    # pheno_dt = pheno_dt.set_index('individual')
    return pheno_dt[['individual', 'trait', 'raw', 'zscore', 'INT', 'new_INT']]


def load_data(trait, embedding_type=None, use_prs=True, normalize_covariates=True, version='filteredv3', genotype="deeprvat", test_split_size=0.4, val_split_size=0.1, split_seed=0, small_gene_list=None):

    if genotype in genotype_path_dict.keys():
        genotype_path = genotype_path_dict[genotype]
    else:
        print("defaulting to pLOF genotype in dataloader")
        genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_vep_plof.parquet"

    # Read GIS matrix
    if small_gene_list:     # Want to run it on a subset of genes
        GT_df = pd.read_parquet(genotype_path, columns=small_gene_list).reset_index()
    else:
        GT_df = pd.read_parquet(genotype_path).reset_index()
        
    gene_list = sorted(list(GT_df.columns))

    # Read embedding vectors
    if embedding_type:
        gene_embeddings_df = pd.read_csv(EMBEDDINGS[embedding_type], sep="\t").rename(columns={"gene_id": "gene"}).sort_values("gene")
        gene_list = sorted(set(gene_embeddings_df.sort_values("gene")["gene"]).intersection(set(gene_list)))
        gene_embeddings_df = gene_embeddings_df.set_index("gene").loc[gene_list].reset_index()
        
        gene_embeddings_df = gene_embeddings_df.set_index("gene")
        gene_embeddings_df = gene_embeddings_df[gene_embeddings_df.index.isin(gene_list)]
        gene_embeddings_df = gene_embeddings_df[~gene_embeddings_df.index.duplicated(keep='first')]
        emb = gene_embeddings_df.values

    else:
        emb = None

    # Filter GIS to include only thos genes whose embedding is availabe
    GT_df = GT_df[['individual'] + gene_list]

    # Read covariates and PRS
    covs_dt = pd.read_parquet('/s/project/uk_biobank/processed/derived_datasets/ukbb_500k_covariates_filteredv3.pq')
    prs_df = pd.read_parquet('/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_500k_filteredv3.parquet', columns=['individual', f'{trait}_PRS', f'{trait}_common_resid']).dropna()

    # Read common phenotype
    raw_pheno = trait_INT_zscore(trait)
    raw_pheno = raw_pheno.merge(prs_df[['individual', f'{trait}_common_resid']], on='individual', how='inner')

    # Get individuals for who all data is available
    if use_prs:
        covs_dt = covs_dt.merge(prs_df[['individual', f'{trait}_PRS']], on='individual', how='inner')
        inds = raw_pheno.merge(covs_dt[['individual']], how='inner').merge(GT_df[['individual']], how='inner')[['individual']]
    else:
        inds = raw_pheno.merge(covs_dt[['individual']], how='inner').merge(GT_df[['individual']], how='inner')[['individual']]
    
    # Individuals to be used only in the train split
    train_inds = pd.read_parquet('/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_filteredv3.parquet', columns=['individual'])

    # Assign 'train' to rows where the 'group' column is in train_groups
    inds['split'] = np.nan
    inds.loc[inds['individual'].isin(train_inds.individual), 'split'] = 'train'
    
    # Compute train/val/test proportions based on some criteria
    test_proportion = test_split_size*inds.shape[0]/(inds.shape[0] - inds[inds.split=='train'].shape[0])
    val_proportion = val_split_size*inds.shape[0]/(inds.shape[0] - inds[inds.split=='train'].shape[0])                     # Keep 30k samples in val
    remaining_train_proportion = 1 - (test_proportion + val_proportion)

    # Set a seed for split reproducibility
    np.random.seed(split_seed)
    
    # Randomly assign splits to the remaining rows
    remaining_rows = inds['split'].isna()
    inds.loc[remaining_rows, 'split'] = np.random.choice(['train', 'val', 'test'], size=remaining_rows.sum(), p=[remaining_train_proportion, val_proportion, test_proportion])

    # Return train/val/test indices
    id_train = inds[inds.split=='train']["individual"]
    id_val = inds[inds.split=='val']["individual"]
    id_test = inds[inds.split=='test']["individual"]

    # Convert all dataframes to numpy arrays for downstream models
    G_train = inds[inds.split=='train'][['individual']].merge(GT_df).set_index('individual').values
    G_val = inds[inds.split=='val'][['individual']].merge(GT_df).set_index('individual').values
    G_test = inds[inds.split=='test'][['individual']].merge(GT_df).set_index('individual').values
    
    Cov_train = inds[inds.split=='train'][['individual']].merge(covs_dt).set_index('individual').values
    Cov_val = inds[inds.split=='val'][['individual']].merge(covs_dt).set_index('individual').values
    Cov_test = inds[inds.split=='test'][['individual']].merge(covs_dt).set_index('individual').values

    # Fix normalization
    if normalize_covariates:
        cov_trainval_means = np.concatenate([Cov_train, Cov_val]).mean(axis=0)
        cov_trainval_sd = np.concatenate([Cov_train, Cov_val]).std(axis=0)
        
        Cov_train = (Cov_train-cov_trainval_means)/cov_trainval_sd
        Cov_val = (Cov_val-cov_trainval_means)/cov_trainval_sd
        Cov_test = (Cov_test-cov_trainval_means)/cov_trainval_sd

    # if use_residuals:
    residuals_train = inds[inds.split=='train'][['individual']].merge(raw_pheno).set_index('individual')[[f'{trait}_common_resid']].values.squeeze()
    residuals_val = inds[inds.split=='val'][['individual']].merge(raw_pheno).set_index('individual')[[f'{trait}_common_resid']].values.squeeze()
    residuals_test = inds[inds.split=='test'][['individual']].merge(raw_pheno).set_index('individual')[[f'{trait}_common_resid']].values.squeeze()
    # else:
    trait_measurement_train = inds[inds.split=='train'][['individual']].merge(raw_pheno).set_index('individual')[['new_INT']].values.squeeze()
    trait_measurement_val = inds[inds.split=='val'][['individual']].merge(raw_pheno).set_index('individual')[['new_INT']].values.squeeze()
    trait_measurement_test = inds[inds.split=='test'][['individual']].merge(raw_pheno).set_index('individual')[['new_INT']].values.squeeze()

    return (G_train, G_val, G_test), (residuals_train, residuals_val, residuals_test), emb, gene_list, (id_train, id_val, id_test), (trait_measurement_train, trait_measurement_val, trait_measurement_test), (Cov_train, Cov_val, Cov_test)
    



def load_phenocode_dict():
    genebass_phenocode_dict = pd.read_csv("/s/project/rep/processed/trait_associations_v5/ukbb_wes_200k/PRS_score_mapping.csv")[["phenotype", "genebass_phenocode"]].set_index("genebass_phenocode")["phenotype"].to_dict()
    genebass_phenocode_dict = {v: k for k, v in genebass_phenocode_dict.items()}
    genebass_phenocode_dict["LDL_direct"] = "30780"
    genebass_phenocode_dict["Albumin"] = "30600"
    return genebass_phenocode_dict

def load_replication_set(traits, key="phenocode"):
    genebass_phenocode_dict = pd.read_csv("/s/project/rep/processed/trait_associations_v5/ukbb_wes_200k/PRS_score_mapping.csv")[["phenotype", "genebass_phenocode"]].set_index("genebass_phenocode")["phenotype"].to_dict()
    genebass_phenocode_dict = {v: k for k, v in genebass_phenocode_dict.items()}
    genebass_phenocode_dict["LDL_direct"] = "30780"
    genebass_phenocode_dict["Albumin"] = "30600"

    #genebass_results = pd.read_parquet("/s/project/rep/processed/genebass/500k/results.parquet", columns=["annotation", "gene_id", "phenocode", "BETA_Burden", "SE_Burden", "Pvalue_Burden"]).query(f"phenocode=='{genebass_phenocode_dict[trait]}'").rename(columns={"Pvalue_Burden": "PVAL_Burden"})

    genebass_df = pd.read_parquet('/s/project/uk_biobank/processed/g2p/modelling/replication_set/genebass_pvals_500k_selected_traits.parquet', engine='pyarrow').query("(significant_burden and keep_gene_burden) or (significant_skato and keep_gene_skato)")#.query("annotation=='pLoF'")
    genebass_df = genebass_df[genebass_df["phenocode"].isin([genebass_phenocode_dict[trait] for trait in traits])]
    #genebass_df = pd.merge(genebass_df, genebass_results, how="left", on=["annotation", "gene_id", "phenocode"])
    genebass_df["association"] = genebass_df["gene_id"] + "_" + genebass_df["phenocode"]

    backman_df = pd.read_excel('/s/project/uk_biobank/processed/g2p/modelling/replication_set/41586_2021_4103_MOESM5_ESM.xlsx', sheet_name='SD2', engine='openpyxl')[['Marker type', 'Gene', 'UKB detailed trait name', 'Marker']]
    backman_df = backman_df[backman_df['Marker type'] == 'Burden']
    backman_df["phenocode"] = backman_df["UKB detailed trait name"].str.split("_", expand=True)[0]
    backman_pheno = [genebass_phenocode_dict[trait] for trait in traits]
    markers = ['M1.1', 'M3.1']
    backman_df = backman_df.query('phenocode == @backman_pheno and Marker in @markers')
    protein_coding_genes_df = pd.read_parquet("/s/project/uk_biobank/processed/g2p/modelling/replication_set/protein_coding_genes.parquet")

    backman_df = pd.merge(backman_df, protein_coding_genes_df, left_on="Gene", right_on="gene_name")
    backman_df["gene_id"] = backman_df["gene"].str.split(".", expand=True)[0]
    backman_df["association"] = backman_df["gene_id"] + "_" + backman_df["phenocode"]

    replication_association_set = set(genebass_df["association"].to_list() + backman_df["association"].to_list())
    replication_dict = pd.concat([genebass_df[["phenocode", "gene_id"]], backman_df[["phenocode", "gene_id"]]]).groupby("phenocode")["gene_id"].apply(set).to_dict()
    
    if key=="trait":
        genebass_phenocode_dict_flipped = {v: k for k, v in genebass_phenocode_dict.items()}
        replication_dict = {genebass_phenocode_dict_flipped[phenocode]:replication_dict[phenocode] for phenocode in replication_dict}

    return replication_dict

def get_pairwise_distances(embedding_type, gene_list=None, metric="euclidean"):
    gene_embeddings_df = pd.read_csv(EMBEDDINGS[embedding_type], sep="\t").rename(columns={"gene_id": "gene"}).sort_values("gene")
    if gene_list:
        pairwise_dist = pairwise_distances(gene_embeddings_df.set_index("gene"), gene_embeddings_df.set_index("gene").loc[gene_list], metric=metric)
        return pd.DataFrame(data=pairwise_dist, index=gene_embeddings_df["gene"].to_list(), columns=gene_list)
    else:
        pairwise_dist = pairwise_distances(gene_embeddings_df.set_index("gene"), metric=metric)
        return pd.DataFrame(data=pairwise_dist, index=gene_embeddings_df["gene"].to_list(), columns=gene_embeddings_df["gene"].to_list())
