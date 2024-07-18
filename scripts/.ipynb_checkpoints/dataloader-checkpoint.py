import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import pairwise_distances

EMBEDDINGS = {
    "string": "/s/project/gene_embedding/embedding/final/STRING_d128.tsv",
    "string_exp": "/s/project/gene_embedding/embedding/final/STRING_EXP_d128.tsv",
    "omics": "/s/project/gene_embedding/embedding/final/dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding.tsv",
    "pops": "/s/project/gene_embedding/embedding/final/pops_mat_d256.tsv",
    "pops_exp": "/s/project/gene_embedding/embedding/final/pops_mat_exp_d256.tsv",
    "omics_pops": "/s/project/gene_embedding/embedding/combination/pops_mat_pca256_omics.tsv",
    "omics_pops_exp": "/s/project/gene_embedding/embedding/combination/pops_mat_exp256_omics.tsv",
    "omics_pops_rescaled": "/s/project/gene_embedding/embedding/combination/omics_pops256_rescaled.tsv",
    "frogs": "/s/project/geno2pheno/data/embeddings/frogs_archs4_d256.tsv"
}

def load_data(trait, embedding_type=None, reduce_embedding_dim=512, shuffle_embeddings=False, add_loeuf=False, add_alpha_mis=False, add_gwas_hits=False, only_genebass_genes=False, extend_covariates=False, normalize=False, version=None, genotype="deepRVAT", ukb_version=None):
    
    if version.startswith("filtered"):
        if genotype=="pLoF":
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_ensembl_loftee_plof_binarized_pivot.parquet"
        elif genotype=="pLoF_raw":
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_ensembl_loftee_plof_pivot.parquet"
        elif genotype=="deepRVAT":
            # genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_avg_20_repeats_pivot.parquet"
            # genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_mean_30_01_25.parquet"          # NEW DEEPRVAT SCORES !!!!!!!!!!!!!!!!-----------------!!!!!!!!!!!!!!!!
            if ukb_version=='500':
                genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_500k_DeepRVAT_mean_6_0.parquet"                 # NEWER DEEPRVAT SCORES !!!!!!!!!!!!!!!!-----------------!!!!!!!!!!!!!!!!
            else:
                genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_mean_6_0.parquet"                 # NEWER DEEPRVAT SCORES !!!!!!!!!!!!!!!!-----------------!!!!!!!!!!!!!!!!
        else:
            print("defaulting to pLoF genotype in dataloader")
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_ensembl_loftee_plof_binarized_pivot.parquet"
        
        if ukb_version=='500':
            dataset_df = pd.read_parquet(f"/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_300k_{version}.parquet")
        else:
            dataset_df = pd.read_parquet(f"/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_{version}.parquet")
        
    else:
        if genotype=="pLoF":
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_loftee_plof_pivot.parquet"
        elif genotype=="deepRVAT":
            # genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_avg_20_repeats_pivot.parquet"
            # genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_mean_30_01_25.parquet"          # NEW DEEPRVAT SCORES !!!!!!!!!!!!!!!!-----------------!!!!!!!!!!!!!!!!
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_DeepRVAT_mean_6_0.parquet"               # NEWER DEEPRVAT SCORES !!!!!!!!!!!!!!!!-----------------!!!!!!!!!!!!!!!!
        else:
            print("defaulting to pLoF genotype in dataloader")
            genotype_path = "/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_loftee_plof_pivot.parquet"

        if version:
            dataset_df = pd.read_parquet(f"/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_{version}.parquet")
        else:
            dataset_df = pd.read_parquet("/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals.parquet")
        
    GT_df = pd.read_parquet(genotype_path)
     
    dataset_df = pd.merge(dataset_df, GT_df[[]].reset_index())
        
    dataset_df = dataset_df[list(dataset_df.columns[:23])+[trait, trait+"_PRS", trait+"_split", trait+"_common_resid"]]
    dataset_df = dataset_df.query(f"~{trait}_split.isna()")
    
    train_dataset_df = dataset_df.groupby(f"{trait}_split").get_group("train").reset_index(drop=True)
    val_dataset_df = dataset_df.groupby(f"{trait}_split").get_group("val").reset_index(drop=True)
    test_dataset_df = dataset_df.groupby(f"{trait}_split").get_group("test").reset_index(drop=True)
    
    y_train = train_dataset_df[trait+"_common_resid"].values
    y_val = val_dataset_df[trait+"_common_resid"].values
    y_test = test_dataset_df[trait+"_common_resid"].values
    
    trait_measurement_train = train_dataset_df[trait].values
    trait_measurement_val = val_dataset_df[trait].values
    trait_measurement_test = test_dataset_df[trait].values
    
    covariates_train = train_dataset_df[list(dataset_df.columns[1:23])+[trait+"_PRS"]].copy()
    covariates_val = val_dataset_df[list(dataset_df.columns[1:23])+[trait+"_PRS"]].copy()
    covariates_test = test_dataset_df[list(dataset_df.columns[1:23])+[trait+"_PRS"]].copy()
    
    if extend_covariates:
        covariates_train["age:sex"] = covariates_train["age"] * covariates_train["sex"]
        covariates_val["age:sex"] = covariates_val["age"] * covariates_val["sex"]
        covariates_test["age:sex"] = covariates_test["age"] * covariates_test["sex"]
        
        covariates_train["age^2"] = covariates_train["age"] * covariates_train["age"]
        covariates_val["age^2"] = covariates_val["age"] * covariates_val["age"]
        covariates_test["age^2"] = covariates_test["age"] * covariates_test["age"]
        
        covariates_train["age^2:sex"] = covariates_train["age^2"] * covariates_train["sex"]
        covariates_val["age^2:sex"] = covariates_val["age^2"] * covariates_val["sex"]
        covariates_test["age^2:sex"] = covariates_test["age^2"] * covariates_test["sex"]
    
    covariates_train = covariates_train.values
    covariates_val = covariates_val.values
    covariates_test = covariates_test.values
    
    id_train = train_dataset_df["individual"]
    id_val = val_dataset_df["individual"]
    id_test = test_dataset_df["individual"]
    
    if normalize:
        cov_df = pd.DataFrame(np.concatenate([covariates_train, covariates_val, covariates_test]), index=np.concatenate((id_train, id_val, id_test)))
        cov_trainval_df = pd.DataFrame(np.concatenate([covariates_train, covariates_val]), index=np.concatenate((id_train, id_val)))
        cov_df = (cov_df-cov_trainval_df.mean())/cov_trainval_df.std()
        (covariates_train, covariates_val, covariates_test) = cov_df.loc[id_train].values, cov_df.loc[id_val].values, cov_df.loc[id_test].values
        # measurements_df = pd.DataFrame(np.concatenate([trait_measurement_train, trait_measurement_val, trait_measurement_test]), index=np.concatenate((id_train, id_val, id_test)))
        # measurements_trainval_df = pd.DataFrame(np.concatenate([trait_measurement_train, trait_measurement_val]), index=np.concatenate((id_train, id_val)))
        # measurements_df = (measurements_df-measurements_trainval_df.mean())/measurements_trainval_df.std()
        # (trait_measurement_train, trait_measurement_val, trait_measurement_test) = measurements_df[0].loc[id_train].values, measurements_df[0].loc[id_val].values, measurements_df[0].loc[id_test].values
    
    gt_train = pd.merge(GT_df, train_dataset_df.set_index("individual")[[]], left_index=True, right_index=True).sort_index()
    
    gene_list = sorted(list(gt_train.columns))
    
    if only_genebass_genes:
        test_genes_df = pd.read_csv("/s/project/uk_biobank/processed/derived_datasets/genebass_500k_test_genes.csv").query("annotation=='pLoF' and ttype=='burden'")[["gene_id"]]
        gene_list = sorted(set(gene_list).intersection(set(test_genes_df["gene_id"])))
    
    if embedding_type:
        gene_embeddings_df = pd.read_csv(EMBEDDINGS[embedding_type], sep="\t").rename(columns={"gene_id": "gene"}).sort_values("gene")
        gene_list = sorted(set(gene_embeddings_df.sort_values("gene")["gene"]).intersection(set(gene_list)))
        
        # if add_loeuf or add_alpha_mis or add_gwas_hits:
        #     # loeuf_df = pd.read_csv("/s/raw/gnomad/storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression='gzip', sep="\t")[["gene_id", "oe_lof_upper"]].dropna().rename(columns={"gene_id": "gene"})
        #     loeuf_df = pd.read_csv("/s/project/geno2pheno/data/constrain_scores.tsv", sep="\t")[["gene", "oe_lof_upper"]].dropna()
        #     gene_list = sorted(set(gene_list).intersection(set(loeuf_df["gene"])))
        
        gene_embeddings_df = gene_embeddings_df.set_index("gene").loc[gene_list].reset_index()
        embedding_components = reduce_embedding_dim
        if gene_embeddings_df.shape[1]-1 > embedding_components:
            pca = PCA(n_components=embedding_components)
            embeddings_pcs = pca.fit_transform(gene_embeddings_df.set_index("gene").values)
            print(f"After embedding PCA: {round(sum(pca.explained_variance_ratio_), 4)} explained variance remaining")
            gene_embeddings_df = pd.DataFrame(embeddings_pcs, index=gene_embeddings_df.set_index("gene").index, columns=[f"emb_col_{c}" for c in range(embedding_components)]).reset_index(drop=False)
            
        if add_loeuf:
            loeuf_df = pd.read_csv("/s/project/geno2pheno/data/constrain_scores.tsv", sep="\t")[["gene", "oe_lof_upper"]].dropna()
            gene_embeddings_df = pd.merge(gene_embeddings_df, loeuf_df)
            gene_list = sorted(set(gene_list).intersection(set(loeuf_df["gene"])))
            
        if add_alpha_mis:
            alpha_df = pd.read_csv("/s/project/geno2pheno/data/constrain_scores.tsv", sep="\t")[["gene", "mean_am_pathogenicity"]].dropna()
            gene_embeddings_df = pd.merge(gene_embeddings_df, alpha_df)
            gene_list = sorted(set(gene_list).intersection(set(alpha_df["gene"])))
            
        if add_gwas_hits:
            gwas_scores_df = pd.read_parquet('/s/project/geno2pheno/data/panukb_hits/max_pval_hq.pq', columns=['gene',trait])
            gene_embeddings_df = pd.merge(gene_embeddings_df, gwas_scores_df)
            gene_list = sorted(set(gene_list).intersection(set(gwas_scores_df["gene"])))
        
        gene_embeddings_df = gene_embeddings_df.set_index("gene")
        gene_embeddings_df = gene_embeddings_df[gene_embeddings_df.index.isin(gene_list)]
        gene_embeddings_df = gene_embeddings_df[~gene_embeddings_df.index.duplicated(keep='first')]
        emb = gene_embeddings_df.values
        
    else:
        emb = None
        
    gt_train = gt_train[gene_list]
    gt_train = gt_train.values

    gt_val = pd.merge(GT_df, val_dataset_df.set_index("individual")[[]], left_index=True, right_index=True).sort_index()
    gt_val = gt_val[gene_list]
    gt_val = gt_val.values

    gt_test = pd.merge(GT_df, test_dataset_df.set_index("individual")[[]], left_index=True, right_index=True).sort_index()
    gt_test = gt_test[gene_list]
    gt_test = gt_test.values
    
    if genotype=="deepRVAT":
        gt_mean = gt_train.mean(0)
        gt_std = gt_train.std(0)
        gt_train = (gt_train-gt_mean)#/gt_std
        gt_val = (gt_val-gt_mean)#/gt_std
        gt_test = (gt_test-gt_mean)#/gt_std
    
    return (gt_train, gt_val, gt_test), (y_train, y_val, y_test), emb, gene_list, (id_train, id_val, id_test), (trait_measurement_train, trait_measurement_val, trait_measurement_test), (covariates_train, covariates_val, covariates_test)

def load_common_residuals_statistics(trait, version="v3"):
    dataset_df = pd.read_parquet(f"/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_{version}.parquet").query(f"~{trait}_split.isna()")
    gt = pd.merge(pd.read_parquet("/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_loftee_plof_pivot.parquet"), dataset_df.set_index("individual")[[]], left_index=True, right_index=True).sort_index()
    
    mean_pheno_pLoF = (gt.values * dataset_df[trait+"_common_resid"].values[:, np.newaxis]).sum(0)
    n_carriers = gt.values.sum(0)
    mean_pheno_pLoF = pd.Series(np.divide(mean_pheno_pLoF, n_carriers, out=np.zeros_like(mean_pheno_pLoF), where=n_carriers!=0), index = gt.columns)
    n_carriers = pd.Series(n_carriers, index = gt.columns)
    return pd.DataFrame({"n_carriers": n_carriers, "mean_pheno_pLoF": mean_pheno_pLoF}).reset_index(names="gene_id")

def load_common_residuals_statistics_by_split(trait, version="filteredv3"):
    dataset_df = pd.read_parquet(f"/s/project/uk_biobank/processed/derived_datasets/ukbb_common_traits_and_covariates_with_split_and_residuals_{version}.parquet").query(f"~{trait}_split.isna()")
    gt = pd.merge(pd.read_parquet("/s/project/uk_biobank/processed/derived_datasets/ukbb_wes_200k_ensembl_loftee_plof_binarized_pivot.parquet"), dataset_df.set_index("individual")[[]], left_index=True, right_index=True).sort_index()

    gt_val = pd.merge(gt, dataset_df.query(f"{trait}_split=='val'").set_index("individual")[[]], left_index=True, right_index=True).sort_index()
    gt_test = pd.merge(gt, dataset_df.query(f"{trait}_split=='test'").set_index("individual")[[]], left_index=True, right_index=True).sort_index()

    statistics_list = []

    mean_pheno_pLoF = (gt.values * dataset_df[trait+"_common_resid"].values[:, np.newaxis]).sum(0)
    n_carriers = gt.values.sum(0)

    mean_pheno_pLoF = pd.Series(np.divide(mean_pheno_pLoF, n_carriers, out=np.zeros_like(mean_pheno_pLoF), where=n_carriers!=0), index = gt.columns, name="mean_pheno")
    n_carriers = pd.Series(n_carriers, index = gt.columns, name="n_carriers")

    statistics_list.append(mean_pheno_pLoF)
    statistics_list.append(n_carriers)

    for split in ["train", "val", "test"]:
        gt_split= pd.merge(gt, dataset_df.query(f"{trait}_split=='{split}'").set_index("individual")[[]], left_index=True, right_index=True).sort_index()

        split_mean_pheno_pLoF = (gt_split.values * dataset_df.query(f"{trait}_split=='{split}'")[trait+"_common_resid"].values[:, np.newaxis]).sum(0)
        split_n_carriers = gt_split.values.sum(0)

        split_mean_pheno_pLoF = pd.Series(np.divide(split_mean_pheno_pLoF, split_n_carriers, out=np.zeros_like(split_mean_pheno_pLoF), where=split_n_carriers!=0), index = gt_split.columns, name=f"mean_pheno_{split}")
        split_n_carriers = pd.Series(split_n_carriers, index = gt_split.columns, name=f"n_carriers_{split}")

        statistics_list.append(split_mean_pheno_pLoF)
        statistics_list.append(split_n_carriers)

    mean_pheno_pLoF_df = pd.DataFrame(statistics_list).T
    return mean_pheno_pLoF_df

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
