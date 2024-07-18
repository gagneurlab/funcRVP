import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
import optuna
import dataloader
import os

version = "filteredv3"
genotype = "deepRVAT"
ukb_version = '300'

sign_genes = False
only_test = False
shuffled_phenotype = False


traits = list(pd.read_csv("phenotype_list_41.txt", header=None)[0].values)

genebass_genes_df = pd.read_csv("/s/project/uk_biobank/processed/derived_datasets/genebass_500k_test_genes.csv").query("annotation=='pLoF' and ttype=='burden'")[["gene_id"]]

genebass_phenocode_dict = pd.read_csv("/s/project/rep/processed/trait_associations_v5/ukbb_wes_200k/PRS_score_mapping.csv")[["phenotype", "genebass_phenocode"]].set_index("genebass_phenocode")["phenotype"].to_dict()
genebass_phenocode_dict = {v: k for k, v in genebass_phenocode_dict.items()}
genebass_phenocode_dict["LDL_direct"] = "30780"
genebass_phenocode_dict["Albumin"] = "30600"

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

genebass_phenocode_dict_flipped = {v: k for k, v in genebass_phenocode_dict.items()}


traits = list(pd.read_csv("phenotype_list_41.txt", header=None)[0].values)

genes_df_list = []
pred_df_list = []

for i, trait in enumerate(traits):
    phenocode = genebass_phenocode_dict[trait]
    
    assoc_df = pd.read_parquet(f"/s/project/uk_biobank/processed/g2p/modelling/ols_{genotype}/{version}_ols_{genotype}_cov{'_shuffledpheno' if shuffled_phenotype else ''}{'_test' if only_test else ''}{f'_{ukb_version}' if ukb_version else ''}_full_{trait}.pq")
    if sum(assoc_df["neglog_pval"].isna())>0:
        print("NaN pvalues")
    
    assoc_df['phenocode'] = str(phenocode)
    genes_df_list.append(assoc_df)
    assoc_df["significant"] = assoc_df["pval"] < 0.05/len(assoc_df)

    (gt_train, gt_val, gt_test), (resid_train, resid_val, resid_test), _, gene_list, (id_train, id_val, id_test), (trait_measurement_train, trait_measurement_val, trait_measurement_test), (covariates_train, covariates_val, covariates_test) = dataloader.load_data(trait, embedding_type = None, extend_covariates=True, normalize=True, version=version, genotype=genotype, ukb_version=ukb_version)

    gt_trainval = pd.DataFrame(np.concatenate([gt_train, gt_val]), columns=gene_list)[assoc_df.query("significant").index].values
    gt_test = pd.DataFrame(gt_test, columns=gene_list)[assoc_df.query("significant").index].values

    covariates_trainval = np.concatenate([covariates_train, covariates_val])
    
    for sign_genes in [True, False]:
        print(trait, f"{i+1}/{len(traits)}")

        output_dir = f"/s/project/uk_biobank/processed/g2p/modelling/LM_significant_genes/{genotype}"
        outfile_pred = f"{output_dir}/{trait}_LM{'_sign_genes' if sign_genes else ''}_{version}{'_shuffledpheno' if shuffled_phenotype else ''}{'_test' if only_test else ''}{f'_{ukb_version}' if ukb_version else ''}.pq"
        
        # if os.path.isfile(outfile_pred):
        #     continue
            
        if sign_genes:
            x_trainval = np.concatenate([covariates_trainval, gt_trainval], axis=1)
            x_test = np.concatenate([covariates_test, gt_test], axis=1)
        else:
            x_trainval = covariates_trainval
            x_test = covariates_test
        
        lm = LinearRegression().fit(x_trainval, np.concatenate([trait_measurement_train, trait_measurement_val]))

        # print("test score", lm.score(x_test, trait_measurement_test))
        print("test score trait_measurement", r2_score(trait_measurement_test , lm.predict(x_test)))

        pred_df_list.append(pd.DataFrame({f"{trait}_measurement": trait_measurement_test, 
                                          "common_residual": resid_test, 
                                          "lm_pred": lm.predict(x_test), 
                                          "trait": trait, 
                                          "phenocode": str(phenocode), 
                                          "model": f"lm{'_sign_genes' if sign_genes else ''}_{genotype}"}, 
                                         index=id_test))#.to_parquet(outfile_pred)
    
    
    # lm_sign_genes_df["phenocode"] = str(phenocode)
    # lm_sign_genes_df["model"] = f"lm{'_sign_genes' if sign_genes else ''}_{genotype}"
    
lm_sign_genes_df = pd.concat(pred_df_list)    
lm_sign_genes_df["genotype"] = genotype
    

ols_df = pd.concat(genes_df_list).reset_index(drop=True)
# ols_df["phenocode"] = str(phenocode)
ols_df["var_beta"] = ols_df["std_err"]**2
ols_df["neglog_pval"] = -np.log(ols_df["pval"])
ols_df["significant"] = ols_df["neglog_pval"] > -np.log(0.05/len(ols_df))
ols_df["replicated"] = ols_df.index.isin(replication_dict[phenocode])
ols_df["genotype"] = genotype
ols_df["model"] = ols_df["model"] + "_cov"

    
    
    
# for study_version in study_versions:
#     study_version = study_version.rstrip(',')

ols_df.reset_index(drop=True).to_parquet(f"/s/project/geno2pheno/predictions/bayesian/ols_{version}_{genotype}{f'_{ukb_version}' if ukb_version else ''}_genes_extended.pq")
lm_sign_genes_df.reset_index(drop=True).to_parquet(f"/s/project/geno2pheno/predictions/bayesian/lm_{version}_{genotype}{f'_{ukb_version}' if ukb_version else ''}_predictions_extended.pq")
