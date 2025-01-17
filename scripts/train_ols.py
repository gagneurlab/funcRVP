import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso, LassoCV
import optuna
# import dataloader
import dataloader_clean
import os
import sys
import statsmodels.api as sm
from tqdm import tqdm

def main():
    # Retrieve command-line arguments
    arguments = sys.argv[1:]  # Exclude the script name

    # Process and use the command-line arguments
    trait = arguments[0]
    genotype = arguments[1]
    # ukb_version = arguments[2]
    test_size = float(arguments[2])
    only_test = arguments[3] in ['t', 'T', 'True']
    # only_test = False
    
    version = "filteredv3"    
    use_prs = True
    normalize_covariates = True
    
    covariates = True
    shuffled_phenotype = False
    
    output_dir = f"/s/project/uk_biobank/processed/g2p/modelling/ols_{genotype}/"
    
    outfile = f"{output_dir}{version}_at_{genotype}{'_cov' if covariates else ''}{'_shuffledpheno' if shuffled_phenotype else ''}{f'_{test_size}' if test_size else ''}{'_test' if only_test else ''}_NEWsplit_{trait}.pq"
    
    if os.path.isfile(outfile):
        return
        
    (gt_train, gt_val, gt_test), (residual_train, residual_val, residual_test), _, gene_list, _, (trait_measurement_train, trait_measurement_val, trait_measurement_test), (covariates_train, covariates_val, covariates_test) = dataloader_clean.load_data(trait, embedding_type=None, use_prs=use_prs, normalize_covariates=normalize_covariates, version=version, genotype=genotype, test_split_size=test_size, split_seed=0)

    print('Train size: ', trait_measurement_train.shape[0])
    print('Val size: ', trait_measurement_val.shape[0])
    print('Test size: ', trait_measurement_test.shape[0])
    
    if shuffled_phenotype:
        np_random_seed = 1234
        np.random.seed(np_random_seed)
        trait_measurement_train = np.random.permutation(trait_measurement_train)
        np_random_seed = 1234
        np.random.seed(np_random_seed)
        trait_measurement_val = np.random.permutation(trait_measurement_val)
        np_random_seed = 1234
        np.random.seed(np_random_seed)
        y_train = np.random.permutation(y_train)
        np_random_seed = 1234
        np.random.seed(np_random_seed)
        y_val = np.random.permutation(y_val)
    
    #if covariates:
    #    y_train = trait_measurement_train
    #    y_val = trait_measurement_val
    #    y_test = trait_measurement_test
    
    betas = []
    std_err = []
    pvals = []
    if only_test:
        y = np.concatenate([residual_test])
        gt = np.concatenate([gt_test]).astype(float)
        c = np.concatenate([covariates_test])
    else:
        y = np.concatenate([residual_train, residual_val])
        gt = np.concatenate([gt_train, gt_val]).astype(float)
        c = np.concatenate([covariates_train, covariates_val])
    
    for i in tqdm(range(len(gene_list))):
        if gt[:,i].sum()==0:
            betas.append(0)
            std_err.append(0)
            pvals.append(1)
        else:
            X = sm.add_constant(gt[:,i])
            if covariates:
                X = np.concatenate([X, c], axis=1)
            results = sm.OLS(y, X).fit()
            betas.append(results.params[1])
            std_err.append(results.bse[1])
            pvals.append(results.pvalues[1])

    at_df = pd.DataFrame({"trait": trait, "model": f"ols_{genotype}", "beta": betas, "std_err": std_err, "pval": pvals}, index=gene_list)
    at_df["pval"] = at_df["pval"].fillna(1)
    at_df["neglog_pval"] = -np.log10(at_df["pval"])
    at_df['test_split_size'] = test_size

    at_df.to_parquet(outfile)

if __name__ == "__main__":
    main()
