import yaml
import click
from typing import Optional
import pandas as pd
import numpy as np

from ..utils import dataloader_old, utils
import os
import sys
import statsmodels.api as sm
from tqdm import tqdm

import logging

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

@click.group()
def cli():
    pass

@cli.command()
@click.option("--trait", type=str, required=True, help="Trait to train model for")
@click.option("--config-path", type=click.Path(exists=True), required=True, help="Config file with all details")
@click.option("--trainval-sampling", type=float, required=False, help="How much should we sample from the training and validation set")
@click.option("--only-test", is_flag=True, default=False, help="Perform the association test only on the test data")
@click.option("--bootstrap", is_flag=True, default=False, help="Should we bootstrap the data (use with only-test flag)")
@click.option("--bs-seed", type=int, required=False, help="bootstrap seed")
@click.option("--gene-subset-file", type=str, required=False, help="file with the gene subset to run for a given trait")
def association_test(
    trait: str, 
    config_path: str, 
    trainval_sampling: Optional[float] = None,
    only_test: Optional[bool] = False,
    bootstrap: Optional[bool] = False,
    bs_seed: Optional[int] = 0,
    gene_subset_file: Optional[str] = None,
):
    # --- Run Configuration Loading (load_data parameters) ---
    if not os.path.exists(config_path):
        logger.error(f"Run configuration file not found: {config_path}. Exiting.") # Warning instead of error, defaults will be used.
        sys.exit(1)
    else:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

    # --- Load Run Parameters from Config ---
    dataloader_params = config.get('dataloader_params', None)
    if dataloader_params.get('genotype_path', None):
        genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]

    # --- Create Experiment Directory ---
    if config.get('output_dir_base', None):
        ols_op_dir = os.path.join(config['output_dir_base'], genotype, f"rvat{'_onlytest' if only_test else ''}{'_sampling' if trainval_sampling else ''}{trainval_sampling if trainval_sampling else ''}")
    else:
        logger.error("Output directory not defined in config. Exiting.")
        sys.exit(1)
    os.makedirs(ols_op_dir, exist_ok=True)

    config_copy_path = os.path.join(ols_op_dir, "run_config.yaml")
    with open(config_copy_path, "w") as f:
        yaml.dump(config, f, indent=2)  # Save the entire config
    logger.info(f"Run configuration saved to: {config_copy_path}")

    outfile = os.path.join(ols_op_dir, f"{trait}_rvat{'_bootstrap' if bootstrap else ''}{bs_seed if bootstrap else ''}.pq")
    if os.path.isfile(outfile):
        logger.info(f"Output file already exists: {outfile}. Exiting.")
        return

    gene_subset = dataloader_params.get('gene_subset', None)
    if gene_subset_file:
        logger.info(f"Using gene subset from file {gene_subset_file}")
        gene_subset = list(pd.read_parquet(gene_subset_file).query("trait == @trait")['gene_id'])

    (
        (gt_train, gt_val, gt_test),
        _,
        _,
        gene_list,
        _,
        (trait_measurement_train, trait_measurement_val, trait_measurement_test),
        (covariates_train, covariates_val, covariates_test),
    ) = dataloader_old.load_data(
        trait,
        embedding_path = None,
        genotype_path = dataloader_params.get('genotype_path', None),
        phenotype_dir= dataloader_params.get('phenotype_dir', None),
        covariates_path = dataloader_params.get('covariates_path', None),
        prs_path = dataloader_params.get('prs_path', None),
        train_individuals_path = dataloader_params.get('train_individuals_path', '-1'),
        test_split_size = dataloader_params.get('test_split_size', 0.25),
        val_split_size = dataloader_params.get('val_split_size', 0.1),
        trainval_sampling = trainval_sampling,
        split_seed = dataloader_params.get('split_seed', 0),
        use_prs = dataloader_params.get('use_prs', True),
        normalize_covariates = dataloader_params.get('normalize_covariates', True),
        normalize_embedding=dataloader_params.get('normalize_embedding', False),
        shuffled_phenotype = dataloader_params.get('shuffled_phenotype', False),
        shuffled_embedding = dataloader_params.get('shuffled_embedding', False),
        random_embedding = dataloader_params.get('random_embedding', False),
        gene_subset = gene_subset,
        dataset_version = dataloader_params.get('dataset_version', "filteredv3"),
    )

    logger.info(f"Train size: {trait_measurement_train.shape[0]}")
    logger.info(f"Val size: {trait_measurement_val.shape[0]}")
    logger.info(f"Test size: {trait_measurement_test.shape[0]}")

    if only_test:
        if bootstrap:
            # Generate random indices with replacement
            np.random.seed(bs_seed)
            indices = np.random.choice(gt_test.shape[0], size=gt_test.shape[0], replace=True)

            # Sample rows from each array using the same indices
            gt = gt_test[indices].astype(float)
            c = covariates_test[indices]
            y = trait_measurement_test[indices]
        else:
            y = np.concatenate([trait_measurement_test])
            gt = np.concatenate([gt_test]).astype(float)
            c = np.concatenate([covariates_test])
    else:
        y = np.concatenate([trait_measurement_train, trait_measurement_val])
        gt = np.concatenate([gt_train, gt_val]).astype(float)
        c = np.concatenate([covariates_train, covariates_val])

    logger.info("Starting association test...")
    betas = []
    std_err = []
    pvals = []
    # --- Association Test ---    
    for i in tqdm(range(len(gene_list))):
        if gt[:,i].sum()==0:
            betas.append(0)
            std_err.append(0)
            pvals.append(1)
        else:
            X = sm.add_constant(gt[:,i])
            X = np.concatenate([X, c], axis=1)

            results = sm.OLS(y, X).fit()
            betas.append(results.params[1])
            std_err.append(results.bse[1])
            pvals.append(results.pvalues[1])

    at_df = pd.DataFrame(
        {
            "trait": trait, 
            "beta": betas, 
            "std_err": std_err, 
            "pval": pvals
        }, index=gene_list)
    at_df = at_df.reset_index().rename(columns={"index": "gene_id"})

    at_df["pval"] = at_df["pval"].fillna(1)
    at_df["neglog_pval"] = -np.log10(at_df["pval"])
    at_df['test_split_size'] = dataloader_params.get('test_split_size', 0.25)
    at_df['model'] = 'ols'
    at_df['genotype'] = genotype
    at_df["dataset_version"] = dataloader_params.get('dataset_version', "filteredv3")

    at_df.to_parquet(outfile)
    logger.info(f"Output saved to {outfile}")

@cli.command()
@click.option("--trait", type=str, required=True, help="Trait to train model for")
@click.option("--config-path", type=click.Path(exists=True), required=True, help="Config file with all details")
@click.option("--pval-threshold", type=float, default=0.05, help="p-value threshold for significance")
@click.option("--only-covariates", is_flag=True, default=False, help="Use only covariates for phenotype prediction")
@click.option("--trainval-sampling", type=float, required=False, help="How much should we sample from the training and validation set")
@click.pass_context  # <--- ADD THIS DECORATOR
def lm_phenopred(
    ctx,
    trait: str,
    config_path: str,
    pval_threshold: float = 0.05,
    only_covariates: Optional[bool] = False,
    trainval_sampling: Optional[float] = None,
):    
    # --- Run Configuration Loading (load_data parameters) ---
    if not os.path.exists(config_path):
        logger.error(f"Run configuration file not found: {config_path}. Exiting.") # Warning instead of error, defaults will be used.
        sys.exit(1)
    else:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

    # if trainval_sampling <= 0 or trainval_sampling >= 1:
    #     trainval_sampling = False
    # if trainval_sampling == None:
    #     trainval_sampling = False

    # --- Load Run Parameters from Config ---
    dataloader_params = config.get('dataloader_params', None)
    if dataloader_params.get('genotype_path', None):
        genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]

    # --- Create Experiment Directory ---
    if config.get('output_dir_base', None):
        # ols_op_dir = os.path.join(config['output_dir_base'], genotype, f"rvat")
        ols_op_dir = os.path.join(config['output_dir_base'], genotype, f"rvat{'_sampling' if trainval_sampling else ''}{trainval_sampling if trainval_sampling else ''}{'_onlycov' if only_covariates else ''}")
    else:
        logger.error("Output directory not defined in config. Exiting.")
        sys.exit(1)
    os.makedirs(ols_op_dir, exist_ok=True)
    
    config_copy_path = os.path.join(ols_op_dir, "run_config.yaml")
    with open(config_copy_path, "w") as f:
        yaml.dump(config, f, indent=2)  # Save the entire config
    logger.info(f"Run configuration saved to: {config_copy_path}")

    outfile = os.path.join(ols_op_dir, f"{trait}_phenopred_{pval_threshold}nom.pq")
    if os.path.isfile(outfile):
        logger.info(f"Output file already exists: {outfile}. Exiting.")
        return

    # --- Check if input file exists ---
    infile = os.path.join(ols_op_dir, f"{trait}_rvat.pq")
    if not os.path.isfile(infile) and only_covariates==False:
        logger.info(f"Input file not found: {infile}. Running association test first.")
        # association_test(trait, config_path)
        ctx.invoke(
                association_test,
                trait=trait,
                config_path=config_path,
                only_test=False,  # Explicitly set to False
                trainval_sampling=trainval_sampling
            )
    else:
        logger.info(f"Input file found: {infile}. Proceeding with phenotype prediction.")

    # --- Load RVAT results ---
    signif_gene_list = ['ENSG00000000419']  # Just some gene - we will not use it anyways
    if only_covariates==False:
        rvat_df = pd.read_parquet(infile)
        # rvat_df['significant'] = rvat_df["pval"]*rvat_df.shape[0] < pval_threshold
        rvat_df['significant'] = rvat_df["pval"] < pval_threshold
        signif_gene_list = list(rvat_df[rvat_df['significant']]['gene_id'])
        logger.info(f"Trait: {trait}, p_theshold: {pval_threshold}, # significant genes: {rvat_df['significant'].sum()}")

    (
        (gt_train, gt_val, gt_test),
        (y_train_residual, y_val_residual, y_test_residual),
        _,
        gene_list,
        (id_train, id_val, id_test),
        (trait_measurement_train, trait_measurement_val, trait_measurement_test),
        (covariates_train, covariates_val, covariates_test),
    ) = dataloader_old.load_data(
        trait,
        embedding_path = None,
        genotype_path = dataloader_params.get('genotype_path', None),
        phenotype_dir= dataloader_params.get('phenotype_dir', None),
        covariates_path = dataloader_params.get('covariates_path', None),
        prs_path = dataloader_params.get('prs_path', None),
        train_individuals_path = dataloader_params.get('train_individuals_path', '-1'),
        test_split_size = dataloader_params.get('test_split_size', 0.25),
        val_split_size = dataloader_params.get('val_split_size', 0.1),
        trainval_sampling = trainval_sampling,
        split_seed = dataloader_params.get('split_seed', 0),
        use_prs = dataloader_params.get('use_prs', True),
        normalize_covariates = dataloader_params.get('normalize_covariates', True),
        normalize_embedding=dataloader_params.get('normalize_embedding', False),
        shuffled_phenotype = dataloader_params.get('shuffled_phenotype', False),
        shuffled_embedding = dataloader_params.get('shuffled_embedding', False),
        random_embedding = dataloader_params.get('random_embedding', False),
        gene_subset = signif_gene_list, # Only use significant genes
        dataset_version = dataloader_params.get('dataset_version', "filteredv3"),
    )

    logger.info(f"Train size: {trait_measurement_train.shape[0]}")
    logger.info(f"Val size: {trait_measurement_val.shape[0]}")
    logger.info(f"Test size: {trait_measurement_test.shape[0]}")

    y = np.concatenate([trait_measurement_train, trait_measurement_val])
    c = np.concatenate([covariates_train, covariates_val])

    logger.info("Starting phenotype prediction on test...")

    if only_covariates:
        logger.info("Using only covariates for phenotype prediction")
        X = sm.add_constant(c)
        X_test = sm.add_constant(covariates_test)

    else:
        logger.info("Using significant genes and covariates for phenotype prediction")
        gt = np.concatenate([gt_train, gt_val]).astype(float)
        X = sm.add_constant(gt)
        X = np.concatenate([X, c], axis=1)
        X_test = sm.add_constant(gt_test)
        X_test = np.concatenate([X_test, covariates_test], axis=1)

    lm_model = sm.OLS(y, X).fit()

    pred_df = pd.DataFrame(
        {
            "trait_measurement": trait_measurement_test, 
            "common_variant_residual": y_test_residual, 
            "best_prediction": lm_model.predict(X_test), 
            "trait": trait, 
            "model": f"lm_sign_genes",
            "genotype": genotype,
            "dataset_version": dataloader_params.get('dataset_version', "filteredv3"),
            "only_covariates": only_covariates,
        }, index=id_test).reset_index()

    logger.info("Finished phenotype prediction on test...")

    pred_df.to_parquet(outfile)
    logger.info(f"Output saved to {outfile}")

if __name__ == "__main__":
    cli()