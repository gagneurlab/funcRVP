import logging
import sys
import os
import yaml
import click
from tqdm import tqdm
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.preprocessing import quantile_transform

import torch
import optuna

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def trait_INT_zscore(trait, phenotype_data_dir):
    pheno_dt = pd.read_parquet(
        f"{phenotype_data_dir}/{trait}/data.parquet"
    )
    trait_mean = pheno_dt[f"{trait}_raw"].mean()
    trait_sd = pheno_dt[f"{trait}_raw"].std()
    pheno_dt["zscore"] = (pheno_dt[f"{trait}_raw"] - trait_mean) / trait_sd
    pheno_dt["raw"] = pheno_dt[f"{trait}_raw"]
    pheno_dt["INT"] = pheno_dt[f"{trait}"]
    pheno_dt["new_INT"] = quantile_transform(
        pheno_dt[["zscore"]], output_distribution="normal", random_state=0, copy=True
    )
    pheno_dt["trait"] = trait

    pheno_dt = pheno_dt.rename(columns={"eid": "individual"})
    # pheno_dt["individual"] = pheno_dt["individual"].astype("str")
    # pheno_dt = pheno_dt.set_index('individual')
    return pheno_dt[["individual", "trait", "raw", "zscore", "INT", "new_INT"]]


def get_best_arch(
    trait,
    study_version="v108cov_deepRVAT",
    embedding_type="omics_pops",
    storage_path="/s/project/geno2pheno/hyperopt/optuna/journal.log",  # Added storage_path as argument and default
):
    """
    Retrieves the best architecture parameters from an Optuna study.

    Args:
        trait (str): Trait name.
        study_version (str): Study version name.
        embedding_type (str): Embedding type.
        storage_path (str): Path to Optuna storage journal log file.

    Returns:
        dict: Dictionary containing best trial parameters.
    """
    study_name = f"study_{study_version}_{trait}_{embedding_type}"
    storage = optuna.storages.JournalStorage(
        optuna.storages.JournalFileStorage(storage_path),  # Use storage_path argument
    )
    try:
        study = optuna.load_study(study_name=study_name, storage=storage)
    except Exception as e:  # Catch potential errors during study loading
        logger.error(f"Error loading Optuna study '{study_name}': {e}")
        sys.exit(1)  # Exit if study loading fails

    if study._is_multi_objective():
        best_trial = max(study.best_trials, key=lambda t: t.values[1])
    else:
        best_trial = study.best_trial

    best_trial_params = {
        "trial_number": best_trial.number,
        "r2": best_trial.values[0],
        "n_total_trials": len(study.trials_dataframe().query("state=='COMPLETE'")),
        "n_hidden": best_trial.params["n_hidden"],
        "hidden_dim": best_trial.params["hidden_dim"],
        "last_layer_bias": best_trial.params["last_layer_bias"],
        "alpha_L1_fE": best_trial.params["alpha_L1_fE"],
        "normalize_embeddings": best_trial.params["normalize_embeddings"],
    }

    return best_trial_params

def save_model_betas(
    g2p_cov_model,
    gene_list,
    trait,
    output_dir,
    config,
):
    """
    Saves model outputs (mean betas and bayes predictions) to parquet files
    and saves a copy of the run configuration.

    Args:
        g2p_cov_model: Trained G2P model.
        gene_list: List of gene names.
        trait: Trait name.
        output_dir: Directory to save outputs in.
        config: Run configuration dictionary.
    """

    logger.info("Saving model outputs...")

    experiment_name = config.get("experiment_name", "default")
    dataloader_params = config.get("dataloader_params", None)
    genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]
    embedding = dataloader_params.get('embedding_path', None).split("/")[-1].split(".")[0]
    dataset_version = dataloader_params.get("dataset_version", "filteredv3")

    config_copy_path = os.path.join(output_dir, "run_config.yaml")
    with open(config_copy_path, "w") as f:
        yaml.dump(config, f, indent=2)  # Save the entire config
    logger.info(f"Run configuration saved to: {config_copy_path}")

    # Save posterior betas
    with torch.no_grad(): # Disable gradient calculations
        betas_df = pd.DataFrame(
            {
                "posterior_beta": g2p_cov_model.best_posterior_mean_beta,
                "posterior_beta_se": np.sqrt(g2p_cov_model.best_posterior_var_beta),
                "prior_var": g2p_cov_model.best_prior_var.flatten(),
                "intercept": g2p_cov_model.best_intercept,
                "y_var": g2p_cov_model.best_var,
                "last_layer_bias": g2p_cov_model.best_last_layer_bias,
            },
            index=gene_list,
        )
        betas_df = betas_df.reset_index().rename(columns={"index": "gene_id"})

        gene_names_path = config["hgnc_gene_names"]
        if not os.path.exists(gene_names_path):
            logger.error(f"Gene names file not found at {gene_names_path}")
            logger.error(f"Skipping gene name merge")
        else:
            logger.info(f"Loading gene names from: {gene_names_path}")
            gene_names = (
                pd.read_csv(gene_names_path, sep="\t")[
                    ["Ensembl gene ID", "Approved symbol"]
                ]
                .drop_duplicates()
                .rename(
                    columns={"Ensembl gene ID": "gene_id", "Approved symbol": "gene_name"}
                )
            )
            betas_df = betas_df.merge(gene_names, on="gene_id")

        betas_df["pd"] = np.maximum(
            scipy.stats.norm.cdf(0, betas_df["posterior_beta"], betas_df["posterior_beta_se"]),
            scipy.stats.norm.sf(0, betas_df["posterior_beta"], betas_df["posterior_beta_se"]),
        )
        betas_df["neglog_pval"] = -np.log10(1 - betas_df["pd"])
        betas_df["significant"] = betas_df["pd"] > config.get("pd_signif_threshold", 0.999)
        betas_df["trait"] = trait
        betas_df["embedding"] = embedding
        betas_df["genotype"] = genotype
        betas_df["model"] = 'funcrvp'
        betas_df["dataset_version"] = dataset_version
        betas_df["experiment_name"] = experiment_name

        output_mean_betas_path = os.path.join(output_dir, f"{trait}_betas.pq")
        logger.info(f"Saving betas to: {output_mean_betas_path}")
        betas_df.to_parquet(output_mean_betas_path)


def save_model_predictions(
    trait_measurement_test,
    y_test_residual,
    best_model_test_pred,
    test_ids,
    trait,
    output_dir,
    config,
):
    """
    Saves model outputs (mean betas and bayes predictions) to parquet files
    and saves a copy of the run configuration.

    Args:
        g2p_cov_model: Trained G2P model.
        gene_list: List of gene names.
        trait: Trait name.
        output_dir: Directory to save outputs in.
        config: Run configuration dictionary.
    """

    logger.info("Saving model outputs...")

    experiment_name = config.get("experiment_name", "default")
    dataloader_params = config.get("dataloader_params", None)
    genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]
    embedding = dataloader_params.get('embedding_path', None).split("/")[-1].split(".")[0]
    dataset_version = dataloader_params.get("dataset_version", "filteredv3")

    config_copy_path = os.path.join(output_dir, "run_config.yaml")
    with open(config_copy_path, "w") as f:
        yaml.dump(config, f, indent=2)  # Save the entire config
    logger.info(f"Run configuration saved to: {config_copy_path}")

    # Save model predictions on test
    with torch.no_grad(): # Disable gradient calculations
        # Save bayes predictions
        phenopred_df = pd.DataFrame(
            {
                "trait_measurement": trait_measurement_test,
                "common_variant_residual": y_test_residual,
                "best_prediction": best_model_test_pred,
            },
            index=test_ids,
        )
        phenopred_df["trait"] = trait
        phenopred_df["model"] = 'funcrvp'
        phenopred_df["embedding"] = embedding
        phenopred_df["genotype"] = genotype
        phenopred_df["dataset_version"] = dataset_version
        phenopred_df["experiment_name"] = experiment_name

        output_bayes_pred_path = os.path.join(output_dir, f"{trait}_phenopred.pq")
        logger.info(f"Saving phenotype predictions to: {output_bayes_pred_path}")
        phenopred_df.to_parquet(output_bayes_pred_path)

        logger.info("All model outputs saved.")


def save_model_outputs(
    g2p_cov_model,
    trait_measurement_test,
    y_test_residual,
    best_model_test_pred,
    test_ids,
    gene_list,
    trait,
    output_dir,
    config,
):
    """
    Saves model outputs (mean betas and bayes predictions) to parquet files
    and saves a copy of the run configuration.

    Args:
        g2p_cov_model: Trained G2P model.
        gene_list: List of gene names.
        trait: Trait name.
        output_dir: Directory to save outputs in.
        config: Run configuration dictionary.
    """

    logger.info("Saving model outputs...")

    experiment_name = config.get("experiment_name", "default")
    dataloader_params = config.get("dataloader_params", None)
    genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]
    embedding = dataloader_params.get('embedding_path', None).split("/")[-1].split(".")[0]
    dataset_version = dataloader_params.get("dataset_version", "filteredv3")

    config_copy_path = os.path.join(output_dir, "run_config.yaml")
    with open(config_copy_path, "w") as f:
        yaml.dump(config, f, indent=2)  # Save the entire config
    logger.info(f"Run configuration saved to: {config_copy_path}")

    # Save posterior betas
    with torch.no_grad(): # Disable gradient calculations
        betas_df = pd.DataFrame(
            {
                "posterior_beta": g2p_cov_model.best_posterior_mean_beta,
                "posterior_beta_se": np.sqrt(g2p_cov_model.best_posterior_var_beta),
                "prior_var": g2p_cov_model.best_prior_var.flatten(),
                "intercept": g2p_cov_model.best_intercept,
                "y_var": g2p_cov_model.best_var,
                "last_layer_bias": g2p_cov_model.best_last_layer_bias,
            },
            index=gene_list,
        )
        betas_df = betas_df.reset_index().rename(columns={"index": "gene_id"})

        gene_names_path = config["hgnc_gene_names"]
        if not os.path.exists(gene_names_path):
            logger.error(f"Gene names file not found at {gene_names_path}")
            logger.error(f"Skipping gene name merge")
        else:
            logger.info(f"Loading gene names from: {gene_names_path}")
            gene_names = (
                pd.read_csv(gene_names_path, sep="\t")[
                    ["Ensembl gene ID", "Approved symbol"]
                ]
                .drop_duplicates()
                .rename(
                    columns={"Ensembl gene ID": "gene_id", "Approved symbol": "gene_name"}
                )
            )
            betas_df = betas_df.merge(gene_names, on="gene_id")

        betas_df["pd"] = np.maximum(
            scipy.stats.norm.cdf(0, betas_df["posterior_beta"], betas_df["posterior_beta_se"]),
            scipy.stats.norm.sf(0, betas_df["posterior_beta"], betas_df["posterior_beta_se"]),
        )
        betas_df["neglog_pval"] = -np.log10(1 - betas_df["pd"])
        betas_df["significant"] = betas_df["pd"] > config.get("pd_signif_threshold", 0.999)
        betas_df["trait"] = trait
        betas_df["embedding"] = embedding
        betas_df["genotype"] = genotype
        betas_df["model"] = 'funcrvp'
        betas_df["dataset_version"] = dataset_version
        betas_df["experiment_name"] = experiment_name

        output_mean_betas_path = os.path.join(output_dir, f"{trait}_betas.pq")
        logger.info(f"Saving betas to: {output_mean_betas_path}")
        betas_df.to_parquet(output_mean_betas_path)

        # Save bayes predictions
        phenopred_df = pd.DataFrame(
            {
                "trait_measurement": trait_measurement_test,
                "common_variant_residual": y_test_residual,
                "best_prediction": best_model_test_pred,
            },
            index=test_ids,
        )
        phenopred_df["trait"] = trait
        phenopred_df["model"] = 'funcrvp'
        phenopred_df["embedding"] = embedding
        phenopred_df["genotype"] = genotype
        phenopred_df["dataset_version"] = dataset_version
        phenopred_df["experiment_name"] = experiment_name

        output_bayes_pred_path = os.path.join(output_dir, f"{trait}_phenopred.pq")
        logger.info(f"Saving phenotype predictions to: {output_bayes_pred_path}")
        phenopred_df.to_parquet(output_bayes_pred_path)

        logger.info("All model outputs saved.")


@click.group()
def cli():
    pass

@cli.command()
@click.option(
    "--config-path",
    type=click.Path(exists=True),
    required=True,
    help="Config file with all details",
)
def read_results(
    config_path: Path,
):
    if not os.path.exists(config_path):
        logger.warning(
            f"Run configuration file not found: {config_path}. Exiting..."
        )
        sys.exit(1)
    else:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

    TRAITS = config.get("traits", None)
    experiment_name = config.get("experiment_name", "default")

    dataloader_params = config.get("dataloader_params", None)
    dataset_version = dataloader_params.get("dataset_version", "filteredv3")
    embedding = dataloader_params.get("embedding", None)
    genotype = dataloader_params.get("genotype", "plof")

    output_dir_base = config.get("output_dir_base", None)
    if not output_dir_base:
        logger.error("Output directory not provided in config file.")
        sys.exit(1)

    # Define output directory
    output_dir_name = f"{experiment_name}_{dataset_version}/{genotype}_{embedding}"
    output_dir = os.path.join(config["output_dir_base"], output_dir_name)

    genes_dt_list = []
    pheno_dt_list = []
    skip_list = []

    logger.info(f"Starting to consolidate results for all traits...")
    for trait in tqdm(TRAITS):
        try:
            # TODO: Add phenocode mapping
            # phenocode = genebass_phenocode_dict[trait]

            betas_df = pd.read_parquet(
                f"{output_dir}/{trait}_betas.pq"
            )  # .reset_index()
            # betas_df["phenocode"] = str(phenocode)
            genes_dt_list.append(betas_df)

            phenopred_df = pd.read_parquet(
                f"{output_dir}/{trait}_phenopred.pq"
            )  # .reset_index()
            # phenopred_df["phenocode"] = str(phenocode)
            pheno_dt_list.append(phenopred_df)
        except:
            skip_list.append(trait)
            continue

    logger.info(f"Skipped the traits: {skip_list}")

    genes_dt = pd.concat(genes_dt_list)
    pheno_dt = pd.concat(pheno_dt_list)

    logger.info(f"Writing consolidated results files to: {output_dir}")
    genes_dt.to_parquet(f"{output_dir}/all_trait_betas.pq", index=False)
    pheno_dt.to_parquet(f"{output_dir}/all_trait_phenopred.pq", index=False)
    logger.info(f"Consolidated results saved successfully.")


def load_phenocode_dict(prs_score_mapping_path):
    if not os.path.exists(prs_score_mapping_path):
        logger.error(f"PRS score mapping file not found: {prs_score_mapping_path}")
        sys.exit(1)
    logger.info(f"Loading phenocode dictionary from: {prs_score_mapping_path}")
    genebass_phenocode_dict = (
        pd.read_csv(prs_score_mapping_path)[["phenotype", "genebass_phenocode"]]
        .set_index("genebass_phenocode")["phenotype"]
        .to_dict()
    )
    genebass_phenocode_dict = {v: k for k, v in genebass_phenocode_dict.items()}
    genebass_phenocode_dict["LDL_direct"] = "30780"
    genebass_phenocode_dict["Albumin"] = "30600"
    return genebass_phenocode_dict