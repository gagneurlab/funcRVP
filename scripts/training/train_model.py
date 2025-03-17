import os
import sys
import yaml
import logging

import click
from typing import Optional

import pandas as pd
import numpy as np

from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression, ElasticNetCV
import matplotlib.pyplot as plt
import torch, gc

from pathlib import Path
import wandb
import optuna

from ..utils import dataloader, utils
from ..models import g2p_bayes_cov_skipcon
g2p_bayes_model = g2p_bayes_cov_skipcon


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

def trainer_(
    trait: str,
    emb: np.ndarray,
    gt_train: np.ndarray,
    covariates_train: np.ndarray,
    trait_measurement_train: np.ndarray,
    gt_val: np.ndarray,
    covariates_val: np.ndarray,
    trait_measurement_val: np.ndarray,
    n_repeats: Optional[int]=1,
    best_n_hidden: Optional[int]=1,
    best_model_arch: Optional[dict]=None,
    y_var_init: Optional[float]=1e-3,
    base_var_init: Optional[float]=1e-5,
    training_params: Optional[dict]=None,
):
    """
    Trains the G2P model for n_repeats with a fixed architecture.

    Args:
        trait (str): Trait name.
        emb (np.ndarray): Gene embedding matrix.
        gt_train (np.ndarray): Training genotype matrix.
        covariates_train (np.ndarray): Training covariate matrix.
        trait_measurement_train (np.ndarray): Training trait measurements.
        gt_val (np.ndarray): Validation genotype matrix.
        covariates_val (np.ndarray): Validation covariate matrix.
        trait_measurement_val (np.ndarray): Validation trait measurements.
        gene_list (list): List of gene names.
        best_n_hidden (int): Number of hidden layers in f(E).
        best_model_arch (dict): Dictionary containing best architecture parameters.
        base_var_init (float): Base variance initialization.
        training_params (dict): Training parameters from config.
        config (dict): Entire run configuration.

    Returns:
        g2p_bayes_model.G2P_Model: Trained G2P model.
    """
    for rep in range(n_repeats):

        logger.info(f"Training model for {trait} - Repetition {rep+1}/{n_repeats}")
        # Initialize model and send to device
        g2p_cov_model = g2p_bayes_model.G2P_Model(
            emb.shape[1], # embedding_dim
            covariates_train.shape[1],
            n_hidden=best_n_hidden,
            hiddem_dim=best_model_arch["hidden_dim"],
            last_layer_bias=best_model_arch["last_layer_bias"],
            nonlinearity=training_params['nonlinearity'],
            y_var_init=y_var_init,
            base_var_init=base_var_init,
            base_var_const=0, 
            n_genes=emb.shape[0],
            device=device,
        ).to(device)

        # Initialize weights and biases
        logger.info(f"Writing wandb logs to {training_params['wandb_dir']}")

        wandb.init(
            project=training_params['wandb_project_name'],
            # name=training_params['wandb_project_name'],
            config=training_params,
            dir=training_params['wandb_dir'],
            settings=wandb.Settings(_service_wait=600),
        )

        wandb.watch(g2p_cov_model)

        logger.info("Starting model training...")
        g2p_cov_model.fit_model(
            training_params, # Pass the entire config dictionary to fit_model
            emb,
            gt_train,
            covariates_train,
            trait_measurement_train,
            G_val=gt_val,
            C_val=covariates_val,
            y_val=trait_measurement_val,
            logging=True,
            device=device,
        )
        logger.info("Model training finished.")

        wandb.finish()

    return g2p_cov_model


@click.group()
def cli():
    pass

@cli.command()
@click.option("--trait", type=str, required=True, help="Trait to train model for")
@click.option("--config-path", type=click.Path(exists=True), required=True, help="Config file with all details")
def trainer(
    trait: str,
    config_path: Path,
):
    # TODO! Define a Minimum config file.
    # --- Run Configuration Loading (load_data parameters) ---
    if not os.path.exists(config_path):
        logger.warning(f"Run configuration file not found: {config_path}. Using default parameters.") # Warning instead of error, defaults will be used.
        config = {} # Empty dict if run_config.yaml is missing
    else:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

    # --- Load Run Parameters from Config ---
    dataloader_params = config.get('dataloader_params', None)
    embedding = dataloader_params.get('embedding', None)
    genotype = dataloader_params.get('genotype', 'plof')

    (
        (gt_train, gt_val, gt_test),
        (y_train_residual, y_val_residual, y_test_residual),
        emb,
        gene_list,
        (id_train, id_val, id_test),
        (trait_measurement_train, trait_measurement_val, trait_measurement_test),
        (covariates_train, covariates_val, covariates_test),
    ) = dataloader.load_data(
        trait,
        embedding_type = embedding,
        genotype = genotype,
        test_split_size = dataloader_params.get('test_split_size', 0.25),
        val_split_size = dataloader_params.get('val_split_size', 0.1),
        split_seed = dataloader_params.get('split_seed', 0),
        gene_subset = dataloader_params.get('gene_subset', None),
        use_prs = dataloader_params.get('use_prs', True),
        normalize_covariates = dataloader_params.get('normalize_covariates', True),
        normalize_embedding=dataloader_params.get('normalize_embedding', False),
        shuffled_phenotype = dataloader_params.get('shuffled_phenotype', False),
        shuffled_embedding = dataloader_params.get('shuffled_embedding', False),
        random_embedding = dataloader_params.get('random_embedding', False),
        dataset_version = dataloader_params.get('dataset_version', "filteredv3"),
    )

    # Read Training and hyperopt params from config
    training_params = config.get('training_params', None)
    training_params['trait'] = trait
    training_params['embedding'] = embedding
    training_params['genotype'] = genotype

    hpopt_params = config.get('hpopt_params', None)

    # If no hyperopt, get best model architecture
    hpopt_params = None
    if hpopt_params is None:
        n_repeats = training_params.get('n_repeats', 1)

        best_model_arch = utils.get_best_arch(
            trait=trait,
            study_version=training_params['old_model_version'], # Access model_params
            embedding_type=training_params['old_model_embedding'], # Access model_params
            storage_path=training_params['old_optuna_journal_log'] # Access paths through 'config'
        )

        # TODO! should we specify below numbers in the config file?
        # To get constant regularization
        if embedding == None:
            best_n_hidden = -1
            base_var_init = 5e-3
        else:
            best_n_hidden = best_model_arch["n_hidden"]  # number of hidden layers in f(E)
            base_var_init = 5e-5

        g2p_cov_model = trainer_( # Call trainer_ to train with fixed arch
            trait=trait,
            emb=emb,
            gt_train=gt_train,
            covariates_train=covariates_train,
            trait_measurement_train=trait_measurement_train,
            gt_val=gt_val,
            covariates_val=covariates_val,
            trait_measurement_val=trait_measurement_val,
            n_repeats=1,
            best_n_hidden=best_n_hidden,
            best_model_arch=best_model_arch,
            y_var_init=y_train_residual.var(),
            base_var_init=base_var_init,
            training_params=training_params,
        )

    # TODO fix optuna suggetions for hyperparams
    else:
        # TODO Get optuna suggestions for hyperparams
        logger.info("Hyperparameter optimization is under development.")
        exit()

        logger.info("Hyperparameter optimization enabled.")

        # Hyperparameter optimization settings from config
        n_repeats = hpopt_params.get('n_repeats', 1)
        n_trials = hpopt_params.get('n_trials', 1)
        objective_direction = hpopt_params.get("direction", "maximize")
        hpopt_file = hpopt_params.get("hpopt_file", None)

        sampler = optuna.samplers.TPESampler(multivariate = True, n_startup_trials=10)

        # Not in a journal
        if hpopt_file is None:
            study = optuna.create_study(
                direction=objective_direction,
                sampler=sampler,
                )
        # In an existing journal
        else:
            # Optuna journal path
            storage = optuna.storages.JournalStorage(
                optuna.storages.journal.JournalFileBackend(hpopt_file)
                )
            study = optuna.create_study(
                study_name=Path(hpopt_file).stem,
                storage=storage,
                load_if_exists=True,
                direction=objective_direction,
                sampler=sampler,
                )

        trials_done = len(study.trials)
        for i in range(n_trials):
            trial_num = trials_done + i
            logger.info(f"Starting optuna trial: {trial_num}")

            study.optimize(
                trainer_( # Call trainer_ to train with fixed arch
                    trait=trait,
                    emb=emb,
                    gt_train=gt_train,
                    covariates_train=covariates_train,
                    trait_measurement_train=trait_measurement_train,
                    gt_val=gt_val,
                    covariates_val=covariates_val,
                    trait_measurement_val=trait_measurement_val,
                    n_repeats=n_repeats,
                    best_n_hidden=best_n_hidden,
                    best_model_arch=best_model_arch,
                    y_var_init=y_train_residual.var(),
                    base_var_init=base_var_init,
                    training_params=training_params,
                ),  
            n_trials=1,
            gc_after_trial=True,
            )
        
        # TODO Get best model from the hyperopt then save model outputs, or retrain model through a recursive call and then save model outputs...

        logger.info(f"Number of finished trials: {len(study.trials)}")

    # --- Create Experiment Directory ---
    output_dir_name = f"{config.get('experiment_name', 'default')}_{dataloader_params.get('dataset_version', 'filteredv3')}"
    output_dir = os.path.join(config['output_dir_base'], output_dir_name)
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Experiment directory created: {output_dir}")

    # Get prediciton on the test set
    best_pred = (gt_test @ g2p_cov_model.best_posterior_mean_beta) + (covariates_test@g2p_cov_model.best_gamma) + g2p_cov_model.best_intercept

    logger.info(f"--- Saving Model Outputs to {output_dir} ---")
    utils.save_model_outputs(g2p_cov_model, trait_measurement_test, y_test_residual, best_pred, id_test, gene_list, trait, output_dir, config)

if __name__ == "__main__":
    cli()