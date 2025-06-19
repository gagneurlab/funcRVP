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

from ..utils import dataloader_old, utils
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
    best_model_arch: Optional[dict]=None,
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
    training_params['best_model_arch'] = best_model_arch
    for rep in range(n_repeats):

        logger.info(f"Training model for {trait} - Repetition {rep+1}/{n_repeats}")
        # Initialize model and send to device
        g2p_cov_model = g2p_bayes_model.G2P_Model(
            emb.shape[1], # embedding_dim
            covariates_train.shape[1],
            n_hidden=best_model_arch["n_hidden"],
            hiddem_dim=best_model_arch["hidden_dim"],
            last_layer_bias=best_model_arch["last_layer_bias"],
            nonlinearity=best_model_arch['nonlinearity'],
            y_var_init=best_model_arch['y_var_init'],
            base_var_init=best_model_arch['base_var_init'],
            base_var_const=0, 
            n_genes=emb.shape[0],
            device=device,
        ).to(device)

        g2p_cov_model = torch.compile(g2p_cov_model, dynamic=False, mode='max-autotune-no-cudagraphs') #Compile the model for better performance

        # Initialize weights and biases
        if training_params.get('wandb_logging', False):
            logger.info(f"Writing wandb logs to {training_params['wandb_dir']}")
            wandb.init(
                project=training_params['wandb_project_name'],
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
            logging=training_params['wandb_logging'],
            device=device,
        )

        logger.info("Model training finished.")
        if training_params['wandb_logging']:
            wandb.finish()

    return g2p_cov_model


@click.group()
def cli():
    pass

@cli.command()
@click.option("--trait", type=str, required=True, help="Trait to train model for")
@click.option("--config-path", type=click.Path(exists=True), required=True, help="Config file with all details")
@click.option("--trainval-sampling", type=float, required=False, help="How much should we sample from the training and validation set")
def trainer(
    trait: str,
    config_path: Path,
    trainval_sampling: Optional[float] = None,
):
    # --- Run Configuration Loading (load_data parameters) ---
    if not os.path.exists(config_path):
        logger.warning(f"Run configuration file not found: {config_path}. Using default parameters.") # Warning instead of error, defaults will be used.
        config = {} # Empty dict if run_config.yaml is missing
    else:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

    exp_name = config.get('experiment_name', 'default')

    # --- Load Run Parameters from Config ---
    dataloader_params = config.get('dataloader_params', None)
    if dataloader_params.get('embedding_path', None):
        embedding = dataloader_params.get('embedding_path', None).split("/")[-1].split(".")[0]
    if dataloader_params.get('genotype_path', None):
        genotype = dataloader_params.get('genotype_path', None).split("/")[-1].split(".")[0]

    (
        (gt_train, gt_val, gt_test),
        (y_train_residual, y_val_residual, y_test_residual),
        emb,
        gene_list,
        (id_train, id_val, id_test),
        (trait_measurement_train, trait_measurement_val, trait_measurement_test),
        (covariates_train, covariates_val, covariates_test),
    ) = dataloader_old.load_data(
        trait,
        embedding_path = dataloader_params.get('embedding_path', None), 
        genotype_path = dataloader_params.get('genotype_path', None),
        phenotype_dir= dataloader_params.get('phenotype_dir', None),
        covariates_path = dataloader_params.get('covariates_path', None),
        prs_path = dataloader_params.get('prs_path', None),
        train_individuals_path = dataloader_params.get('train_individuals_path', '-1'),
        test_split_size = dataloader_params.get('test_split_size', 0.25),
        val_split_size = dataloader_params.get('val_split_size', 0.1),
        trainval_sampling = trainval_sampling, #dataloader_params.get('trainval_sampling', 1),
        split_seed = dataloader_params.get('split_seed', 0),
        use_prs = dataloader_params.get('use_prs', True),
        normalize_covariates = dataloader_params.get('normalize_covariates', True),
        normalize_embedding=dataloader_params.get('normalize_embedding', False),
        shuffled_phenotype = dataloader_params.get('shuffled_phenotype', False),
        shuffled_embedding = dataloader_params.get('shuffled_embedding', False),
        random_embedding = dataloader_params.get('random_embedding', False),
        gene_subset = dataloader_params.get('gene_subset', None),
        dataset_version = dataloader_params.get('dataset_version', "filteredv3"),
    )

    logger.info(f"Using {gt_train.shape[0]} train samples, {gt_val.shape[0]} val samples, and {gt_test.shape[0]} test samples")

    # Read Training and hyperopt params from config
    training_params = config.get('training_params', None)
    training_params['exp_name'] = exp_name
    training_params['trait'] = trait
    training_params['embedding'] = embedding
    training_params['genotype'] = genotype

    hpopt_params = config.get('hpopt_params', None)
    # If no hyperopt, get best model architecture
    hpopt_params = None
    if hpopt_params is None:
        n_repeats = training_params.get('n_repeats', 1)

        # best_model_arch = utils.get_best_arch(
        #     trait=trait,
        #     study_version=training_params['old_model_version'], # Access model_params
        #     embedding_type=training_params['old_model_embedding'], # Access model_params
        #     storage_path=training_params['old_optuna_journal_log'] # Access paths through 'config'
        # )
        hparam_df = pd.read_csv(training_params.get('hparams_df', None), sep='\t')
        best_model_arch = hparam_df.query("trait==@trait").iloc[0].to_dict()

        best_model_arch["y_var_init"] = y_train_residual.var()
        
        # TODO! should we specify below numbers in the config file?
        # To get constant regularization
        if embedding == None:
            # best_n_hidden = -1
            best_model_arch["n_hidden"] = -1 # number of hidden layers in f(E)
            best_model_arch["base_var_init"] = 5e-3
        else:
            # best_n_hidden = best_model_arch["n_hidden"]  # number of hidden layers in f(E)
            best_model_arch["base_var_init"] = 5e-5

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
            best_model_arch=best_model_arch,
            training_params=training_params,
        )

    # TODO fix optuna suggetions for hyperparams
    else:
        # TODO Get optuna suggestions for hyperparams
        logger.info("Hyperparameter optimization is under development.")
        sys.exit(1)

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
    if config.get('output_dir_base', None):
        run_op_dir = os.path.join(config['output_dir_base'], genotype, embedding)
    else:
        logger.error("Output directory not defined in config. Exiting.")
        sys.exit(1)

    # output_dir_name = f"{exp_name}_{dataloader_params.get('dataset_version', 'filteredv3')}"
    output_dir_name = f"{exp_name}_{dataloader_params.get('dataset_version', 'filteredv3')}{'_sampling' if trainval_sampling else ''}{trainval_sampling if trainval_sampling else ''}"
    output_dir = os.path.join(run_op_dir, output_dir_name)
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Experiment directory created: {output_dir}")

    # --- Move model to cpu ---
    g2p_cov_model = g2p_cov_model.to('cpu')

    logger.info("Generating predictions on the test set...")
    with torch.no_grad(): # Disable gradient calculations
        best_model_test_pred = ((gt_test @ g2p_cov_model.best_posterior_mean_beta) + (covariates_test@g2p_cov_model.best_gamma) + g2p_cov_model.best_intercept)

        logger.info(f"--- Saving Model Outputs to {output_dir} ---")
        # utils.save_model_outputs(g2p_cov_model, trait_measurement_test, y_test_residual, best_model_test_pred, id_test, gene_list, trait, output_dir, config)
        utils.save_model_betas(g2p_cov_model, gene_list, trait, output_dir, config)
        if len(id_test) > 0:
            utils.save_model_predictions(trait_measurement_test, y_test_residual, best_model_test_pred, id_test, trait, output_dir, config)
        

if __name__ == "__main__":
    cli()