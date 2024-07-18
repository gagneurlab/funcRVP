import joblib
import os.path
import argparse

import pandas as pd
import numpy as np
import scipy
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression, ElasticNetCV
import matplotlib.pyplot as plt

import torch, gc
import plotnine as pn

import wandb
import optuna
import dataloader
import g2p_bayes_cov_skipcon as g2p_bayes_cov


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def get_best_arch(trait, study_version="v108cov_deepRVAT", embedding_type="omics_pops"):

    best_trial_params = {}
    
    study_name = f"study_{study_version}_{trait}_{embedding_type}"
    study = optuna.load_study(study_name=study_name, storage=storage)

    if study._is_multi_objective():
        best_trial = max(study.best_trials, key=lambda t: t.values[1])
    else:
        best_trial = study.best_trial

    best_trial_params = {
        "trial_number": best_trial.number,
        "r2": best_trial.values[0],
        "n_total_trials": len(study.trials_dataframe().query("state=='COMPLETE'")),
        "n_hidden": best_trial.params['n_hidden'], 
        "hidden_dim": best_trial.params['hidden_dim'], 
        "last_layer_bias": best_trial.params['last_layer_bias'], 
        "alpha_L1_fE": best_trial.params['alpha_L1_fE'],
        "normalize_embeddings": best_trial.params['normalize_embeddings'],
    }
        
    return best_trial_params


# Create ArgumentParser object
parser = argparse.ArgumentParser(description='model parameters')

# Add arguments
parser.add_argument('--trait', '-t', type=str, help='Trait')
parser.add_argument('--embedding', '-e', type=str, default='omics_pops', help='Embedding')
parser.add_argument('--genotype', '-g', type=str, default='deepRVAT', help='Genotype matrix type')
parser.add_argument('--version', '-v', type=str, help='Model version')
parser.add_argument('--shuffled_pheno', action='store_true', default=False, help='Shuffle phenotype')
parser.add_argument('--shuffled_emb', action='store_true', default=False, help='Shuffle embedding-gene map')
parser.add_argument('--random_seed', type=int, default=1234, help='Random seed to shuffle phenotype')

# Parse command-line arguments
args = parser.parse_args()

trait = args.trait
embedding_type = args.embedding
genotype = args.genotype
version = args.version
shuffled_phenotype = args.shuffled_pheno
shuffled_embedding = args.shuffled_emb
np_random_seed = args.random_seed

# shuffled_embedding = False

study_version = f"{version}{f'_shuffledpheno{np_random_seed}' if shuffled_phenotype else ''}{f'_shuffledemb{np_random_seed}' if shuffled_embedding else ''}_{genotype}"
study_name = f"study_{study_version}_{embedding_type}"

    
only_genebass_genes = False #restrict everything to genebass genes
extend_covariates = True #extend covariates with age**2, age:sex, age**2:sex 
normalize = True #z-score normalize covariates
dataset_version = "filteredv3" #dataset to use (filteredv2 has removed non-Caucasians)

(gt_train, gt_val, gt_test), (y_train_residual, y_val_residual, y_test_residual), emb, gene_list, (id_train, id_val, id_test), (trait_measurement_train, trait_measurement_val, trait_measurement_test), (covariates_train, covariates_val, covariates_test) = dataloader.load_data(trait, embedding_type, add_loeuf=False, add_alpha_mis=False, add_gwas_hits=False, only_genebass_genes=only_genebass_genes, extend_covariates=extend_covariates, normalize=normalize, version=dataset_version, genotype=genotype)

if embedding_type is None:
    emb = np.zeros(len(gene_list))
    
if shuffled_phenotype:
    # np_random_seed = 1225
    np.random.seed(np_random_seed)
    trait_measurement_train = np.random.permutation(trait_measurement_train)
    trait_measurement_val = np.random.permutation(trait_measurement_val)

if shuffled_embedding:
    np.random.seed(np_random_seed)
    emb = np.random.permutation(emb)


storage = optuna.storages.JournalStorage(
    optuna.storages.JournalFileStorage(f"/s/project/uk_biobank/processed/g2p/optuna/journal.log"),
)
    
# best_model_arch = get_best_arch(trait=trait, study_version="v108cov_deepRVAT", embedding_type=embedding_type)
best_model_arch = get_best_arch(trait=trait, study_version="v108cov_deepRVAT", embedding_type="omics_pops")


args = {
    "trait": trait,
    "genotype": genotype,
    "study_name": study_name,
    "study_version": study_version,
    "dataset_version": dataset_version,
    "extend_covariates": extend_covariates,
    "np_seed": np_random_seed,
    "embedding_type": embedding_type,
    "shuffled_embedding": shuffled_embedding,
    "shuffled_phenotype": shuffled_phenotype,
    "add_loeuf": False,
    "add_alpha_mis": False,
    "add_gwas_hits": False,
    "skip_con": True,
    "normalize_embeddings": best_model_arch['normalize_embeddings'], #z-score normalize embeddings
    "embedding_dim": emb.shape[1], #dimension of embedding
    "n_genes": len(gene_list), #number of genes
    "n_samples": len(y_train_residual), #number of training samples
    "learning_rate": 0.001, #lr for Adam optimizer 
    "epochs": 25, #Number of epochs to train
    "nonlinearity": "softplus", #choose from ['relu', 'silu', 'elu', 'gelu', 'softplus']
    "weight_decay": 0, #L2 penalty on trainable parameters
    "batch_size": 16*1024,
    "n_hidden": best_model_arch['n_hidden'], #number of hidden layers in f(E)
    "hidden_dim": best_model_arch['hidden_dim'], #number of neurons in each hidden layer of f(E)
    "last_layer_bias": best_model_arch['last_layer_bias'], #initialization of bias term of last layer in f(E)
    "early_stopping": False, #If True, stops training if val r2 is lower than covariate r2 for >5 epochs
    "y_var_init": y_train_residual.var(), #initialization of variance of y
    "alpha_L1_fE": best_model_arch['alpha_L1_fE'], #L1 penalty weighting of f(E)
    "base_var_init": 5e-5, #Initialization of additional bias term to variance (to be added to f(E))
    "base_var_const": 0, #Constant variance value to be added to f(E)
}

g2p_cov_model = g2p_bayes_cov.G2P_Model(args["embedding_dim"], covariates_train.shape[1], n_hidden=args["n_hidden"], hiddem_dim=args["hidden_dim"], last_layer_bias=args["last_layer_bias"], nonlinearity=args["nonlinearity"], y_var_init=args["y_var_init"], base_var_init=args["base_var_init"], base_var_const=args["base_var_const"], n_genes=args["n_genes"], device=device).to(device)

# Initialize weights and biases
wandb.init(project="g2p_cov_hyperopt", name=f"{trait}_{embedding_type}", config=args, dir="/s/project/uk_biobank/processed/g2p/wandb/", settings=wandb.Settings(_service_wait=300))
wandb.watch(g2p_cov_model)

# print("Number of trainable parameters in model:", sum(p.numel() for p in g2p_cov_model.parameters() if p.requires_grad))

g2p_cov_model.fit_model(gt_train, covariates_train, trait_measurement_train, 
                        emb, 
                        args, 
                        G_val=gt_val, C_val=covariates_val, y_val=trait_measurement_val,
                        gene_list=gene_list,
                        fast=False, #Compute posterior distribution after every epoch. Necessary for r2. 
                        logging=True,
                        device=device
)

# Finish weights and biases log
wandb.finish()

output_dir = f"/s/project/uk_biobank/processed/g2p/modelling/hyperopt/fixed_arch/{study_name}"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Log the parameters of the model 
mean_betas_df = pd.DataFrame({"best_r2_mean_beta": g2p_cov_model.best_r2_trainval_mean_beta, 
                              "best_r2_var_beta": g2p_cov_model.best_r2_trainval_var_beta, 
                              "best_r2_intercept": g2p_cov_model.best_r2_intercept, 
                              "best_r2_base_var": g2p_cov_model.best_r2_base_var,
                              "best_r2_last_layer_bias": g2p_cov_model.best_r2_last_layer_bias,
                              "best_r2_epoch": g2p_cov_model.best_r2_epoch,
                              "best_r2_fE": g2p_cov_model.best_r2_fE.flatten(),
                              "best_loss_mean_beta": g2p_cov_model.best_loss_trainval_mean_beta, 
                              "best_loss_var_beta": g2p_cov_model.best_loss_trainval_var_beta, 
                              "best_loss_intercept": g2p_cov_model.best_loss_intercept, 
                              "best_loss_base_var": g2p_cov_model.best_loss_base_var,
                              "best_loss_last_layer_bias": g2p_cov_model.best_loss_last_layer_bias,
                              "best_loss_epoch": g2p_cov_model.best_loss_epoch,
                              "best_loss_fE": g2p_cov_model.best_loss_fE.flatten()}, index=gene_list)

mean_betas_df["std_err"] = np.sqrt(mean_betas_df["best_r2_var_beta"])
mean_betas_df["mean_beta"] = mean_betas_df["best_r2_mean_beta"]

mean_betas_df["neglog_pval"] = -np.minimum(np.log(2) + scipy.stats.norm.logcdf(0, mean_betas_df["mean_beta"], mean_betas_df["std_err"]), np.log(2) + scipy.stats.norm.logsf(0, mean_betas_df["mean_beta"], mean_betas_df["std_err"]))
mean_betas_df["pd"] = np.maximum(scipy.stats.norm.cdf(0, mean_betas_df["mean_beta"], mean_betas_df["std_err"]), scipy.stats.norm.sf(0, mean_betas_df["mean_beta"], mean_betas_df["std_err"]))
mean_betas_df["significant"] = mean_betas_df["pd"] > 0.999
mean_betas_df["trait"] = trait
mean_betas_df["base_var_const"] = args["base_var_const"]
mean_betas_df["version"] = study_version.rstrip(',')
mean_betas_df = mean_betas_df.reset_index().rename(columns={"index": "gene_id"})
mean_betas_df.to_parquet(f"{output_dir}/{trait}_mean_betas.pq")   # Write to file

# Log the model predictions on the test data
bayes_pred_df = pd.DataFrame({f"{trait}_measurement": trait_measurement_test, 
                              "common_residual": y_test_residual,
                              "best_r2_pred": (gt_test @ g2p_cov_model.best_r2_trainval_mean_beta) + (covariates_test@g2p_cov_model.best_r2_gamma) + g2p_cov_model.best_r2_intercept, 
                              "best_loss_pred": (gt_test @ g2p_cov_model.best_loss_trainval_mean_beta) + (covariates_test@g2p_cov_model.best_loss_gamma) + g2p_cov_model.best_loss_intercept}, index=id_test) 
bayes_pred_df["trait"] = trait
bayes_pred_df["model"] = embedding_type + f"_bayesian_{study_version.rstrip(',')}"
bayes_pred_df["version"] = study_version.rstrip(',')
bayes_pred_df["genotype"] = "deepRVAT" if "deepRVAT" in study_version else "pLoF"
bayes_pred_df.to_parquet(f"{output_dir}/{trait}_bayes_pred.pq")    # Write to file