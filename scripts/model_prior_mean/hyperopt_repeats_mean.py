import sys
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.metrics import r2_score
import optuna
import joblib
import os.path
import g2p_bayes_cov_mean as g2p_bayes_cov
import wandb
import torch
from pathlib import Path
import dataloader


EMBEDDINGS = {
    "string": "/s/project/gene_embedding/embedding/final/STRING_d128.tsv",
    "string_exp": "/s/project/gene_embedding/embedding/final/STRING_EXP_d128.tsv",
    "omics": "/s/project/gene_embedding/embedding/final/dtf_gtex_depMap_portT5_lunar-snowflake-239_3787637_embedding.tsv",
    "pops": "/s/project/gene_embedding/embedding/final/pops_mat_d256.tsv",
    "pops_exp": "/s/project/gene_embedding/embedding/final/pops_mat_exp_d256.tsv",
    "omics_pops": "/s/project/gene_embedding/embedding/combination/pops_mat_pca256_omics.tsv",
    "omics_pops_exp": "/s/project/gene_embedding/embedding/combination/pops_mat_exp256_omics.tsv"
}

def main():
    if not torch.cuda.is_available():
        print("Cuda not available")
        return
    
    # Retrieve command-line arguments
    arguments = sys.argv[1:]  # Exclude the script name

    # Process and use the command-line arguments
    trait = arguments[0]
    embedding_type = arguments[1]
    genotype = arguments[2]
    version = arguments[3]
    
    np_random_seed = 1234
    np.random.seed(np_random_seed)
    
    print(trait, embedding_type, genotype)
    
    add_loeuf = False
    dataset_version = "filteredv3"
    extend_covariates = True
    (gt_train, gt_val, gt_test), (y_train, y_val, y_test), emb, gene_list, (id_train, id_val, id_test), (trait_measurement_train, trait_measurement_val, trait_measurement_test), (covariates_train, covariates_val, covariates_test) = dataloader.load_data(trait, embedding_type, add_loeuf=add_loeuf, only_genebass_genes=False, extend_covariates=extend_covariates, normalize=True, version=dataset_version, genotype=genotype)
    
    n_genes = emb.shape[0]
    n_samples = gt_train.shape[0]
    print(n_genes, n_samples)
    
    # study_version = f"v97cov_{genotype}"
    study_version = f"{version}_{genotype}"
    study_basedir = "/s/project/uk_biobank/processed/g2p/mean_model/optuna/"
    study_name = f"study_{study_version}_{trait}_{embedding_type}"
    output_dir = f"/s/project/uk_biobank/processed/g2p/mean_model/hyperopt/{study_name}"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    storage = optuna.storages.JournalStorage(
       optuna.storages.JournalFileStorage(f"{study_basedir}/journal.log"),
    )
    sampler = optuna.samplers.TPESampler(multivariate = True, n_startup_trials=10)
    study = optuna.create_study(directions=["maximize"], study_name=study_name, storage=storage, load_if_exists=True, sampler = sampler)
    
    def objective(trial):
        shuffled_embedding = False
        normalize_embeddings = False
        args = {
            "trait": trait,
            "genotype": genotype,
            "study_name": study_name,
            "study_version": study_version,
            "dataset_version": dataset_version,
            "extend_covariates": extend_covariates,
            "trial_id": trial.number,
            "np_seed": np_random_seed,
            "embedding_type": embedding_type,
            "shuffled_embedding": shuffled_embedding,
            "normalize_embeddings": normalize_embeddings,
            "embedding_dim": emb.shape[1],
            "add_loeuf": add_loeuf,
            "n_genes": n_genes,
            "n_samples": n_samples,
            "learning_rate": 0.001,#trial.suggest_float("learning_rate", 0.00001, 0.01, log=True),
            "epochs": 50,
            "nonlinearity": trial.suggest_categorical("nonlinearity", ["softplus"]),
            "weight_decay": 0,#trial.suggest_float("weight_decay", 1e-10, 1e-2, log=True),
            "batch_size": 1024*8, #trial.suggest_int("batch_size", 1024*16, 1024*16),
            "n_hidden": trial.suggest_int("n_hidden", 1, 3),
            "hidden_dim": trial.suggest_int("hidden_dim", 2, 16),
            "last_layer_bias": trial.suggest_float("last_layer_bias", -3, 3),
            "early_stopping": True,
            "y_var_init": y_train.var(),
            "beta_var_init": trial.suggest_float("base_var", 1e-5, 1e-1, log=True),
            "alpha_L1_fE": trial.suggest_float("alpha_L1_fE", 1e-10, 10, log=True),
            "base_var": 0,#trial.suggest_float("base_var", 1e-10, 1, log=True)
        }
        
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        repetitions = 1
        best_r2_list = []
        final_r2_list = []
        best_loss_list = []
        for repetition in range(repetitions):
        
            # Generate the model.
            g2p_model = g2p_bayes_cov.G2P_Model(args["embedding_dim"], covariates_train.shape[1], n_hidden=args["n_hidden"], hiddem_dim=args["hidden_dim"], last_layer_bias=args["last_layer_bias"], y_var_init=args["y_var_init"], nonlinearity=args["nonlinearity"], beta_var_init=args['beta_var_init'], base_var=args["base_var"], device=device).to(device)

            wandb.init(project="g2p_cov_mean_hyperopt", name=f"{trait}_{embedding_type}_{trial.number}_{repetition}", config=args, dir="/s/project/uk_biobank/processed/g2p/mean_model/wandb/", settings=wandb.Settings(_service_wait=300))
            wandb.run.summary["repetition"] = repetition
            
            run_id = wandb.run.id
            wandb.watch(g2p_model)

            g2p_model.fit_model(gt_train, covariates_train, trait_measurement_train, emb, args, G_val=gt_val, C_val=covariates_val, y_val=trait_measurement_val, gene_list=gene_list, logging=True, trial=None, device=device)

            wandb.finish()

            val_r2 = r2_score(y_val, (gt_val@g2p_model.mean_beta) + (covariates_val@g2p_model.gamma) + g2p_model.intercept)
            best_r2 = g2p_model.best_r2
            val_loss = g2p_model.best_loss

            pd.DataFrame({"best_r2_mean_beta": g2p_model.best_r2_trainval_mean_beta, "best_r2_var_beta": g2p_model.best_r2_trainval_var_beta, "best_r2_intercept": g2p_model.best_r2_intercept, "best_r2_epoch": g2p_model.best_r2_epoch, "best_loss_mean_beta": g2p_model.best_loss_trainval_mean_beta, "best_loss_var_beta": g2p_model.best_loss_trainval_var_beta, "best_loss_intercept": g2p_model.best_loss_intercept, "best_loss_epoch": g2p_model.best_loss_epoch}, index=gene_list).to_parquet(f"{output_dir}/{trial.number}_{repetition}_mean_betas.pq")
            pd.DataFrame({f"{trait}_measurement": trait_measurement_test, "common_residual": y_test, "best_r2_pred": (gt_test @ g2p_model.best_r2_trainval_mean_beta) + (covariates_test@g2p_model.best_r2_gamma) + g2p_model.best_r2_intercept, "best_loss_pred": (gt_test @ g2p_model.best_loss_trainval_mean_beta) + (covariates_test@g2p_model.best_loss_gamma) + g2p_model.best_loss_intercept}, index=id_test).to_parquet(f"{output_dir}/{trial.number}_{repetition}_bayes_pred.pq")
            pd.DataFrame({"var": g2p_model.var}, index=["var"]).to_parquet(f"{output_dir}/{trial.number}_{repetition}_var.pq")
            pd.DataFrame({"run_id": run_id}, index=["run_id"]).to_parquet(f"{output_dir}/{trial.number}_{repetition}_run_id.pq")
            
            pd.DataFrame({"prior_var": np.squeeze(g2p_model.prior_var)}, index=gene_list).to_parquet(f"{output_dir}/{trial.number}_{repetition}_prior_var.pq")
            pd.DataFrame({"prior_mean": np.squeeze(g2p_model.prior_mean)}, index=gene_list).to_parquet(f"{output_dir}/{trial.number}_{repetition}_prior_mean.pq")

            torch.save(g2p_model.state_dict(), f"{output_dir}/{trial.number}_{repetition}_model_state.m")
            
            best_r2_list.append(best_r2)
            best_loss_list.append(val_loss)
            final_r2_list.append(val_r2)

        trial.set_user_attr("best_r2_list", best_r2_list)
        trial.set_user_attr("best_loss_list", best_loss_list)
        trial.set_user_attr("final_r2_list", final_r2_list)

        value = sorted(best_r2_list)[1] if len(best_r2_list)>1 else best_r2_list[0]
        return value

    def get_num_trials_done(study):
        trials_df = study.trials_dataframe()
        if trials_df.empty:
            return 0
        else:
            return len(trials_df.query("state == 'COMPLETE'"))
    
    while get_num_trials_done(study) < 100:
        study.optimize(objective, n_trials=1)

if __name__ == "__main__":
    main()
