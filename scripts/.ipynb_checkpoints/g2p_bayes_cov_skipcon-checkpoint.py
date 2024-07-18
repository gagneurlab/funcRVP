import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lr_scheduler
from torch.distributions.multivariate_normal import MultivariateNormal
from sklearn.metrics import r2_score
import scipy.sparse as sp
import wandb
import time
import warnings
import sys
import yaml
import copy
import statsmodels.formula.api as smf
import sys
import optuna
import dataloader

from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LassoCV, LinearRegression
from sklearn.metrics import PrecisionRecallDisplay, precision_recall_curve, average_precision_score
from sklearn.decomposition import PCA

# Ignoring sparse warning
warnings.filterwarnings("ignore")
torch.set_float32_matmul_precision('high')
torch.backends.cuda.matmul.allow_tf32 = True

class EarlyStopper:
    def __init__(self, patience=5, baseline=-1):
        self.patience = patience
        self.baseline = baseline
        self.counter = 0

    def early_stop(self, validation_r2):
        if validation_r2 >= self.baseline:
            self.counter = 0
        elif validation_r2 < self.baseline:
            self.counter += 1
            if self.counter >= self.patience:
                return True
        return False

class ConstantModule(nn.Module):
    def __init__(self, constant_value):
        super(ConstantModule, self).__init__()
        # Create a parameter with the specified constant value
        self.constant = nn.Parameter(torch.tensor(float(constant_value)), requires_grad=True)

    def forward(self, x):
        return self.constant.expand(x.shape[0]).unsqueeze(1)

# Class that models variance of beta as a function of the embedding
class VarPredModel(nn.Module):

    def __init__(self, emb_dim, n_cov, n_hidden=0, hiddem_dim=None, last_layer_bias = None, nonlinearity="relu", final_nonlinearity="softplus", y_var_init=1e-3, gene_var_init=None, base_var_init=5e-5, base_var_const=0, batch_norm="linear", n_genes=None, skip_con=False, device=None):
        super(VarPredModel, self).__init__()
        
        if nonlinearity == 'relu':
            self.activation = nn.ReLU()
        elif nonlinearity == 'silu':
            self.activation = nn.SiLU()
        elif nonlinearity == 'gelu':
            self.activation = nn.GELU()
        elif nonlinearity == 'softplus':
            self.activation = nn.Softplus()
        elif nonlinearity == 'elu':
            self.activation = nn.ELU()
        elif nonlinearity == 'celu':
            self.activation = nn.CELU()
        else:
            self.activation = nn.Identity()
        
        self.layers = nn.ModuleList()
        if n_hidden>=0:
            if n_hidden>0:
                self.layers.append(nn.Linear(emb_dim, hiddem_dim or emb_dim))
                if batch_norm == "linear":
                    self.layers.append(nn.BatchNorm1d(hiddem_dim or emb_dim))
                for _ in range(n_hidden - 1):
                    self.layers.append(self.activation)
                    if batch_norm == "activation":
                        self.layers.append(nn.BatchNorm1d(hiddem_dim or emb_dim))
                    self.layers.append(nn.Linear(hiddem_dim or emb_dim, hiddem_dim or emb_dim))
                    if batch_norm == "linear":
                        self.layers.append(nn.BatchNorm1d(hiddem_dim or emb_dim))
                self.layers.append(self.activation)
                if batch_norm == "activation":
                        self.layers.append(nn.BatchNorm1d(hiddem_dim or emb_dim))
                self.layers.append(nn.Linear(hiddem_dim or emb_dim, 1))
            else:
                self.layers.append(nn.Linear(emb_dim, 1))
            if last_layer_bias:
                self.layers[-1].bias.data.fill_(last_layer_bias)
            
        else:
            self.layers.append(ConstantModule(last_layer_bias if last_layer_bias else 0))
            
        if final_nonlinearity=="softplus":
            self.layers.append(nn.Softplus())
        elif final_nonlinearity=="sigmoid":
            self.layers.append(nn.Sigmoid())
        else:
            self.layers.append(nn.Softplus())
            
       
        # Linear Residual connection
        self.skip_con = skip_con
        if skip_con:
            self.shortcut = nn.Sequential(
                nn.Linear(emb_dim, 1),
                nn.Softplus()
            )
            self.shortcut[0].bias.data.fill_(last_layer_bias)
        
        # Initialize variance of y|b
        self.log_var = nn.Parameter(torch.FloatTensor([np.log(y_var_init)]))
        
        # Initialize intercept parameter
        self.intercept = nn.Parameter(torch.FloatTensor([0]))
        
        # Initialize covariate effects
        self.gamma = nn.Parameter(torch.zeros(n_cov).to(device))
        
        if n_genes and gene_var_init:
            self.gene_var = nn.Parameter(gene_var_init*torch.ones(n_genes, 1).to(device))
        else:
            self.gene_var = None
        
        # If you want to add some additional variance to f(E)
        # self.base_var = base_var_init
        self.log_base_var = nn.Parameter(torch.FloatTensor([np.log(base_var_init)]).to(device))
        self.base_var_const = base_var_const
        
    def forward(self, emb):
        if self.skip_con:
            emb_copy = emb.detach().clone()
        
        for layer in self.layers:
            emb = layer(emb)
        
        if self.skip_con:
            emb = emb + self.shortcut(emb_copy)
            
        if not (self.gene_var is None):
            emb = emb + torch.exp(self.gene_var)
            
        emb = emb + torch.exp(self.log_base_var) + self.base_var_const
        
        return emb, self.log_var, self.intercept, self.gamma, torch.exp(self.log_base_var)


# Class that fits the genotype to phenotype model
class G2P_Model(nn.Module):
    
    def __init__(self, emb_dim, n_cov, n_hidden=0, hiddem_dim=None, last_layer_bias=None, nonlinearity="relu", final_nonlinearity="softplus", y_var_init=1e-3, gene_var_init=None, base_var_init=1e-7, base_var_const=0, batch_norm="linear", n_genes=None, skip_con=False, device=None):
        super(G2P_Model, self).__init__()
        
        # Initialize the f(E) model
        self.var_pred_model = VarPredModel(emb_dim, n_cov, n_hidden=n_hidden, hiddem_dim=hiddem_dim, last_layer_bias=last_layer_bias, nonlinearity=nonlinearity, final_nonlinearity=final_nonlinearity, y_var_init=y_var_init, gene_var_init=gene_var_init, base_var_init=base_var_init, base_var_const=base_var_const, batch_norm=batch_norm, n_genes=n_genes, skip_con=skip_con, device=None)
        
        self.total_epochs_trained = 0
        
        # Initialize/Record mean and variance of beta
        self.mean_beta = None
        self.var_beta = None
        
        # Initialize metrics
        self.best_r2 = -np.inf
        self.best_r2_mean_beta = None
        self.best_r2_intercept = None
        self.best_r2_var_beta = None
        self.best_r2_epoch = None
        self.best_r2_base_var = None
        
        self.best_loss = np.inf
        self.best_loss_mean_beta = None
        self.best_loss_intercept = None
        self.best_loss_var_beta = None
        self.best_loss_epoch = None
        self.best_loss_base_var = None
        
        self.epoch_loss_list = []
        self.loss_list = []
        self.val_list = []
        self.train_r2_list = []
        self.val_r2_list = []
        self.gE_by_epoch = []
        self.var_by_epoch = []
        self.intercept_by_epoch = []
        self.gamma_by_epoch = []
        self.base_var_by_epoch = []
        
        
    # Compute likelihood of y
    def forward(self, G, emb, C, device=None):
        gE, var, b, gamma, base_var = self.var_pred_model(emb)
        cov = (G * gE.squeeze(1)) @ G.transpose(1,0) + torch.diag(torch.exp(var).expand(G.shape[0]))
        pred = MultivariateNormal(loc = (C @ gamma) + b.expand(G.shape[0]), scale_tril = torch.linalg.cholesky(cov)) #torch.linalg???
        return pred
    
    # Returns variance of prior
    def get_prior_var(self, emb):
        gE, _, _, _, _ = self.var_pred_model(emb)
        self.prior_var = gE.detach().cpu().numpy()
        return self.prior_var
    
    # Returns variance of y
    def get_var(self, emb):
        _, var, _, _, _ = self.var_pred_model(emb)
        self.var = var.detach().cpu().numpy()
        return self.var
    
    # Returns intercept (y)
    def get_intercept(self, emb):
        _, _, intercept, _, _ = self.var_pred_model(emb)
        self.intercept = intercept.detach().cpu().numpy()[0]
        return self.intercept
    
    # Returns effect of covariates    
    def get_gamma(self, emb):
        _, _, _, gamma, _ = self.var_pred_model(emb)
        self.gamma = gamma.detach().cpu().numpy()
        return self.gamma

    def get_base_var(self, emb):
        _, _, _, _, base_var = self.var_pred_model(emb)
        self.base_var = base_var.detach().cpu().numpy()
        return self.base_var
        
    # Returns mean of posterior
    def get_posterior_mean(self, G, emb, C, y, device=None):
        with torch.no_grad():
            n_genes = G.shape[1]
            gE, var, b, gamma, base_var = self.var_pred_model(emb)
            sigma_inv = torch.diag((1/gE).squeeze()) + (1/torch.exp(var))[0] * (torch.transpose(G, 0, 1) @ G)
            sigma = torch.linalg.inv(sigma_inv)
            mean_beta = (sigma @ ((1/torch.exp(var))[0] * torch.transpose(G, 0, 1) @ (y-((C @ gamma)+b)))).detach().cpu().numpy()
        self.mean_beta = mean_beta
        self.var_beta = np.diag(sigma.detach().cpu().numpy())
        return mean_beta
    
    # Compute and return mean and variance of posterior (Redundant - above function)
    def recompute_posterior(self, G, gE, var, b, gamma, C, y, device=None):
        with torch.no_grad():
            sigma_inv = torch.diag((1/gE).squeeze()) + (1/torch.exp(var))[0] * (torch.transpose(G, 0, 1) @ G)
            sigma = torch.linalg.inv(sigma_inv)
            mean_beta = (sigma @ ((1/torch.exp(var))[0] * torch.transpose(G, 0, 1) @ (y-((C @ gamma)+b)))).detach().cpu().numpy()
            var_beta = np.diag(sigma.detach().cpu().numpy())
        return mean_beta, var_beta
    
    # Phenotype prediction using posterior beta and covariates
    def predict(self, G, C, mean_beta=None, gamma=None):
        if mean_beta and gamma:
            return (G @ mean_beta) + (C @ gamma) + self.intercept
        else:
            return (G @ self.mean_beta) + (C @ self.gamma) + self.intercept
    
    # Function that runs training loop and logs parameters and metrics
    def fit_model(self, G, C, y, emb, args, G_val=None, C_val=None, y_val=None, gene_list=None, logging=False, fast=False, trial=None, device=None):
        
        # If you want to normalize embeddings along each dimension
        if args["normalize_embeddings"]:
            emb_torch = torch.from_numpy((emb / np.linalg.norm(emb, axis=0)).astype("float32")).to(device)
        else:
            emb_torch = torch.from_numpy(emb.astype("float32")).to(device)
        
        #if "batch_size_schedule" not in args:
        #    args["batch_size_schedule"] = {0: args["batch_size"]}
        
        if "learning_rate_schedule" not in args:
            args["learning_rate_schedule"] = {0: args["learning_rate"]}
        if "early_stopping" not in args:
            args["early_stopping"] = False
        
        # Fit Linear model on the covariates
        lm = LinearRegression().fit(C, y)
        
        # Get R^2 of the linear model on the validation set
        if (G_val is not None) and (C_val is not None) and (y_val is not None):
            common_r2 = r2_score(y_val, lm.predict(C_val))
            print("Covariates r2 on val:", common_r2)
        
        # Initialize covariate effects in FuncRVAT = effect of covariates from the LM
        self.var_pred_model.gamma.data = torch.from_numpy(lm.coef_.astype("float32")).to(device)

        # Initialize intercept in FuncRVAT = effect of covariates from the LM
        self.var_pred_model.intercept.data = torch.FloatTensor([lm.intercept_]).to(device)
        
        # Send data to GPU and make it a Tensor Dataset
        G_torch = torch.from_numpy(G.astype("float32")).to(device)
        C_torch = torch.from_numpy(C.astype("float32")).to(device)
        y_torch = torch.from_numpy(y.astype("float32")).to(device)
        dataset_train = torch.utils.data.TensorDataset(G_torch, C_torch, y_torch) 
        n_train_samples = len(dataset_train)
        
        if (G_val is not None) and (C_val is not None) and (y_val is not None):
            G_val_torch = torch.from_numpy(G_val.astype("float32")).to(device)
            C_val_torch = torch.from_numpy(C_val.astype("float32")).to(device)
            y_val_torch = torch.from_numpy(y_val.astype("float32")).to(device)
            dataset_val = torch.utils.data.TensorDataset(G_val_torch, C_val_torch, y_val_torch)
            n_val_samples = len(dataset_val)
        else:
            n_val_samples = 0
        
        # Initialize the Adam optimizer with learning rate and weight decay as input arguments from the user
        optimizer = torch.optim.Adam(self.parameters(), lr=args["learning_rate_schedule"][0], weight_decay=args["weight_decay"])
        batch_size = args["batch_size"]
        
        if logging:
            wandb.run.summary["num_parameters"] = sum(p.numel() for p in self.parameters() if p.requires_grad)
            wandb.run.summary["covariate_r2"] = common_r2
        
        # PyTorch speed-up trick
        scaler = torch.cuda.amp.GradScaler()
        
        # Instantiate early stopper
        if args["early_stopping"]:
            early_stopper = EarlyStopper(patience=5, baseline=common_r2)
        
        # START training
        for epoch in range(0, args["epochs"]):
            self.train()
            #if epoch in args["batch_size_schedule"]:
            #    optimizer = torch.optim.Adam(self.parameters(), lr=args["learning_rate"], weight_decay=args["weight_decay"])
            #    scheduler = lr_scheduler.LinearLR(optimizer, start_factor=1.0, end_factor=0.1, total_iters=10)
            #    batch_size = args["batch_size_schedule"][epoch]
            
            if epoch in args["learning_rate_schedule"]:
                optimizer = torch.optim.Adam(self.parameters(), lr=args["learning_rate_schedule"][epoch], weight_decay=args["weight_decay"])
                #print("New lr:", optimizer.param_groups[0]['lr'])
            
            # Record time
            epoch_start_time = time.time()
            epoch_loss = 0
            
            # Permute train samples (we don't want same batches in every epoch)
            permutation = torch.randperm(n_train_samples)
            
            # Loop through all batches
            for i in range(0, n_train_samples, batch_size):
                # Set gradient to zero
                optimizer.zero_grad()
                
                # Get indices to create batch
                indices = permutation[i:i+batch_size]
                G_batch, C_batch, y_batch = dataset_train[indices]
                
                # PyTorch speed-up trick
                with torch.autocast(device_type='cuda', dtype=torch.float16):
                    # Compute -log(likelihood). Here "pred" is an object of the class MultivariateNormal
                    pred = self.forward(G_batch, emb_torch, C_batch, device)
                    loss = -pred.log_prob(y_batch).mean()
                    
                    # Add extra regularization to the loss
                    if "alpha_L1_fE" in args:
                        loss += args["alpha_L1_fE"] * torch.norm(self.var_pred_model(emb_torch)[0], p=1)
                    if "alpha_L2_fE" in args:
                        loss += args["alpha_L2_fE"] * torch.norm(self.var_pred_model(emb_torch)[0], p=2)
                    if "alpha_L1_gene_var" in args and (not (self.var_pred_model.gene_var is None)):
                        loss += args["alpha_L1_gene_var"] * torch.norm(self.var_pred_model.gene_var, p=1)
                
                # PyTorch speed-up trick
                #loss.backward()
                scaler.scale(loss).backward()
                #optimizer.step()
                scaler.step(optimizer)
                scaler.update()
                
                # Record epoch loss
                epoch_loss += loss.detach().item() * len(y_batch)
                self.loss_list.append(loss.detach().item()/len(indices))
            
            self.epoch_loss_list.append(epoch_loss/n_train_samples)

            # Initiate evaluation mode (no gradients)
            self.eval()
            val_loss = 0
            
            # Loop over all batches
            for i in range(0, n_val_samples, batch_size):

                G_batch, C_batch, y_batch = dataset_val[i:i+batch_size]
                
                # Compute prediction and -log(likelihood) on validation data
                with torch.no_grad():
                    with torch.autocast(device_type='cuda', dtype=torch.float16):
                        pred = self.forward(G_batch, emb_torch, C_batch, device)
                        val_loss += -pred.log_prob(y_batch).mean().detach().item() * len(y_batch)
            
            # Record model parameters
            self.gE_by_epoch.append(self.get_prior_var(emb_torch))
            self.var_by_epoch.append(self.get_var(emb_torch))
            self.intercept_by_epoch.append(self.get_intercept(emb_torch))
            self.gamma_by_epoch.append(self.get_gamma(emb_torch))
            self.base_var_by_epoch.append(self.get_base_var(emb_torch))
            
            # Compute posterior to obtain R^2 on the validation set 
            if (not fast) or (epoch+1==args["epochs"]):
                mean_beta = self.get_posterior_mean(G_torch, emb_torch, C_torch, y_torch, device)
                if logging:
                    if self.total_epochs_trained==0:
                        mean_beta_df = pd.DataFrame({0: mean_beta}, index=gene_list)
                        mean_beta_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_mean_beta.pq")
                        var_beta_df = pd.DataFrame({0: self.var_beta}, index=gene_list)
                        var_beta_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_var_beta.pq")
                    else:
                        mean_beta_df = pd.read_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_mean_beta.pq")
                        mean_beta_df[epoch] = mean_beta
                        mean_beta_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_mean_beta.pq")
                        
                        var_beta_df = pd.read_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_var_beta.pq")
                        var_beta_df[epoch] = self.var_beta
                        var_beta_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_var_beta.pq")
                        
                    gE_df = pd.DataFrame({e: np.squeeze(self.gE_by_epoch[e]) for e in range(len(self.gE_by_epoch))}, index=gene_list)
                    gE_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_gE.pq")
                    var_df = pd.DataFrame({e: self.var_by_epoch[e] for e in range(len(self.var_by_epoch))}, index=["var"])
                    var_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_var.pq")
                    intercept_df = pd.DataFrame({e: self.intercept_by_epoch[e] for e in range(len(self.intercept_by_epoch))}, index=["intercept"])
                    intercept_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_intercept.pq")
                    gamma_df = pd.DataFrame({e: np.squeeze(self.gamma_by_epoch[e]) for e in range(len(self.gamma_by_epoch))})
                    gamma_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_gamma.pq")
                    # base_var_df = pd.DataFrame({e: np.squeeze(self.base_var_by_epoch[e]) for e in range(len(self.base_var_by_epoch))})
                    # base_var_df.to_parquet(f"/s/project/uk_biobank/processed/g2p/wandb/{wandb.run.id}_base_var.pq")
            
            if n_val_samples:
                self.val_list.append(val_loss/n_val_samples)
                if self.val_list[-1] < self.best_loss:
                    self.best_loss = self.val_list[-1]
                    self.best_loss_mean_beta = self.mean_beta
                    self.best_loss_intercept = self.intercept
                    self.best_loss_gamma = self.gamma
                    self.best_loss_var_beta = self.var_beta
                    self.best_loss_base_var = self.base_var
                    self.best_loss_last_layer_bias = self.var_pred_model.layers[-2].constant.item() if hasattr(self.var_pred_model.layers[-2], "constant") else self.var_pred_model.layers[-2].bias.item()
                    self.best_loss_epoch = self.total_epochs_trained
                    self.best_loss_fE = self.prior_var
            
            # Compute R^2 on train and validation set
            if (not fast) or (epoch+1==args["epochs"]):
                train_r2 = r2_score(y, (G@mean_beta)+(C@self.gamma)+self.intercept)
                self.train_r2_list.append(train_r2)
                if n_val_samples:
                    val_r2 = r2_score(y_val, (G_val@mean_beta)+(C_val@self.gamma)+self.intercept)
                    self.val_r2_list.append(val_r2)
                    if trial is not None:
                        trial.report(val_r2, self.total_epochs_trained)
                    if val_r2 > self.best_r2:
                        self.best_r2 = val_r2
                        self.best_r2_mean_beta = self.mean_beta
                        self.best_r2_intercept = self.intercept
                        self.best_r2_gamma = self.gamma
                        self.best_r2_var_beta = self.var_beta
                        self.best_r2_base_var = self.base_var
                        self.best_r2_last_layer_bias = self.var_pred_model.layers[-2].constant.item() if hasattr(self.var_pred_model.layers[-2], "constant") else self.var_pred_model.layers[-2].bias.item()
                        self.best_r2_epoch = self.total_epochs_trained
                        self.best_r2_fE = self.prior_var
                        
                        # Log many things
                        if logging:
                            wandb.run.summary["best_r2"] = self.best_r2
                print(f"After epoch {self.total_epochs_trained}: Train loss: {round(self.epoch_loss_list[-1], 5)}, Validation loss: {round(self.val_list[-1], 5) if n_val_samples else ''}, Train r2: {round(train_r2, 5)}, Val r2: {round(val_r2, 5) if n_val_samples else ''}, Epoch time: {round(time.time()-epoch_start_time,2)}s")
                if logging:
                    wandb.log({"epoch_train_loss": self.epoch_loss_list[-1], "epoch_val_loss": self.val_list[-1] if n_val_samples else None, "epoch_mean_beta": mean_beta, "epoch_var_prior": self.gE_by_epoch[-1], "epoch_intercept": self.intercept_by_epoch[-1], "epoch_train_r2": train_r2, "epoch_val_r2": val_r2 if n_val_samples else None, "epoch_val_r2_delta": (val_r2-common_r2)/common_r2 if n_val_samples else None, "epoch_var": np.exp(self.var), "epoch_last_layer_bias": self.var_pred_model.layers[-2].constant.item() if hasattr(self.var_pred_model.layers[-2], "constant") else self.var_pred_model.layers[-2].bias.item(), "epoch_base_var":self.base_var})
            
            # Skip computing posterior and R^2 on validation
            else:
                print(f"After epoch {self.total_epochs_trained}: Train loss: {round(self.epoch_loss_list[-1], 5)}, Validation loss: {round(self.val_list[-1], 5) if n_val_samples else ''}, Epoch time: {round(time.time()-epoch_start_time,2)}s")
            self.total_epochs_trained += 1
            
            
            if n_val_samples and args["early_stopping"]:
                if early_stopper.early_stop(val_r2):
                    print("Early stop")
                    break
            
            #if (trial is not None) and (trial.should_prune()):
            #    if logging:
            #        wandb.finish()
            #    raise optuna.TrialPruned()
        
        # Compute posterior using train and validation data to evaluate model on test set.
        if n_val_samples:
            mean_beta, var_beta = self.recompute_posterior(torch.cat([G_torch, G_val_torch]), torch.from_numpy(self.gE_by_epoch[self.best_r2_epoch].astype("float32")).to(device), torch.from_numpy(self.var_by_epoch[self.best_r2_epoch].astype("float32")).to(device), torch.FloatTensor([self.intercept_by_epoch[self.best_r2_epoch]]).to(device), torch.from_numpy(self.best_r2_gamma.astype("float32")).to(device), torch.cat([C_torch, C_val_torch]), torch.cat([y_torch, y_val_torch]), device=None)
            self.best_r2_trainval_mean_beta = mean_beta
            self.best_r2_trainval_var_beta = var_beta

            mean_beta, var_beta = self.recompute_posterior(torch.cat([G_torch, G_val_torch]), torch.from_numpy(self.gE_by_epoch[self.best_loss_epoch].astype("float32")).to(device), torch.from_numpy(self.var_by_epoch[self.best_loss_epoch].astype("float32")).to(device), torch.FloatTensor([self.intercept_by_epoch[self.best_loss_epoch]]).to(device), torch.from_numpy(self.best_loss_gamma.astype("float32")).to(device), torch.cat([C_torch, C_val_torch]), torch.cat([y_torch, y_val_torch]), device=None)
            self.best_loss_trainval_mean_beta = mean_beta
            self.best_loss_trainval_var_beta = var_beta

if __name__ == "__main__":
    pass
