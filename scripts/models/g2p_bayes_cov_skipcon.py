import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lr_scheduler
from torch.distributions.multivariate_normal import MultivariateNormal
from sklearn.metrics import r2_score
import wandb
import time
import warnings
import logging
import sys
from typing import Optional
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression



# Ignoring sparse warning
warnings.filterwarnings("ignore")
torch.set_float32_matmul_precision("high")
torch.backends.cuda.matmul.allow_tf32 = True

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

def tensor_to_numpy(tensor):
    """Safely convert any tensor to numpy array"""
    return tensor.detach().cpu().to(torch.float32).numpy()

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
        self.constant = nn.Parameter(
            torch.tensor(float(constant_value)), requires_grad=True
        )

    def forward(self, x):
        return self.constant.expand(x.shape[0]).unsqueeze(1)


# Class that models variance of beta as a function of the embedding
class VarPredModel(nn.Module):
    def __init__(
        self,
        emb_dim: int,
        n_cov: int,
        n_hidden: Optional[int]=0,
        hiddem_dim: Optional[int]=None,
        last_layer_bias: Optional[float]=None,
        nonlinearity: Optional[str]="relu",
        final_nonlinearity: Optional[str]="softplus",
        y_var_init: Optional[float]=1e-3,
        gene_var_init: Optional[float]=None,
        base_var_init: Optional[float]=5e-5,
        base_var_const: Optional[float]=0,
        batch_norm: Optional[str]="linear",
        n_genes: Optional[int]=None,
        skip_con: Optional[bool]=False,
        device: Optional[str]=None,
    ):
        super(VarPredModel, self).__init__()

        # Define activation functions using a dictionary for conciseness
        activation_functions = {
            "relu": nn.ReLU,
            "silu": nn.SiLU,
            "gelu": nn.GELU,
            "softplus": nn.Softplus,
            "elu": nn.ELU,
            "celu": nn.CELU,
            "identity": nn.Identity,
        }
        self.activation = activation_functions.get(nonlinearity, nn.Identity)()

        final_activation_functions = {
            "softplus": nn.Softplus,
            "sigmoid": nn.Sigmoid,
        }
        FinalActivation = final_activation_functions.get(final_nonlinearity, nn.Softplus)

        # Determine hidden dimension, default to embedding dimension if not provided
        hidden_dim = hiddem_dim or emb_dim

        # Define layers using a more concise structure
        self.layers = nn.ModuleList()
        if n_hidden >= 0:
            if n_hidden > 0:
                # Input layer
                self.layers.append(nn.Linear(emb_dim, hidden_dim))
                if batch_norm == "linear":
                    self.layers.append(nn.BatchNorm1d(hidden_dim))

                # Hidden layers
                for _ in range(n_hidden - 1):
                    self.layers.append(self.activation)
                    if batch_norm == "activation":
                        self.layers.append(nn.BatchNorm1d(hidden_dim))
                    self.layers.append(nn.Linear(hidden_dim, hidden_dim))
                    if batch_norm == "linear":
                        self.layers.append(nn.BatchNorm1d(hidden_dim))

                # Output layer
                self.layers.append(self.activation)
                if batch_norm == "activation":
                    self.layers.append(nn.BatchNorm1d(hidden_dim))
                self.layers.append(nn.Linear(hidden_dim, 1))

            else: # n_hidden == 0, single linear layer
                self.layers.append(nn.Linear(emb_dim, 1))

            if last_layer_bias is not None: # Use 'is not None' for boolean/float optionals
                self.layers[-1].bias.data.fill_(last_layer_bias)
        
        # n_hidden < 0, Constant Variance
        else:
            class ConstantModule(nn.Module): # Define ConstantModule locally if it's simple
                def __init__(self, value):
                    super().__init__()
                    self.value = nn.Parameter(torch.FloatTensor([value]))
                def forward(self, input):
                    return self.value.expand(input.size(0), 1) # expand to batch size

            self.layers.append(
                ConstantModule(last_layer_bias if last_layer_bias is not None else 0)
            )

        self.layers.append(FinalActivation()) # Final non-linearity layer

        # Linear Residual connection
        self.skip_con = skip_con
        if skip_con:
            self.shortcut = nn.Sequential(nn.Linear(emb_dim, 1), FinalActivation())
            if last_layer_bias is not None:
                self.shortcut[0].bias.data.fill_(last_layer_bias)

        # Initialize parameters
        self.log_var = nn.Parameter(torch.FloatTensor([np.log(y_var_init)]))
        self.intercept = nn.Parameter(torch.FloatTensor([0]))
        self.gamma = nn.Parameter(torch.zeros(n_cov, device=device)) # Directly set device here

        self.gene_var = nn.Parameter(gene_var_init * torch.ones(n_genes, 1, device=device)) if n_genes and gene_var_init else None

        self.log_base_var = nn.Parameter(torch.FloatTensor([np.log(base_var_init)], device=device))
        self.base_var_const = base_var_const

    def forward(self, emb):
        if self.skip_con:
            emb_copy = emb.clone().detach() # More explicit detach

        for layer in self.layers:
            emb = layer(emb)

        if self.skip_con:
            emb = emb + self.shortcut(emb_copy)

        if self.gene_var is not None:
            emb = emb + torch.exp(self.gene_var)

        emb = emb + torch.exp(self.log_base_var) + self.base_var_const

        return (
            emb,
            self.log_var,
            self.intercept,
            self.gamma,
            torch.exp(self.log_base_var),
        )


# Class that fits the genotype to phenotype model
class G2P_Model(nn.Module):

    def __init__(
        self,
        emb_dim: int,
        n_cov: int,
        n_hidden: Optional[int]=0,
        hiddem_dim: Optional[int]=None,
        last_layer_bias: Optional[float]=None,
        nonlinearity: Optional[str]="relu",
        final_nonlinearity: Optional[str]="softplus",
        y_var_init: Optional[float]=1e-3,
        gene_var_init: Optional[float]=None,
        base_var_init: Optional[float]=1e-7,
        base_var_const: Optional[float]=0,
        batch_norm: Optional[str]="linear",
        n_genes: Optional[int]=None,
        skip_con: Optional[bool]=False,
        device: Optional[str]=None,
    ):
        super(G2P_Model, self).__init__()

        # Initialize the f(E) model
        self.var_pred_model = VarPredModel(
            emb_dim,
            n_cov,
            n_hidden=n_hidden,
            hiddem_dim=hiddem_dim,
            last_layer_bias=last_layer_bias,
            nonlinearity=nonlinearity,
            final_nonlinearity=final_nonlinearity,
            y_var_init=y_var_init,
            gene_var_init=gene_var_init,
            base_var_init=base_var_init,
            base_var_const=base_var_const,
            batch_norm=batch_norm,
            n_genes=n_genes,
            skip_con=skip_con,
            device=None,
        )

        self.total_epochs_trained = 0
        self.best_r2 = -np.inf
        self.best_loss = np.inf

        self.train_loss_list = []
        self.val_loss_list = []
        self.train_r2_list = []
        self.val_r2_list = []

    # Compute likelihood of y, and get outputs of the variance prediction model as self attributes
    def forward(self, G, emb, C):
        # import ipdb; ipdb.set_trace()
        self.prior_var, self.var, self.intercept, self.gamma, self.base_var = (
            self.var_pred_model(emb)
        )
        cov = (G * self.prior_var.squeeze(1)) @ G.transpose(1, 0) + torch.diag(
            torch.exp(self.var).expand(G.shape[0])
        )
        pred = MultivariateNormal(
            loc=(C @ self.gamma) + self.intercept.expand(G.shape[0]),
            scale_tril=torch.linalg.cholesky(cov),
        )
        return pred

    # Phenotype prediction using posterior beta and covariates
    def predict(self, G, C, beta=None, gamma=None, intercept=None):
        if beta and gamma and intercept:
            return (G @ beta) + (C @ gamma) + intercept
        else:
            return (G @ self.posterior_beta) + (C @ self.gamma) + self.intercept


    # Faster computation of posterior mean
    def _get_posterior_mean_faster(self, G, C, y, GT_G, gE, gamma, intercept, var):
        with torch.no_grad():
            sigma_inv = torch.diag((1/gE).squeeze()) + (1/torch.exp(var))[0] * GT_G
            x = (1/torch.exp(var))[0] * torch.matmul(torch.transpose(G, 0, 1), (y - ((C @ gamma) + intercept)))
            mu_hat_unscaled = torch.linalg.solve(sigma_inv, x)
            mean_beta = mu_hat_unscaled  
        return mean_beta
    
    def _get_posterior_params(self, G, C, y, GT_G, gE, gamma, intercept, var, faster=True):
        if faster:
            return self._get_posterior_mean_faster(G, C, y, GT_G, gE, gamma, intercept, var), None
        else:
            with torch.no_grad():
                sigma_inv = torch.diag((1 / gE).squeeze()) + (1 / torch.exp(var))[0] * GT_G
                sigma = torch.linalg.inv(sigma_inv)

                mean_beta = sigma @ (
                    (1 / torch.exp(var))[0] * torch.transpose(G, 0, 1)
                    @ (y - ((C @ gamma) + intercept))
                )
                var_beta = torch.diag(sigma)
            return tensor_to_numpy(mean_beta), tensor_to_numpy(var_beta)

    def _update_best_params(self):
        """Updates best params if current validation metric is better."""
        # TODO, call recompute posterior to get mean_beta
        self.best_mean_beta = tensor_to_numpy(self.posterior_beta)
        #self.best_var_beta = self.posterior_beta
        self.best_var = tensor_to_numpy(self.var)
        self.best_gamma = tensor_to_numpy(self.gamma)
        self.best_intercept = tensor_to_numpy(self.intercept)
        self.best_base_var = tensor_to_numpy(self.base_var)
        self.best_last_layer_bias = tensor_to_numpy(
            self.var_pred_model.layers[-2].constant
            if hasattr(self.var_pred_model.layers[-2], "constant")
            else self.var_pred_model.layers[-2].bias
        )
        self.best_epoch = self.total_epochs_trained
        self.best_prior_var = tensor_to_numpy(self.prior_var)

    def _log_epoch_wandb(self, epoch, train_r2, val_r2, common_r2, n_val_samples):
        """Logs epoch metrics to WandB."""
        wandb.log(
            {   
                "epoch": epoch,
                "epoch_train_loss": self.train_loss_list[-1],
                "epoch_val_loss": (self.val_loss_list[-1] if n_val_samples else None),
                "epoch_posterior_beta": tensor_to_numpy(self.posterior_beta),
               # "epoch_posterior_var_beta": np.squeeze(tensor_to_numpy(self.posterior_var_beta)),
                "epoch_prior_var": tensor_to_numpy(self.prior_var),
                "epoch_gamma": tensor_to_numpy(self.gamma),
                "epoch_intercept": tensor_to_numpy(self.intercept),
                "epoch_train_r2": train_r2,
                "epoch_val_r2": val_r2 if n_val_samples else None,
                "epoch_val_r2_delta": (
                    ((val_r2 - common_r2) / common_r2) if n_val_samples else None
                ),
                "epoch_var": np.exp(tensor_to_numpy(self.var)),
                "epoch_last_layer_bias": (
                    self.var_pred_model.layers[-2].constant.item()
                    if hasattr(self.var_pred_model.layers[-2], "constant")
                    else self.var_pred_model.layers[-2].bias.item()
                ),
                "epoch_base_var": tensor_to_numpy(self.base_var),
            }
        )

    # Function that runs training loop and logs parameters and metrics
    def fit_model(
        self,
        args: dict,
        emb: np.ndarray,
        G: np.ndarray,
        C: np.ndarray,
        y: np.ndarray,
        G_val: Optional[np.ndarray]=None,
        C_val: Optional[np.ndarray]=None,
        y_val: Optional[np.ndarray]=None,
        logging: Optional[bool]=False,
        device: Optional[str]=None,
    ):

        # If you want to normalize embeddings along each dimension
        follow_metric = args.get("follow_metric", "loss")

        # TODO Define defaults for params read from the config
        if "batch_size_schedule" not in args:
            args["batch_size_schedule"] = {0: args["batch_size"]}

        if "learning_rate_schedule" not in args:
            args["learning_rate_schedule"] = {0: args["learning_rate"]}
        if "early_stopping" not in args:
            args["early_stopping"] = False

        logger.info("Initializing model params with linear regression on covariates")
        # Fit Linear model on the covariates
        lm = LinearRegression().fit(C, y)

        # Get R^2 of the linear model on the validation set
        if (G_val is not None) and (C_val is not None) and (y_val is not None):
            common_r2 = r2_score(y_val, lm.predict(C_val))
            logger.info(f"Covariates r2 on val: {common_r2}")

        # Initialize covariate effects in FuncRVP = effect of covariates from the LM
        self.var_pred_model.gamma.data = torch.Tensor(lm.coef_).to(device)

        # Initialize intercept in FuncRVP = effect of covariates from the LM
        self.var_pred_model.intercept.data = torch.Tensor([lm.intercept_]).to(device)

        logger.info("Uploading data to GPU")
        emb_torch = torch.Tensor(emb).to(device)

        #TODO seems redundant
        G = torch.Tensor(G)
        C = torch.Tensor(C)
        y = torch.Tensor(y)
        G_val = torch.Tensor(G_val)
        C_val = torch.Tensor(C_val)
        y_val = torch.Tensor(y_val)

        G_torch = G.to(device)
        C_torch = C.to(device)
        y_torch = y.to(device)
        dataset_train = torch.utils.data.TensorDataset(G_torch, C_torch, y_torch) 
        n_train_samples = len(dataset_train)

        if (G_val is not None) and (C_val is not None) and (y_val is not None):
            G_val_torch = G_val.to(device)
            C_val_torch = C_val.to(device)
            y_val_torch = y_val.to(device)
            dataset_val = torch.utils.data.TensorDataset(G_val_torch, C_val_torch, y_val_torch)
            n_val_samples = len(dataset_val)
        else:
            n_val_samples = 0

        # Computer GT^G for posterior computation
        # GT_G = (torch.transpose(G_torch, 0, 1) @ G_torch).cpu()
        GT_G = (torch.transpose(G, 0, 1) @ G)

        # Initialize the Adam optimizer with learning rate and weight decay as input arguments from the user
        optimizer = torch.optim.Adam(
            self.parameters(),
            lr=args["learning_rate_schedule"][0],
            weight_decay=args["weight_decay"],
        )
        batch_size = args["batch_size"]

        if logging:
            wandb.run.summary["num_parameters"] = sum(
                p.numel() for p in self.parameters() if p.requires_grad
            )
            wandb.run.summary["covariate_r2"] = common_r2

        # PyTorch speed-up trick
        scaler = torch.cuda.amp.GradScaler()

        # Instantiate early stopper
        if args["early_stopping"]:
            early_stopper = EarlyStopper(patience=5, baseline=common_r2)

        # START training
        logger.info("---------------- Starting training ----------------")
        for epoch in range(0, args["epochs"]):
            self.train()

            # if epoch in args["batch_size_schedule"]:
            #     optimizer = torch.optim.Adam(
            #         self.parameters(),
            #         lr=args["learning_rate"],
            #         weight_decay=args["weight_decay"],
            #     )
            #     scheduler = lr_scheduler.LinearLR(
            #         optimizer, start_factor=1.0, end_factor=0.1, total_iters=10
            #     )
            #     batch_size = args["batch_size_schedule"][epoch]

            # if epoch in args["learning_rate_schedule"]:
            #     optimizer = torch.optim.Adam(
            #         self.parameters(),
            #         lr=args["learning_rate_schedule"][epoch],
            #         weight_decay=args["weight_decay"],
            #     )

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
                indices = permutation[i : i + batch_size]
                G_batch, C_batch, y_batch = dataset_train[indices]

                # PyTorch speed-up trick
                with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                # Compute -log(likelihood). Here "pred" is an object of the class MultivariateNormal
                    pred = self.forward(G_batch, emb_torch, C_batch)
                    loss = -pred.log_prob(y_batch).mean()

                    # Add extra regularization to the loss
                    if "alpha_L1_fE" in args:
                        loss += args["alpha_L1_fE"] * torch.norm(
                                self.var_pred_model(emb_torch)[0], p=1
                            )
                    if "alpha_L2_fE" in args:
                        loss += args["alpha_L2_fE"] * torch.norm(
                                self.var_pred_model(emb_torch)[0], p=2
                            )
                    if "alpha_L1_gene_var" in args and (not (self.var_pred_model.gene_var is None)):
                        loss += args["alpha_L1_gene_var"] * torch.norm(
                                self.var_pred_model.gene_var, p=1
                            )

                # PyTorch speed-up trick
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()

                # Record epoch loss
                epoch_loss += loss.detach().item() * len(y_batch)
                self.train_loss_list.append(loss.detach().item() / len(indices))

            self.train_loss_list.append(epoch_loss / n_train_samples)

            logger.info(f"Epoch {epoch}: Train loss: {self.train_loss_list[-1]}")

            # Initiate evaluation mode (no gradients)
            self.eval()
            val_loss = 0
            # Loop over all validation set batches
            for i in range(0, n_val_samples, batch_size):

                G_batch, C_batch, y_batch = dataset_val[i : i + batch_size]

                # Compute prediction and -log(likelihood) on validation data
                with torch.no_grad(): 
                    pred = self.forward(G_batch, emb_torch, C_batch)
                    val_loss += -pred.log_prob(
                        y_batch
                    ).mean().detach().item() * len(y_batch)

            # Record loss on validation set
            if n_val_samples:
                self.val_loss_list.append(val_loss / n_val_samples)
                if (
                    self.val_loss_list[-1] < self.best_loss
                ) and follow_metric == "loss":
                    #self._update_best_params()
                    self.best_loss = self.val_loss_list[-1]

            logger.info(f"Epoch {epoch}: Val loss: {self.val_loss_list[-1]}")

            # Compute R^2 on train and validation set
            if (follow_metric == "r2") or (epoch + 1 == args["epochs"]):
                self.posterior_beta, _ = (
                    self._get_posterior_params(
                        G,
                        C,
                        y,
                        GT_G,
                        self.prior_var.detach().cpu(),
                        self.gamma.detach().cpu(),
                        self.intercept.detach().cpu(),
                        self.var.detach().cpu(),
                        faster=True,
                    )
                )
                
                self._update_best_params()
                train_r2 = r2_score(
                    y.to(torch.float32).numpy(), ((G @ self.posterior_beta) + (C @ self.gamma.detach().cpu()) + self.intercept.detach().cpu()).to(torch.float32).numpy()
                )
                self.train_r2_list.append(train_r2)

                if n_val_samples:
                    val_r2 = r2_score(
                        y_val.to(torch.float32).numpy(),
                        ((G_val @ self.posterior_beta) + (C_val @ self.gamma.detach().cpu()) + self.intercept.detach().cpu()).to(torch.float32).numpy()
                    )
                    self.val_r2_list.append(val_r2)
                    # if trial is not None:
                    #     trial.report(val_r2, self.total_epochs_trained)

                    if (val_r2 > self.best_r2) and follow_metric == "r2":
                        self.best_r2 = val_r2

                        # wandb logging
                        if logging:
                            wandb.run.summary["best_r2"] = self.best_r2

                logger.info(f"Epoch {epoch}: Train r2: {round(train_r2, 5)}, Val r2: {round(val_r2, 5)}")

                if logging:
                    self._log_epoch_wandb(
                        epoch, train_r2, val_r2, common_r2, n_val_samples
                    )  # Call wandb logging function

            # Skip computing posterior and R^2 on validation
            logger.info(f"Epoch time: {round(time.time()-epoch_start_time,2)}s")

            self.total_epochs_trained += 1

            if n_val_samples and args["early_stopping"]:
                if early_stopper.early_stop(val_r2):
                    logger.info("Early stop")
                    break

        del GT_G # Free up memory

        # Compute posterior using train and validation data to evaluate model on test set.
        if n_val_samples:
            # compute GT_G for the whole matrix
            GT_G_trainval = (torch.transpose(torch.cat([G, G_val]), 0, 1) @ torch.cat([G, G_val]))
            self.best_posterior_mean_beta, self.best_posterior_var_beta = (
                self._get_posterior_params(
                    torch.cat([G, G_val]),
                    torch.cat([C, C_val]),
                    torch.cat([y, y_val]),
                    GT_G_trainval,
                    torch.Tensor(self.best_prior_var),
                    torch.Tensor(self.best_gamma),
                    torch.Tensor(self.best_intercept),
                    torch.Tensor(self.best_var),
                    faster=False,
                )
            )


if __name__ == "__main__":
    pass