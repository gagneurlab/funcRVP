import os
import sys
import pandas as pd
import numpy as np
import yaml # Import yaml
import logging
from typing import Optional
from .utils import trait_INT_zscore

# --- Logging Setup ---
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s",
    level=logging.INFO,
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# import ipdb; ipdb.set_trace()
# --- Configuration Loading ---
data_paths_file = './scripts/utils/data_paths.yaml' # Changed to .yaml
if not os.path.exists(data_paths_file):
    logger.error(f"Configuration file not found: {data_paths_file}\nExiting.")
    sys.exit(1)
with open(data_paths_file, 'r') as f: # Open and load yaml
    data_paths = yaml.safe_load(f)


def load_data(
    trait: str,
    embedding_path: str,
    genotype_path: str,
    phenotype_dir: str,
    test_split_size: Optional[float]=0.25,
    val_split_size: Optional[float]=0.1,
    split_seed: Optional[int]=0,
    train_individuals_path: Optional[str]=None,
    covariates_path: Optional[str]=None,
    prs_path: Optional[str]=None,
    normalize_embedding: Optional[bool]=False,
    shuffled_phenotype: Optional[bool]=False,
    shuffled_embedding: Optional[bool]=False,
    random_embedding: Optional[bool]=False,
    use_prs: Optional[bool]=True,
    normalize_covariates: Optional[bool]=True,
    gene_subset: Optional[list]=None,
    dataset_version: Optional[str]="filteredv3", #TODO: add this to file paths
):

    # Read GIS matrix
    if genotype_path:
        if os.path.exists(genotype_path):
            logger.info(f"Reading genotype data from: {genotype_path}")
            if gene_subset:  # Want to run it on a subset of genes
                logger.info(f"Using a subset of genes: {gene_subset}")
                GT_df = pd.read_parquet(genotype_path, columns=gene_subset).reset_index()
            else:
                GT_df = pd.read_parquet(genotype_path).reset_index()
        else:
            logger.error(f"Genotype path not found: {genotype_path}\nExiting.")
            sys.exit(1)
    else:
        logger.error(f"Genotype path not defined in the config file.\nExiting.")
        sys.exit(1)

    gene_list = sorted(list(GT_df.columns))

    # Read gene embeddings    
    if embedding_path:
        if os.path.exists(embedding_path):
            logger.info(f"Reading gene embeddings from: {embedding_path}")
            gene_embeddings_df = (
                pd.read_csv(embedding_path, sep="\t")
                .rename(columns={"gene_id": "gene"})
                .sort_values("gene")
            )
            gene_list = sorted(
                set(gene_embeddings_df.sort_values("gene")["gene"]).intersection(
                    set(gene_list)
                )
            )
            gene_embeddings_df = (
                gene_embeddings_df.set_index("gene").loc[gene_list].reset_index()
            )
            gene_embeddings_df = gene_embeddings_df.set_index("gene")
            gene_embeddings_df = gene_embeddings_df[
                gene_embeddings_df.index.isin(gene_list)
            ]
            gene_embeddings_df = gene_embeddings_df[
                ~gene_embeddings_df.index.duplicated(keep="first")
            ]
            emb = gene_embeddings_df.values
        else:
            logger.info(f"Embedding path does not exist {embedding_path}\nProceeding without embeddings.")
            emb = None
    else:
        logger.info("No embedding type specified, proceeding without embeddings.")
        emb = None

    # Filter GIS to include only those genes whose embedding is available
    GT_df = GT_df[["individual"] + gene_list]

    # Read phenotype
    if phenotype_dir:
        if os.path.exists(phenotype_dir):
            raw_pheno = trait_INT_zscore(trait, phenotype_dir) # Pass pheno_data_dir
        else:
            logger.error(f"Phenotype data path not found: {phenotype_dir}\nExiting.")
            sys.exit(1)
    else:
        logger.error(f"Phenotype data path not defined in the config file.\nExiting.")
        sys.exit(1)

    # Read covariates and PRS
    covariates_path = data_paths['covariate_paths'].get('covariates_path')
    prs_path = data_paths['covariate_paths'].get('prs_path')
    if covariates_path:
        if os.path.exists(covariates_path):
            logger.info(f"Reading covariates from: {covariates_path}")
            # TODO: Add check for what columns are being used form the covariates file.
            covs_dt = pd.read_parquet(covariates_path).dropna()
            logger.info(f"Using covariates: {covs_dt.columns}")
        else:
            logger.error(f"Covariates path not found: {covariates_path}\nProceeding without covariates.")
            covs_dt = pd.DataFrame(columns=["individual"])
    else:
        logger.info("No covariates path specified, proceeding without covariates.")
        covs_dt = pd.DataFrame(columns=["individual"])

    if use_prs:
        if prs_path:
            if os.path.exists(prs_path):
                logger.info(f"Reading PRS from: {prs_path}")
                prs_df = pd.read_parquet(prs_path,
                    columns=["individual", f"{trait}_PRS", f"{trait}_common_resid"],
                ).dropna()

                raw_pheno = raw_pheno.merge(
                    prs_df[["individual", f"{trait}_common_resid"]], on="individual", how="inner"
                )

                try:
                    covs_dt = covs_dt.merge(
                        prs_df[["individual", f"{trait}_PRS"]], on="individual", how="inner"
                    )
                except:
                    covs_dt = prs_df
            else:
                logger.error(f"PRS path not found: {prs_path}\nProceeding without PRS.")
                use_prs = False
        else:
            logger.error(f"PRS path not defined in the config file.\nProceeding without PRS.")
            use_prs = False

    # Get individuals for who all data is available
    try:
        inds = raw_pheno.merge(covs_dt[["individual"]], how="inner").merge(
            GT_df[["individual"]], how="inner"
        )[["individual"]]
    except:
        inds =  GT_df[["individual"]]

    #TODO: add else condition to handle missing train_individuals_path, just skip this whole block.
    # Individuals to be used only in the train split
    # Set a seed for split reproducibility
    np.random.seed(split_seed)
    if train_individuals_path and os.path.exists(train_individuals_path):
        logger.info(f"Reading train individuals list from: {train_individuals_path}")
        train_inds = pd.read_parquet(train_individuals_path).reset_index()[["individual"]]

        # Assign 'train' to rows where the 'group' column is in train_groups
        inds["split"] = np.nan
        inds.loc[inds["individual"].isin(train_inds.individual), "split"] = "train"

        # Compute train/val/test proportions based on some criteria
        test_proportion = (
            test_split_size
            * inds.shape[0]
            / (inds.shape[0] - inds[inds.split == "train"].shape[0])
        )
        val_proportion = (
            val_split_size
            * inds.shape[0]
            / (inds.shape[0] - inds[inds.split == "train"].shape[0])
        )  # Keep 30k samples in val
        remaining_train_proportion = 1 - (test_proportion + val_proportion)

        # Randomly assign splits to the remaining rows
        remaining_rows = inds["split"].isna()
        inds.loc[remaining_rows, "split"] = np.random.choice(
            ["train", "val", "test"],
            size=remaining_rows.sum(),
            p=[remaining_train_proportion, val_proportion, test_proportion],
        )

    else:
        # This block handles the standard random split when no predefined train list is given
        logger.info(f"Train individuals path not found or not specified: {train_individuals_path}. Performing standard random split.")

        # Define target proportions for train, val, test
        if not (0 < test_split_size < 1 and 0 < val_split_size < 1):
            raise ValueError("test_split_size and val_split_size must be between 0 and 1.")
        if test_split_size + val_split_size >= 1.0:
            raise ValueError("The sum of test_split_size and val_split_size must be less than 1.0")

        train_proportion = 1.0 - test_split_size - val_split_size
        split_proportions = [train_proportion, val_split_size, test_split_size]
        split_categories = ['train', 'val', 'test']

        logger.info(f"Assigning splits randomly to all {inds.shape[0]} individuals with proportions: Train={split_proportions[0]:.4f}, Val={split_proportions[1]:.4f}, Test={split_proportions[2]:.4f}")

        # Assign splits randomly based on the overall proportions
        inds["split"] = np.random.choice(
            split_categories,
            size=inds.shape[0], # Assign to all rows
            p=split_proportions, # Use the overall proportions
            replace=True # Each choice is independent
        )

    # Return train/val/test indices
    id_train = inds[inds.split == "train"]["individual"]
    id_val = inds[inds.split == "val"]["individual"]
    id_test = inds[inds.split == "test"]["individual"]

    # Convert all dataframes to numpy arrays for downstream models
    G_train = (
        inds[inds.split == "train"][["individual"]]
        .merge(GT_df)
        .set_index("individual")
        .values
    )
    G_val = (
        inds[inds.split == "val"][["individual"]]
        .merge(GT_df)
        .set_index("individual")
        .values
    )
    G_test = (
        inds[inds.split == "test"][["individual"]]
        .merge(GT_df)
        .set_index("individual")
        .values
    )

    Cov_train = (
        inds[inds.split == "train"][["individual"]]
        .merge(covs_dt)
        .set_index("individual")
        .values
    )
    Cov_val = (
        inds[inds.split == "val"][["individual"]]
        .merge(covs_dt)
        .set_index("individual")
        .values
    )
    Cov_test = (
        inds[inds.split == "test"][["individual"]]
        .merge(covs_dt)
        .set_index("individual")
        .values
    )

    # Fix normalization
    if normalize_covariates:
        cov_trainval_means = np.concatenate([Cov_train, Cov_val]).mean(axis=0)
        cov_trainval_sd = np.concatenate([Cov_train, Cov_val]).std(axis=0)

        Cov_train = (Cov_train - cov_trainval_means) / cov_trainval_sd
        Cov_val = (Cov_val - cov_trainval_means) / cov_trainval_sd
        Cov_test = (Cov_test - cov_trainval_means) / cov_trainval_sd

    # if use_residuals:
    residuals_train = (
        inds[inds.split == "train"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[[f"{trait}_common_resid"]]
        .values.squeeze()
    )
    residuals_val = (
        inds[inds.split == "val"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[[f"{trait}_common_resid"]]
        .values.squeeze()
    )
    residuals_test = (
        inds[inds.split == "test"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[[f"{trait}_common_resid"]]
        .values.squeeze()
    )
    # else:
    trait_measurement_train = (
        inds[inds.split == "train"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[["new_INT"]]
        .values.squeeze()
    )
    trait_measurement_val = (
        inds[inds.split == "val"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[["new_INT"]]
        .values.squeeze()
    )
    trait_measurement_test = (
        inds[inds.split == "test"][["individual"]]
        .merge(raw_pheno)
        .set_index("individual")[["new_INT"]]
        .values.squeeze()
    )

    if shuffled_phenotype:
        np.random.seed(split_seed)
        trait_measurement_train = np.random.permutation(trait_measurement_train)
        trait_measurement_val = np.random.permutation(trait_measurement_val)

    if normalize_embedding:
        emb = emb / np.linalg.norm(emb, axis=0)
        
    if shuffled_embedding:
        np.random.seed(split_seed)
        emb = np.random.permutation(emb)

    if random_embedding:
        emb = np.random.normal(loc=0, scale=1, size=emb.shape)

    return (
        (G_train, G_val, G_test),
        (residuals_train, residuals_val, residuals_test),
        emb,
        gene_list,
        (id_train, id_val, id_test),
        (trait_measurement_train, trait_measurement_val, trait_measurement_test),
        (Cov_train, Cov_val, Cov_test),
    )

