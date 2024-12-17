# --- Configuration ---
configfile: "config.yaml"

# --- Rules ---

# --- One Rule to Rule Them All ---
rule all:
    input:
        f"{config['OUTPUT_DIR']}/{config['STUDY_NAME']}_genes_extended.pq",
        f"{config['OUTPUT_DIR']}/{config['STUDY_NAME']}_predictions_extended.pq"
        # expand(f"{config['OUTPUT_DIR']}/{{trait}}_mean_betas.pq", trait=config.get("TRAITS", [])),
        # expand(f"{config['OUTPUT_DIR']}/{{trait}}_bayes_pred.pq", trait=config.get("TRAITS", [])),

# --- Train Model ---
rule train_model_trait:
    output:
        result_betas=f"{config['OUTPUT_DIR']}/{{trait}}_mean_betas.pq",
        result_preds=f"{config['OUTPUT_DIR']}/{{trait}}_bayes_pred.pq"
    params:
        epochs=25 # Default parameters
    threads: 16
    resources:
        mem_mb=128000,
        gpu="a40:1",
        exclude="ouga12"
    log:
        f"{config['WANDB_LOGS_DIR']}/{{trait}}_train_model.log"
    shell:
        """
        scripts/train_model.py --trait {wildcards.trait} \
                               --epochs {params.epochs} \
                               --study_name {config['STUDY_NAME']} \
                               --output_dir {config['OUTPUT_DIR']} \
                               --hyperparam_version {config['HO_VERSION']} \
                               --optuna_dir {config['OPTUNA_DIR']} \
                               --logs_dir {config['WANDB_LOGS_DIR']} \
                               > {log} 2>&1
        """

# --- Compile Results ---
rule compile_results:
    input:
        all_beta_results=expand(f"{config['OUTPUT_DIR']}/{{trait}}_mean_betas.pq", trait=config.get("TRAITS", [])),
        all_pred_results=expand(f"{config['OUTPUT_DIR']}/{{trait}}_bayes_pred.pq", trait=config.get("TRAITS", [])),
    output:
        compiled_genes_results=f"{config['OUTPUT_DIR']}/{config['STUDY_NAME']}_genes_extended.pq",
        compiled_pheno_results=f"{config['OUTPUT_DIR']}/{config['STUDY_NAME']}_predictions_extended.pq"
    log:
        f"{config['LOGS_DIR']}/compile_results.log"
    shell:
        """
        scripts/consolidate_results.py --output_genes {output.compiled_genes_results} \
                                       --output_pheno {output.compiled_pheno_results} \
                                       --study_name {config['STUDY_NAME']} \
                                       --embedding {config['EMBEDDING']} \
                                       --genotype {config['GENOTYPE']} \
                                       --test_split_size {config['TEST_SPLIT']} \
                                       --logs_dir {config['LOGS_DIR']} \
                                       --all_beta_results "{input.all_beta_results}" \
                                       --all_pred_results "{input.all_pred_results}" \
                                       --traits '{config["TRAITS"]}' \
                                       > {log} 2>&1
        """


# --- Create Plots ---
# --- Figure 3 ---
# rule create_plot_1:
#     input:
#         model = f"{config['TRAINING_DIR']}/{{trait}}/model.pkl" #Update if hyperopt was used. Select the best model.
#     output:
#         plot1=f"{config['PLOTS_DIR']}/plot1.pdf"
#     script:
#         "scripts/create_plot_1.py"
