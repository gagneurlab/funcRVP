library(data.table)
library(ggplot2)
library(ggpubr)

# # Compute R2 -------------------------------------------------
r2_general <-function(preds,actual){
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

# dt_lm <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_deepRVAT_0.25_predictions_NEWsplit.pq') %>% as.data.table()
dt_lm <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/ukbb_wes_500k_DeepRVAT_final_090924_medianshifted/rvat/all_traits_phenopred_0.05.pq') %>% as.data.table()
dt_lm <- dt_lm[model=='lm_sign_genes', .(r2=r2_general(best_prediction, trait_measurement)), by=trait]
dt_lm <- dt_lm[, model:='lm_sign_genes_deepRVAT']

dt_cov <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_traits_covariates_only_phenopred_filteredv3.pq') %>% as.data.table()
dt_cov <- dt_cov[model == 'lm_cov', .(cov_r2=r2_general(best_prediction, trait_measurement)), by=trait]
# dt_cov[, model:='lm_cov']

# dt_bayes <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_bayes <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_bayes <- dt_bayes[, .(r2=r2_general(best_prediction, trait_measurement)), by=trait]
dt_bayes[, model:='FuncRVP']

## Phenotype correlations
# traits list to remove high correlation traits:
# drop_list <- c('Apolipoprotein_B', 'Creatinine', 'Apolipoprotein_A', 'Cholesterol', 'Reticulocyte_count', 'Total_bilirubin', 'Platelet_count', 'Mean_corpuscular_volume')
drop_list <- c('Calcium', 'Cholesterol', 'Creatinine', 'Cystatin_C', 'HDL_cholesterol', 'High_light_scatter_reticulocyte_count', 'High_light_scatter_reticulocyte_percentage', 'LDL_direct', 'Lymphocyte_percentage', 'Mean_corpuscular_haemoglobin', 'Mean_corpuscular_volume', 'Mean_sphered_cell_volume', 'Platelet_crit','Reticulocyte_count','Total_bilirubin')

# Generated random Embeddings
# dt_rand <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_randEmb_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_rand <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_randEmb_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_rand <- dt_rand[, .(r2=r2_general(best_prediction, trait_measurement)), by=trait]
dt_rand[, model:='rand_emb']

rand_plot <- rbindlist(list(dt_bayes, dt_lm, dt_rand), , use.names=TRUE, fill = TRUE)
rand_plot <- merge(rand_plot, dt_cov, by = 'trait')
rand_plot[, delta_r2 := r2 - cov_r2]
rand_plot[, rel_delta_r2 := delta_r2/cov_r2]

pval_rand_a <- wilcox_test(rand_plot[(trait %in% drop_list) & (model %in% c('lm_sign_genes_deepRVAT', 'rand_emb'))], rel_delta_r2 ~ model, paired = T)$p
pval_rand_b <- wilcox_test(rand_plot[(trait %in% drop_list) & (model %in% c('FuncRVP', 'rand_emb'))], rel_delta_r2 ~ model, paired = T)$p
plot_rand <- dcast(rand_plot[, .(trait, model, rel_delta_r2)], ...~model, value.var = 'rel_delta_r2')

rand_a <- ggplot(plot_rand) +
  geom_point(aes(x=lm_sign_genes_deepRVAT, y=rand_emb)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.02, y=0.04), label=paste0('Wilcoxon p-value\n',pval_rand_a), check_overlap = T, color='black') +
  xlim(rand_plot[, min(rel_delta_r2)], rand_plot[, max(rel_delta_r2)]) +
  ylim(rand_plot[, min(rel_delta_r2)], rand_plot[, max(rel_delta_r2)]) +
  xlab(TeX("Linear model on significant genes (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, Random Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()

rand_b <- ggplot(plot_rand) +
  geom_point(aes(x=FuncRVP, y=rand_emb)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.02, y=0.04), label=paste0('Wilcoxon p-value\n',pval_rand_b), check_overlap = T, color='black') +
  xlim(rand_plot[, min(rel_delta_r2)], rand_plot[, max(rel_delta_r2)]) +
  ylim(rand_plot[, min(rel_delta_r2)], rand_plot[, max(rel_delta_r2)]) +
  xlab(TeX("FuncRVP, Omics+PoPS (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, Random Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()

ggarrange(rand_a, rand_b, nrow = 1, labels = c("A", "B"), align = 'h')

# No Embedding
# dt_ridge <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_noemb_deepRVAT_testsplit0.25_noemb_predictions_extended.pq') %>% as.data.table()
dt_ridge <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_noemb_deepRVAT_testsplit0.25_noemb_predictions_extended.pq') %>% as.data.table()
dt_ridge <- dt_ridge[, .(r2=r2_general(best_prediction, trait_measurement)), by=trait]
dt_ridge[, model:='no_emb']

ridge_plot <- rbindlist(list(dt_bayes, dt_lm, dt_ridge), , use.names=TRUE, fill = TRUE)
ridge_plot <- merge(ridge_plot, dt_cov, by = 'trait')
ridge_plot[, delta_r2 := r2 - cov_r2]
ridge_plot[, rel_delta_r2 := delta_r2/cov_r2]
ridge_plot

pval_ridge_a <- wilcox_test(ridge_plot[model %in% c('lm_sign_genes_deepRVAT', 'no_emb')], rel_delta_r2 ~ model, paired = T)$p
pval_ridge_b <- wilcox_test(ridge_plot[model %in% c('FuncRVP', 'no_emb')], rel_delta_r2 ~ model, paired = T)$p
plot_ridge <- dcast(ridge_plot[, .(trait, model, rel_delta_r2)], ...~model, value.var = 'rel_delta_r2')

ridge_a <- ggplot(plot_ridge) +
  geom_point(aes(x=lm_sign_genes_deepRVAT, y=no_emb)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.01, y=0.04), label=paste0('Wilcoxon p-value\n',pval_ridge_a), check_overlap = T, color='black') +
  xlim(ridge_plot[, min(rel_delta_r2)], ridge_plot[, max(rel_delta_r2)]) +
  ylim(ridge_plot[, min(rel_delta_r2)], ridge_plot[, max(rel_delta_r2)]) +
  xlab(TeX("Linear model on significant genes (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, No Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()

ridge_b <- ggplot(plot_ridge) +
  geom_point(aes(x=FuncRVP, y=no_emb)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.01, y=0.04), label=paste0('Wilcoxon p-value\n',pval_ridge_b), check_overlap = T, color='black') +
  xlim(ridge_plot[, min(rel_delta_r2)], ridge_plot[, max(rel_delta_r2)]) +
  ylim(ridge_plot[, min(rel_delta_r2)], ridge_plot[, max(rel_delta_r2)]) +
  xlab(TeX("FuncRVP, Omics+PoPS (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, No Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()

ggarrange(ridge_a, ridge_b, nrow = 1, labels = c("A", "B"), align = 'h')


# Permuted Embeddings
# path_prefix <- '/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_shuffEmb_shuffledemb' 
path_prefix <- '/s/project/geno2pheno/funcrvp/paper_results/paper_funcrvp/v1NEWsplit_shuffEmb_shuffledemb' 
path_sufffix <- '_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq'
perm_ls <- lapply(0:9, FUN = function (i) {
  dt_perm <- read_parquet(paste0(path_prefix, i, path_sufffix)) %>% as.data.table()
  a <- dt_perm[, .(r2=r2_general(best_r2_pred, trait_measurement)), by=trait]
  a[, model:=paste0('perm_emb_',i)]
  })
dt_perm_all <- rbindlist(perm_ls)

dt_perm <- dt_perm_all[, .(r2 = mean(r2), sd_r2 = sd(r2)), by=trait]
dt_perm[, model:='perm_emb']

perm_plot <- rbindlist(list(dt_bayes, dt_lm, dt_perm), , use.names=TRUE, fill = TRUE)
perm_plot <- merge(perm_plot, dt_cov, by = 'trait')
perm_plot[, delta_r2 := r2 - cov_r2]
perm_plot[, rel_delta_r2 := delta_r2/cov_r2]
perm_plot
pval_perm_a <- wilcox_test(perm_plot[model %in% c('lm_sign_genes_deepRVAT', 'perm_emb')], rel_delta_r2 ~ model, paired = T)$p
pval_perm_b <- wilcox_test(perm_plot[model %in% c('FuncRVP', 'perm_emb')], rel_delta_r2 ~ model, paired = T)$p
plot_perm <- dcast(perm_plot[, .(trait, model, rel_delta_r2)], ...~model, value.var = 'rel_delta_r2')
plot_perm <- merge(plot_perm, dt_perm[, .(trait, sd_r2)], by='trait')

perm_a <- ggplot(plot_perm, aes(x=lm_sign_genes_deepRVAT, y=perm_emb)) +
  geom_point(size=0.75) +
  geom_linerange(aes(ymin=perm_emb-1.96*sd_r2, ymax=perm_emb+1.96*sd_r2)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.02, y=0.04), label=paste0('Wilcoxon p-value\n',pval_perm_a), check_overlap = T, color='black') +
  xlim(perm_plot[, min(rel_delta_r2)], perm_plot[, max(rel_delta_r2)]) +
  ylim(perm_plot[, min(rel_delta_r2)], perm_plot[, max(rel_delta_r2)]) +
  xlab(TeX("Linear model on significant genes (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, Permuted Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()
perm_a

perm_b <- ggplot(plot_perm, aes(x=FuncRVP, y=perm_emb)) +
  geom_point(size=0.75) +
  geom_linerange(aes(ymin=perm_emb-1.96*sd_r2, ymax=perm_emb+1.96*sd_r2)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_text(aes(x=0.02, y=0.04), label=paste0('Wilcoxon p-value\n',pval_perm_b), check_overlap = T, color='black') +
  xlim(perm_plot[, min(rel_delta_r2)], perm_plot[, max(rel_delta_r2)]) +
  ylim(perm_plot[, min(rel_delta_r2)], perm_plot[, max(rel_delta_r2)]) +
  xlab(TeX("FuncRVP, Omics+PoPS (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, Permuted Embedding (Relative $\\Delta R^2$)")) +
  theme_cowplot()

ggarrange(perm_a, perm_b, nrow = 1, labels = c("A", "B"), align = 'h')


# All box plots
dt <- rbindlist(list(dt_bayes, dt_lm, dt_ridge, dt_rand, dt_perm), , use.names=TRUE, fill = TRUE)
dt <- merge(dt, dt_cov, by = 'trait')
dt[, delta_r2 := r2 - cov_r2]
dt[, rel_delta_r2 := delta_r2/cov_r2]
dt

dt %>%
  wilcox_test(rel_delta_r2 ~ model, paired = T) %>%
  add_significance("p")

ggplot(dt, aes(x=reorder(model, rel_delta_r2, func_c='median'), y=rel_delta_r2)) +
  scale_y_sqrt(na.value=0) +
  geom_boxplot() +
  guides(x = guide_axis(angle=45)) +
  xlab('Models') +
  ylab(TeX("Relative $\\Delta R^2$")) +
  theme_cowplot()

