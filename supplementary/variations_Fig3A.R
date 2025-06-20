library(ggplot2)
library(data.table)
library(cowplot)
library(ggrepel)
library(magrittr)
library(BiocParallel)
library(latex2exp)
library(arrow)
library(rstatix)
library(ggpubr)
library(dplyr)

# # Compute R2 -------------------------------------------------
r2_general <-function(preds,actual){
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}


theme_pub <- theme(axis.text=element_text(size=12),
                   text = element_text(size=12))

cols <- palette.colors(8, palette = "Okabe-Ito")

# Compute R2 for model
dt_bayes <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
# dt_bayes <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/cv/ukbb_wes_500k_DeepRVAT_final_090924_medianshifted/pops_mat_pca256_omics/allcv_traits_phenopred.pq') %>% as.data.table()
dt_bayes[, model:='funcrvp_omics_pops']

op_base_dir <- '/s/project/geno2pheno/funcrvp/paper_revisions/predictions/'
# op_base_dir <- '/s/project/geno2pheno/funcrvp/paper_revisions/cv/'
geno_dir <- 'ukbb_wes_500k_DeepRVAT_final_090924_medianshifted/'
rvat_dir <- 'rvat/'
p_thresh <- 0.5
dt_ols <- read_parquet(paste0(op_base_dir, geno_dir, rvat_dir, 'all_traits_phenopred_', p_thresh, '.pq')) %>% as.data.table()
dt_ols[, embedding := 'None']
if("pred" %in% colnames(dt_ols)){
  dt_ols[, best_prediction := pred]
  dt_ols[, pred := NULL]
}

dt_cov <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_traits_covariates_only_phenopred_filteredv3.pq') %>% as.data.table()
# dt_cov <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/cv/ukbb_wes_500k_DeepRVAT_final_090924_medianshifted/rvat_onlycov/allcv_traits_phenopred.pq') %>% as.data.table()
dt_cov[, model := 'lm_cov']

dt <- rbindlist(list(dt_bayes, dt_ols, dt_cov), use.names=TRUE, fill = TRUE)
dt[, `:=` (version=NULL, dataset_version=NULL, embedding=NULL)]

comp_r2 <- dcast(dt[, .(r2 = r2_general(best_prediction, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
comp_r2 <- melt(comp_r2, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
comp_r2[, delta_r2 := r2 - lm_cov]
comp_r2[, rel_delta_r2 := delta_r2/lm_cov]

# write_parquet(comp_r2, '/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_deepRVAT_testsplit0.25_r2.pq')

# Fix this
ols_model <- c('lm_sign_genes')
model <- c('funcrvp_omics_pops')

dt_r2_rel <- dcast(comp_r2[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')
dt_r2_rel <- melt(dt_r2_rel, id.vars = c("trait", ols_model), variable.name = "model", value.name = "model_r2")
dt_r2_rel <- melt(dt_r2_rel, id.vars = c("trait", "model", "model_r2"), variable.name = "baseline_model", value.name = "baseline_model_r2")

# Bootstraps file
r2_bt_multi <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_models_phenopred_10k_bootstrap_rvat.pq') %>% as.data.table()
r2_bt_multi <- r2_bt_multi[model %in% c('funcrvp', paste0('rvat_', p_thresh))]
r2_bt_multi[, rel_delta_r2 := (r2-cov_r2)/cov_r2]

base_dt <- r2_bt_multi[model==paste0('rvat_', p_thresh), .(bootstrap_iteration, trait, r2, rel_delta_r2)]
# base_dt <- r2_bt_multi[model==paste0('lm_',p_thresh), .(bootstrap_iteration, trait, r2, rel_delta_r2)]
colnames(base_dt) <- c('trial', 'trait', 'baseline_r2', 'rel_delta_r2_baseline')
base_dt[, baseline_model := 'lm_sign_genes']

model_dt <- r2_bt_multi[model=='funcrvp', .(bootstrap_iteration, model, trait, r2, rel_delta_r2)]
model_dt[, model:='funcrvp_omics_pops']
colnames(model_dt) <- c('trial', 'model', 'trait', 'model_r2','rel_delta_r2_model')

b_model_dt <- merge(model_dt, base_dt, by = c('trial', 'trait'))
b_model_dt[, r2_diff := model_r2-baseline_r2]

N_rep <- b_model_dt[, uniqueN(trial)]
N_trait <- b_model_dt[, uniqueN(trait)]
stats_plot <- b_model_dt[, .(N_greater = sum(r2_diff > 0), N_lesser = sum(r2_diff < 0)), by=.(model, trait, baseline_model)]
stats_plot[, pval:= 2 * pmin((N_greater + 1)/(N_rep + 1), (N_lesser + 1)/(N_rep + 1))]
stats_plot[, pval_fdr := p.adjust(stats_plot$pval, method = 'BY')]
# stats_plot[, pval_bonf:= pmin(1, pval*N_trait)]

dt_plot_r2_rel <- merge(dt_r2_rel, stats_plot, by = c("model", "trait", "baseline_model"))
dt_plot_r2_rel[, point_color := ifelse(pval_fdr<=0.05, ifelse(model_r2>=baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]

# drop_list <- c('Calcium', 'Cholesterol', 'Creatinine', 'Cystatin_C', 'HDL_cholesterol', 'High_light_scatter_reticulocyte_count', 'High_light_scatter_reticulocyte_percentage', 'LDL_direct', 'Lymphocyte_percentage', 'Mean_corpuscular_haemoglobin', 'Mean_corpuscular_volume', 'Mean_sphered_cell_volume', 'Platelet_crit','Reticulocyte_count','Total_bilirubin')
# dt_plot_r2_rel <- dt_plot_r2_rel[!(trait %in% drop_list)]

a <- rbindlist(list(dt_plot_r2_rel[, .(model, trait, model_r2)], setNames(dt_plot_r2_rel[, .(baseline_model, trait, baseline_model_r2)], names(dt_plot_r2_rel[, .(model, trait, model_r2)]))))
r2_pval <- wilcox_test(a, model_r2 ~ model, paired = T) %>% pull(p)

cols <- palette.colors(8, palette = "Dark2")
color_vec <- c('Significantly\nbetter' = cols[[1]], 'No significant\ndifference' = "gray", 'Significantly\nworse' = cols[[2]])
dt_plot_r2_rel[['point_color']] <- factor(dt_plot_r2_rel[['point_color']], levels = c('Significantly\nworse', 'No significant\ndifference', 'Significantly\nbetter'))

lm_plot <- ggplot(data = dt_plot_r2_rel, 
                  aes(x=baseline_model_r2, y=model_r2, color=point_color)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point(size=2) +
  geom_text(aes(x=0.04, y=0.005), label=paste0('Wilcoxon p-value\n',r2_pval), check_overlap = T, color='black', size=3.5) +
  geom_text_repel(data = dt_plot_r2_rel[pval<0.05],
                  aes(label = trait,
                      segment.size = 0.25,
                      segment.linetype=3),
                  # max.overlaps=0,
                  min.segment.length = 0,
                  box.padding = 1.25,
                  point.padding = 0,
                  force = 2,
                  size = 3.5,
                  show_guide  = FALSE,
                  arrow = arrow(length = unit(0.0015, "npc"))) +
  xlab(TeX("Linear model on significant genes, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  scale_y_sqrt(na.value=0) +
  scale_x_sqrt(na.value=0) +
  theme_cowplot() +
  scale_color_manual(name = "Significant\n(FDR<5%)",
                     values = color_vec,
                     labels = c('Significantly\nbetter' = paste0('Significantly\nbetter (', dt_plot_r2_rel[point_color=='Significantly\nbetter', .N],')'),
                                'No significant\ndifference' = paste0('No significant\ndifference (', dt_plot_r2_rel[point_color=='No significant\ndifference', .N],')'),
                                'Significantly\nworse' = paste0('Significantly\nworse (', dt_plot_r2_rel[point_color=='Significantly\nworse', .N],')'))) +
  guides(color = guide_legend(textsize=9, keyheight = 2.5)) +
  theme(legend.direction = "vertical",
        # legend.spacing.y = unit(2.0, 'cm'),
        legend.position = c(0.05, 0.8))

lm_plot
