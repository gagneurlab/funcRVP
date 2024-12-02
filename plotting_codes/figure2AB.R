library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(arrow)
library(rstatix)
library(ggpubr)
library(latex2exp)
library(dplyr)

# Compute R2 -------------------------------------------------    
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

theme_pub <- theme(axis.text=element_text(size=9),
                   text = element_text(size=9))

# Other embeddings
dt_mean <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v5cov_deepRVAT_mean_predictions_extended.pq') %>% as.data.table()
dt_mean[, .(r2=r2_general(pred, trait_measurement)), by=trait]

ggplot(dt_mean[trait %in% c('Vitamin_D', 'Glucose')]) +
  geom_hex(aes(x=pred, y=trait_measurement)) +
  facet_grid(~trait) +
  theme_minimal()

# Omics_pops
dt_bayes <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v108cov_deepRVAT_predictions_extended.pq') %>% as.data.table()

# dt_plof <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v98cov_pLoF_predictions_extended.pq') %>% as.data.table()
# dt_plof <- dt_plof[model=='omics_pops_bayesian_v98cov_pLoF']

# dt_flat <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v107cov_flat_deepRVAT_predictions_extended.pq') %>% as.data.table()
# dt_flat

dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_deepRVAT_predictions_extended.pq') %>% as.data.table()

dt <- rbindlist(list(dt_mean, dt_bayes, dt_ols), use.names = TRUE, fill = TRUE)
dt[, `:=` (version = NULL, phenocode = NULL)]

a <- dcast(dt[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
a <- melt(a, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
a[, delta_r2 := r2 - lm_cov]
a[, rel_delta_r2 := delta_r2/lm_cov]

a[order(r2), by=model]

ap_val <- a[model %in% c('omics_pops_bayesian_v5cov_deepRVAT', 'omics_pops_bayesian_v108cov_deepRVAT')] %>% wilcox_test(rel_delta_r2 ~ model, paired = T)

# model_plot_dt <- dcast(a[,.(trait, model, r2)], ... ~ model, value.var = 'r2')
model_plot_dt <- dcast(a[,.(trait, model, r2)], ... ~ model, value.var = 'r2')

model_plot <- ggplot(data = model_plot_dt, aes(x=omics_pops_bayesian_v5cov_deepRVAT, y=omics_pops_bayesian_v108cov_deepRVAT)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  geom_text_repel(data = model_plot_dt,
                  aes(label = trait,
                      segment.size = 0.25,
                      segment.linetype=3),
                  # max.overlaps=0,
                  min.segment.length = 0,
                  box.padding = 1.25,
                  point.padding = 0,
                  force = 2,
                  size = 3,
                  show_guide  = FALSE,
                  arrow = arrow(length = unit(0.0015, "npc"))) +
  geom_text(aes(x=0.25, y=-0.75), label=paste0('Wilcoxon p-value\n', ap_val$p), check_overlap = T, size=3) +
  xlim(-1.2, 0.9) +
  ylim(-1.2, 0.9) +
  # xlab(TeX("FuncRVP modeling expectation (Relative $\\Delta R^2$)")) +
  # ylab(TeX("FuncRVP modeling variance (Relative $\\Delta R^2$)")) +
  xlab(TeX("$R^2$, FuncRVP modeling expectation")) +
  ylab(TeX("$R^2$, FuncRVP modeling variance")) +
  theme_cowplot() +
  theme_pub

model_plot
ggsave('/s/project/geno2pheno/figures/resub_figures/figure2A.png', model_plot, width=10, height=10, units = "cm", bg = 'white')


# Emb vs Permuted embeddings
dt_emb <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_perm <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_noemb_deepRVAT_testsplit0.25_noemb_predictions_extended.pq') %>% as.data.table()
# dt_perm <-  read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_shuffEmb_shuffledemb0_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
dt_lm_cov <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_deepRVAT_0.25_predictions_NEWsplit.pq') %>% as.data.table()
dt_lm <- dt_lm_cov[model=='lm_sign_genes_deepRVAT']
dt_lm_cov <- dt_lm_cov[model=='lm_cov']

bayes <- rbindlist(list(dt_emb, dt_perm, dt_lm, dt_lm_cov), use.names = TRUE, fill = TRUE)
bayes[is.na(pred), pred := best_r2_pred]
bayes[, `:=` (version=NULL, phenocode=NULL, genotype=NULL, best_loss_pred=NULL, best_r2_pred=NULL)]
bayes

b <- dcast(bayes[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
b <- melt(b, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
b[, delta_r2 := r2 - lm_cov]
b[, rel_delta_r2 := delta_r2/lm_cov]
emb_plot_dt <- dcast(b[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')

bp_val <- b %>% wilcox_test(rel_delta_r2 ~ model, paired = T) 
bp_val

emb_plot <- ggplot(data = emb_plot_dt, aes(x=lm_sign_genes_deepRVAT, y=omics_pops_bayesian_v1NEWsplit_noEmb_noemb_deepRVAT)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  # geom_text(aes(x=0.05, y=0.005), label=paste0('Wilcoxon p-value\n',bp_val), check_overlap = T, size=3) +
  xlab(TeX("LM (Relative $\\Delta R^2$)")) + 
  ylab(TeX("FuncRVP, No embeddings (Relative $\\Delta R^2$)")) +
  # xlim(b[, min(rel_delta_r2)], r2[, max(rel_delta_r2)]) +
  # ylim(b[, min(rel_delta_r2)], r2[, max(rel_delta_r2)]) +
  # scale_y_sqrt(na.value=0, limits=c(0, 0.09)) +
  # scale_x_sqrt(na.value=0, limits=c(0, 0.09)) +
  theme_cowplot() +
  theme_pub

emb_plot

fig_2ab <- ggarrange(model_plot, emb_plot, nrow=2, labels = c("A", "B"))
fig_2ab

ggsave('/s/project/geno2pheno/figures/resub_figures/figure2AB.png', fig_2ab, width=9, height=19, units = "cm", bg = 'white')

# OLD plot -------------------------------!!!!!!!!!!!!!!!!!!!!---------------------------------
# bayes_plof <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_pLoF_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
# bayes_dprvt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_predictions_extended.pq') %>% as.data.table()
# dt_lm_cov <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_pLoF_0.25_predictions_NEWsplit.pq') %>% as.data.table()
# dt_lm_cov <- dt_lm_cov[model=='lm_cov']
# 
# bayes <- rbindlist(list(bayes_dprvt, bayes_plof, dt_lm_cov), use.names = TRUE, fill = TRUE)
# bayes[is.na(pred), pred := best_r2_pred]
# bayes[, `:=` (version=NULL, phenocode=NULL, genotype=NULL, best_loss_pred=NULL, best_r2_pred=NULL)]
# 
# b <- dcast(bayes[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
# b <- melt(b, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
# b[, delta_r2 := r2 - lm_cov]
# b[, rel_delta_r2 := delta_r2/lm_cov]
# geno_plot_dt <- dcast(b[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')
# 
# geno_plot <- ggplot(data = geno_plot_dt, aes(x=omics_pops_bayesian_v1NEWsplit_pLoF, y=omics_pops_bayesian_v1NEWsplit_deepRVAT)) +
#   geom_abline(color="gray", linetype="dashed") +
#   geom_point() +
#   xlab(TeX("FuncRVP using pLoF (Relative $\\Delta R^2$)")) + 
#   ylab(TeX("FuncRVP using DeepRVAT (Relative $\\Delta R^2$)")) + 
#   scale_y_sqrt(na.value=0, limits=c(0, 0.09)) +
#   scale_x_sqrt(na.value=0, limits=c(0, 0.09)) +
#   theme_cowplot() +
#   guides(color = guide_legend(textsize=5, keyheight = 2.5))
# 
# geno_plot
# 
# ggarrange(model_plot, geno_plot, nrow=1, labels = c("A", "C"))





# reg_plot <- ggplot(data = geno_plot_dt, aes(x=omics_pops_bayesian_v107cov_flat_deepRVAT, y=omics_pops_bayesian_v108cov_deepRVAT)) +
#   geom_abline(color="gray", linetype="dashed") +
#   geom_point() +
#   xlab(TeX("Bayesian model with constant prior (Relative $\\Delta R^2$)")) + 
#   ylab(TeX("FuncRVP using embeddings (Relative $\\Delta R^2$)")) + 
#   scale_y_sqrt(na.value=0, limits=c(0, 0.1)) +
#   scale_x_sqrt(na.value=0, limits=c(0, 0.1)) +
#   theme_cowplot() +
#   guides(color = guide_legend(textsize=5, keyheight = 2.5))
# 
# reg_plot
