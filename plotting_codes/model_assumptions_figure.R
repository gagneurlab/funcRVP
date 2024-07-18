library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(arrow)
library(rstatix)
library(ggpubr)
library(latex2exp)


# Other embeddings
dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v5cov_deepRVAT_mean_predictions_extended.pq') %>% as.data.table()
dt[, unique(model)]

# Omics_pops
dt2 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v108cov_deepRVAT_predictions_extended.pq') %>% as.data.table()

dt_plof <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v98cov_pLoF_predictions_extended.pq') %>% as.data.table()
dt_plof <- dt_plof[model=='omics_pops_bayesian_v98cov_pLoF']

dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_deepRVAT_predictions_extended.pq') %>% as.data.table()

dt_flat <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v107cov_flat_deepRVAT_predictions_extended.pq') %>% as.data.table()
dt_flat

dt <- rbindlist(list(dt, dt2, dt_ols, dt_plof, dt_flat), use.names = TRUE, fill = TRUE)

dt[, unique(model)]
dt[, uniqueN(trait), by=.(model)]
dt[, `:=` (version = NULL, phenocode = NULL)]


# Compute R2 -------------------------------------------------    
r2_general <-function(preds,actual){ 
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

a <- dcast(dt[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
a <- melt(a, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
a[, delta_r2 := r2 - lm_cov]
a[, rel_delta_r2 := delta_r2/lm_cov]
a

a[, uniqueN(trait)]

model_plot_dt <- dcast(a[,.(trait, model, r2)], ... ~ model, value.var = 'r2')

model_plot <- ggplot(data = model_plot_dt, aes(x=omics_pops_bayesian_v5cov_deepRVAT, y=omics_pops_bayesian_v108cov_deepRVAT)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  xlim(-1.2, 0.85) +
  ylim(-1.2, 0.85) +
  xlab(TeX("$R^2$, FuncRVP modeling expectation")) +
  ylab(TeX("$R^2$, FuncRVP modeling variance")) +
  # scale_y_sqrt(na.value=0) +
  # scale_x_sqrt(na.value=0) +
  theme_cowplot() +
  guides(color = guide_legend(textsize=5, keyheight = 2.5))

model_plot


geno_plot_dt <- dcast(a[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')

geno_plot <- ggplot(data = geno_plot_dt, aes(x=omics_pops_bayesian_v98cov_pLoF, y=omics_pops_bayesian_v108cov_deepRVAT)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  xlab(TeX("FuncRVP using pLoF (Relative $\\Delta R^2$)")) + 
  ylab(TeX("FuncRVP using DeepRVAT (Relative $\\Delta R^2$)")) + 
  # xlab(TeX("Relative $ \\Delta R^2_{FuncRVAT \\; using \\; pLoF} $")) +
  # ylab(TeX("Relative $ \\Delta R^2_{FuncRVAT \\; using \\; DeepRVAT} $")) +
  # xlab(TeX("$R^2$, FuncRVAT using pLoF")) +
  # ylab(TeX("$R^2$, FuncRVAT using DeepRVAT")) +
  # xlim(0, 0.1) +
  # ylim(0, 0.1) +
  scale_y_sqrt(na.value=0, limits=c(0, 0.1)) +
  scale_x_sqrt(na.value=0, limits=c(0, 0.1)) +
  theme_cowplot() +
  guides(color = guide_legend(textsize=5, keyheight = 2.5))

geno_plot

ggarrange(model_plot, geno_plot, nrow=1, labels = c("A", "C"))





reg_plot <- ggplot(data = geno_plot_dt, aes(x=omics_pops_bayesian_v107cov_flat_deepRVAT, y=omics_pops_bayesian_v108cov_deepRVAT)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  xlab(TeX("Bayesian model with constant prior (Relative $\\Delta R^2$)")) + 
  ylab(TeX("FuncRVP using embeddings (Relative $\\Delta R^2$)")) + 
  # xlab(TeX("Relative $ \\Delta R^2_{FuncRVAT \\; using \\; pLoF} $")) +
  # ylab(TeX("Relative $ \\Delta R^2_{FuncRVAT \\; using \\; DeepRVAT} $")) +
  # xlab(TeX("$R^2$, FuncRVAT using pLoF")) +
  # ylab(TeX("$R^2$, FuncRVAT using DeepRVAT")) +
  # xlim(0, 0.1) +
  # ylim(0, 0.1) +
  scale_y_sqrt(na.value=0, limits=c(0, 0.1)) +
  scale_x_sqrt(na.value=0, limits=c(0, 0.1)) +
  theme_cowplot() +
  guides(color = guide_legend(textsize=5, keyheight = 2.5))

reg_plot
