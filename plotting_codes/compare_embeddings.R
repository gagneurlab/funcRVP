library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(arrow)
library(rstatix)
library(ggpubr)
library(latex2exp)

dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_deepRVAT_predictions_extended.pq') %>% as.data.table()

# omics_pops
dt1 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v108cov_deepRVAT_predictions_extended.pq') %>% as.data.table()

# pops_exp
dt2 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v111arch_deepRVAT_pops_exp_predictions_extended.pq') %>% as.data.table()
dt2[, pred := best_r2_pred]

# pops
dt3 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v112arch_deepRVAT_pops_predictions_extended.pq') %>% as.data.table()
dt3[, pred := best_r2_pred]

# omics
dt4 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v113arch_deepRVAT_omics_predictions_extended.pq') %>% as.data.table()
dt4[, pred := best_r2_pred]

# omics_pops_exp
dt5 <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v114arch_deepRVAT_omics_pops_exp_predictions_extended.pq') %>% as.data.table()
dt5[, pred := best_r2_pred]

dt <- rbindlist(list(dt1, dt2, dt3, dt4, dt5, dt_ols), use.names = TRUE, fill = TRUE)
dt


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

a[, uniqueN(trait)]

model_name_vec <- c(
  "lm_sign_genes_deepRVAT" = "Linear model",
  "omics_pops_bayesian_v108cov_deepRVAT" = "Omics + PoPS",
  "pops_exp_bayesian_v111arch" = "PoPS Exp"
)

emb_color_vec <- c(
  "lm_sign_genes_deepRVAT" = "#FF7F00",
  "omics_pops_bayesian_v108cov_deepRVAT" = "#6A3D9A",
  "pops_exp_bayesian_v111arch_deepRVAT" = "ivory3",
  "pops_bayesian_v112arch_deepRVAT" = "ivory3",
  "omics_bayesian_v113arch_deepRVAT" = "ivory3",
  "omics_pops_exp_bayesian_v114arch_deepRVAT" = "ivory3"
)

a[, model_name := ifelse(model=='lm_sign_genes_deepRVAT', 'Linear model', 
                         ifelse(model=='omics_pops_bayesian_v108cov_deepRVAT', 'Omics + PoPS',
                                ifelse(model=='pops_exp_bayesian_v111arch_deepRVAT', 'PoPS Exp',
                                       ifelse(model=='pops_bayesian_v112arch_deepRVAT', 'PoPS',
                                              ifelse(model=='omics_bayesian_v113arch_deepRVAT', 'Omics',
                                                     ifelse(model=='omics_pops_exp_bayesian_v114arch_deepRVAT', 'Omics + PoPS Exp', 0))))))]


ggplot(a, aes(reorder(model_name,rel_delta_r2, median), rel_delta_r2, fill=model)) +
  geom_boxplot(alpha=0.6) +
  ylab(TeX("Relative $\\Delta R^2$")) +
  scale_fill_manual(name="", values=emb_color_vec) +
  xlab('Embedding used') +
  scale_y_sqrt(na.value=0) +
  guides(x = guide_axis(angle=45)) +
  theme_cowplot() +
  theme(legend.position = 'none')
