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

theme_pub <- theme(axis.text=element_text(size=9),
                   text = element_text(size=9))

cols <- palette.colors(8, palette = "Okabe-Ito")

op_base_dir <- '/s/project/geno2pheno/funcrvp/paper_revisions/predictions/'
geno_dir <- 'ukbb_wes_500k_DeepRVAT_final_090924_medianshifted/'
rvat_dir <- 'rvat/'
p_thresh <- c('0.001', '0.005', '0.01', '0.05', '0.1', '0.5', '1.0', '0.0001nom', '0.001nom',
'0.01nom', '0.05nom', '0.1nom', '0.5nom', '1.0nom')

dt_ols_list <- lapply(p_thresh, FUN = function(pval){
  read_parquet(paste0(op_base_dir, geno_dir, rvat_dir, 'all_traits_phenopred_', pval, '.pq')) %>% as.data.table()
})

names(dt_ols_list) <- p_thresh
dt_ols <- rbindlist(dt_ols_list, idcol = "p_thresh", fill=T)
dt_ols

dt_ols[, embedding := 'None']
dt_ols[, .N, by=.(model, trait)]

dt_cov <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_traits_covariates_only_phenopred_filteredv3.pq') %>% as.data.table()
dt_cov[, phenocode := NULL]
dt_cov[, p_thresh := 'lm_cov']

dt <- rbindlist(list(dt_ols, dt_cov), use.names=TRUE, fill = TRUE)
dt

# compute R2 for plot --------------------
a <- dcast(dt[, .(r2 = r2_general(best_prediction, trait_measurement)), by = .(p_thresh, trait)], ...~p_thresh, value.var = 'r2')
a <- melt(a, id.vars = c('trait', 'lm_cov'), variable.name = "p_thresh", value.name = "r2")
a[, delta_r2 := r2 - lm_cov]
a[, rel_delta_r2 := delta_r2/lm_cov]
a

label_dict_unsorted = c(
  '0.001' = '0.001 bonf.',
  '0.005' = '0.005 bonf.',
  '0.01' = '0.01 bonf.',
  '0.05' = '0.05 bonf.',
  '0.1' = '0.1 bonf.',
  '0.5' = '0.5 bonf.',
  '1.0' = '<1.0 bonf.\n(0.000055 nominal p)',
  '0.0001nom' = '0.0001',
  '0.001nom' = '0.001',
  '0.01nom' = '0.01',
  '0.05nom' = '0.05',
  '0.1nom' = '0.1',
  '0.5nom' = '0.5',
  '1.0nom' = '1.0'
)

label_dt <- data.table(p_thresh_key = names(label_dict_unsorted), pval_label = label_dict_unsorted)
plot_dt <- merge(a, label_dt, by.x = "p_thresh", by.y = "p_thresh_key", all.x = TRUE)

pv_dt <- plot_dt %>%
  wilcox_test(rel_delta_r2 ~ pval_label, paired=T) %>%
  add_significance("p") %>%
  filter(p.signif != 'ns') %>%
  add_y_position(step.increase = 0.005) %>%
  as.data.table()

plot_dt$pval_label <- factor(
  plot_dt$pval_label, 
  levels = c('0.001 bonf.', '0.005 bonf.', '0.01 bonf.', '0.05 bonf.', "0.1 bonf.", "0.5 bonf.", "<1.0 bonf.\n(0.000055 nominal p)", '0.0001', '0.0005', '0.001', '0.005', '0.01', '0.05', '0.1', '0.5', '1.0')
  )
plot_dt

fp <- ggplot(plot_dt, aes(x=pval_label, y=rel_delta_r2)) +
  geom_boxplot() +
  xlab("p-value thresholds") +
  ylab(TeX("LM on significant genes, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  guides(x = guide_axis(angle = 45)) +
  theme_bw()


rem_thresh <- c('0.001 bonf.', '0.001', '0.01', '0.05', '0.1', '0.5', '1.0')

sp <- ggplot(plot_dt[!(pval_label %in% rem_thresh)], aes(x=pval_label, y=rel_delta_r2)) +
  geom_boxplot() +
  xlab("p-value thresholds") +
  ylab(TeX("LM on significant genes, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  stat_pvalue_manual(pv_dt[!(group1 %in% rem_thresh) & !(group2 %in% rem_thresh)], label = "p.signif", tip.length = 0.01) +
  guides(x = guide_axis(angle = 45)) +
  theme_bw()


ggarrange(fp, sp, ncol=1, labels = c("A", "B"), font.label = list(size = 12))
