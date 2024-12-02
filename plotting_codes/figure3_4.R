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

# FuncRVAT predictions
# dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/v108cov_deepRVAT_predictions_extended.pq') %>% as.data.table()
# LM model predictions
# dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/lm_filteredv3_deepRVAT_predictions_extended.pq') %>% as.data.table()
# dt_ols[, `:=` (trait_measurement = measurement,
#                pred = lm_pred,
#                measurement = NULL,
#                lm_pred = NULL)]


test_size <- '0.25'
version <- 'v1NEWsplit'
genotype <- 'deepRVAT'
emb <- 'omics_pops'

dt_bayes <- read_parquet(paste0('/s/project/geno2pheno/predictions/bayesian/fixed_arch/', version, '_', genotype, '_testsplit', test_size, '_', emb, '_predictions_extended.pq')) %>% as.data.table()
# dt_bayes <- read_parquet('/s/project/geno2pheno/predictions/bayesian/fixed_arch/v1NEWsplit_noemb_deepRVAT_testsplit0.25_noemb_predictions_extended.pq') %>% as.data.table()
dt_bayes[, pred := best_r2_pred]
dt_bayes[, best_r2_pred := NULL]
dt_bayes[, best_loss_pred := NULL]
dt_bayes[, .N, by=.(model, trait)]

dt_ols <- read_parquet(paste0('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_', genotype, '_', test_size,'_predictions_NEWsplit.pq')) %>% as.data.table()
# dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_am_plof_0.25_predictions_NEWsplit.pq') %>% as.data.table()
# dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/lm_filteredv3_pLoF_0.25_predictions_NEWsplit.pq') %>% as.data.table()

dt_ols[, .N, by=.(model, trait)]

dt <- rbindlist(list(dt_bayes, dt_ols), use.names=TRUE, fill = TRUE)
dt

dt[, .N, by=.(model, trait)]


bayes_model <- dt_bayes[, unique(model)][[1]]
ols_model <- dt_ols[, unique(model)][dt_ols[, unique(model)] != "lm_cov"]

dt[, version:=NULL]
dt[, genotype:=NULL]
dt[, common_residual:=NULL]

dt[model=='lm_cov', model_version := 'PRS_model_pred']
dt[model==ols_model, model_version := 'baseline_model_pred']
dt[model==bayes_model, model_version := 'model_pred']
dt_plot <- dcast(dt[, .(individual, trait, phenocode, trait_measurement, model_version, pred)], ... ~ model_version, value.var = 'pred')

dt_plot[, `:=` (PRS_model = 'lm_cov',
                model = bayes_model,
                baseline_model = ols_model)]

# Bootstrap -------------------------------------------------
N_rep <- 1000
list_rep <- lapply(1:N_rep, function(x) dt_plot[sample(1:.N, .N, replace=TRUE), .SD, by = .(model, trait, baseline_model)][,.(model_r2 = r2_general(model_pred, trait_measurement),
                                                                                                                              baseline_model_r2 = r2_general(baseline_model_pred, trait_measurement),
                                                                                                                              PRS_model_r2 = r2_general(PRS_model_pred, trait_measurement)), by=.(model, trait, baseline_model)])
names(list_rep) <- 1:length(list_rep)

bs_model_dt <- rbindlist(list_rep, idcol = "trial")
bs_model_dt[, `:=` (r2_diff = model_r2 - baseline_model_r2,
                    rel_delta_r2_model = model_r2 - PRS_model_r2,
                    rel_delta_r2_baseline = baseline_model_r2 - PRS_model_r2)]


std_err_dt <- bs_model_dt[,.(SE_rel_delta_r2_model = sd(rel_delta_r2_model),
               SE_rel_delta_r2_baseline = sd(rel_delta_r2_baseline)), by = .(model, baseline_model, trait)]

stats_plot <- bs_model_dt[, .(N_greater = sum(r2_diff > 0), N_lesser = sum(r2_diff < 0)), by=.(model, trait, baseline_model)]
stats_plot[, pval:= 2 * pmin((N_greater + 1)/(N_rep + 1), (N_lesser + 1)/(N_rep + 1))]

stats_plot <- merge(stats_plot, std_err_dt, by = c('model', 'baseline_model', 'trait'))
stats_plot

# # compute R2 for plot --------------------
a <- dcast(dt[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
a <- melt(a, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
a[, delta_r2 := r2 - lm_cov]
a[, rel_delta_r2 := delta_r2/lm_cov]

# # Relative Delta R2
dt_plot_r2_rel <- dcast(a[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')
dt_plot_r2_rel <- melt(dt_plot_r2_rel, id.vars = c("trait", ols_model), variable.name = "model", value.name = "model_r2")
dt_plot_r2_rel <- melt(dt_plot_r2_rel, id.vars = c("trait", "model", "model_r2"), variable.name = "baseline_model", value.name = "baseline_model_r2")

dt_plot_r2_rel <- merge(dt_plot_r2_rel, stats_plot[baseline_model != 'lm_cov'], by = c("model", "trait", "baseline_model"))

dt_plot_r2_rel[, unique(model)]

dt_plot_r2_rel[, point_color := ifelse(pval<=0.05, ifelse(model_r2>=baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]

dt_plot_r2_rel[, `:=` (x_min = baseline_model_r2 - 1.96*SE_rel_delta_r2_baseline,
                       x_max = baseline_model_r2 + 1.96*SE_rel_delta_r2_baseline,
                       y_min = model_r2 - 1.96*SE_rel_delta_r2_model,
                       y_max = model_r2 + 1.96*SE_rel_delta_r2_model)]


# fwrite(dt_plot_r2_rel, paste0('/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_', version, '_', emb, '_', genotype, '.csv'))
# fwrite(dt_plot_r2_rel, "/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_enformer_small_deepRVAT.csv")
# fwrite(dt_plot_r2_rel, "/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_am_plof.csv")
# fwrite(dt_plot_r2_rel, "/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_noemb_deepRVAT.csv")

# dt_plot_r2_rel <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_enformer_small_deepRVAT.csv")
# dt_plot_r2_rel <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_am_plof.csv")
# dt_plot_r2_rel <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_deepRVAT.csv")

a <- rbindlist(list(dt_plot_r2_rel[, .(model, trait, model_r2)], setNames(dt_plot_r2_rel[, .(baseline_model, trait, baseline_model_r2)], names(dt_plot_r2_rel[, .(model, trait, model_r2)]))))
r2_pval <- wilcox_test(a, model_r2 ~ model, paired = T) %>% pull(p)

cols <- palette.colors(8, palette = "Dark2")
color_vec <- c('Significantly\nbetter' = cols[[1]], 'No significant\ndifference' = "gray", 'Significantly\nworse' = cols[[2]])
# color_vec <- c('Significantly\nbetter' = "chartreuse4", 'No significant\ndifference' = "gray", 'Significantly\nworse'="firebrick3")
dt_plot_r2_rel[['point_color']] <- factor(dt_plot_r2_rel[['point_color']], levels = c('Significantly\nworse', 'No significant\ndifference', 'Significantly\nbetter'))


lm_plot <- ggplot(data = dt_plot_r2_rel, 
                  aes(x=baseline_model_r2, y=model_r2, color=point_color)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point(size=2) +
  geom_text(aes(x=0.04, y=0.005), label=paste0('Wilcoxon p-value\n',r2_pval), check_overlap = T, color='black', size=3) +
  geom_text_repel(data = dt_plot_r2_rel[pval<0.05],
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
  xlab(TeX("Linear model on significant genes, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  ylab(TeX("FuncRVP, DeepRVAT GIS (Relative $\\Delta R^2$)")) +
  # xlab(TeX("Linear model on significant genes, (Alphamissense + pLoF) burden (Relative $\\Delta R^2$)")) +
  # ylab(TeX("FuncRVP, (Alphamissense + pLoF) burden (Relative $\\Delta R^2$)")) +
  scale_y_sqrt(na.value=0) +
  scale_x_sqrt(na.value=0) +
  theme_cowplot() +
  scale_color_manual(name = "Significant (p < 0.05)",
                     values = color_vec,
                     labels = c('Significantly\nbetter' = paste0('Significantly\nbetter (', dt_plot_r2_rel[point_color=='Significantly\nbetter', .N],')'),
                                'No significant\ndifference' = paste0('No significant\ndifference (', dt_plot_r2_rel[point_color=='No significant\ndifference', .N],')'),
                                'Significantly\nworse' = paste0('Significantly\nworse (', dt_plot_r2_rel[point_color=='Significantly\nworse', .N],')'))) +
  guides(color = guide_legend(textsize=5, keyheight = 2.5)) +
  theme_pub +
  theme(legend.direction = "vertical",
        # legend.spacing.y = unit(2.0, 'cm'),
        legend.position = c(0.05, 0.8))

lm_plot

# Supplementary plots
library(readxl)

burden_her_file = "/s/project/geno2pheno/data/burden_heritability.xlsx"

phenotype_key_df <- read_excel(burden_her_file, sheet = "ST5") %>% as.data.table()
burden_herit <- read_excel(burden_her_file, sheet = "ST8") %>% as.data.table()

herit_dt <- merge(phenotype_key_df[, .(phenotype_key, phenocode)], burden_herit[, .(phenotype_key, aggregated_h2, aggregate_h2_se)], by = 'phenotype_key')
herit_dt <- herit_dt %>% mutate(phenocode = as.character(phenocode))

asd <- unique(dt[, .(trait, phenocode)])
asd <- asd %>% mutate(phenocode = as.character(phenocode))

asd <- merge(herit_dt, asd, by = "phenocode")

asd <- merge(dt_plot_r2_rel, asd, by = "trait") #, all.x = TRUE)
asd[, delta_imp := model_r2 - baseline_model_r2]
asd[, relative_imp := delta_imp/baseline_model_r2]

label_herit <- asd[baseline_model=='lm_sign_genes_deepRVAT', .N, by=point_color]
label_herit[, labs := paste0("n=", N)]

pv_herit <- asd %>%
  wilcox_test(aggregated_h2 ~ point_color) %>%
  add_significance("p") %>%
  filter(p.signif != 'ns') %>%
  add_y_position()
pv_herit

herit <- ggplot(asd, aes(x=point_color, y=aggregated_h2)) +
  geom_boxplot(aes(fill=point_color), width=0.5, alpha=0.4) +
  geom_dotplot(aes(fill=point_color), binaxis='y', stackdir='center', dotsize=0.5, alpha=0.8) +
  geom_text(data=label_herit, aes(y=-0.001, label=labs), size=3) + #y=-0.005, hjust=-0.01) +
  ylab('Burden heritability') +
  scale_fill_manual(values = color_vec) +
  theme_cowplot() +
  # guides(x = guide_axis(angle = 45)) +
  ylim(-0.005, 0.04) +
  theme_pub +
  theme(legend.position = 'none') +
  # stat_pvalue_manual(pv_herit, label = "p.signif", tip.length = 0, coord.flip = TRUE) +
  coord_flip() +
  xlab('')

herit

  
# Number of gene associations in three categories
cat_dt <- dt_plot_r2_rel[, list(trait, category = ifelse(pval<0.05, ifelse(model_r2>baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference'))]

genes_dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v2cleansplitHO_deepRVAT_omics_pops_genes_CLEANsplit.pq') %>% as.data.table()
genes_dt <- genes_dt[replicated==TRUE, .N, by=trait]

cat_dt <- merge(cat_dt, genes_dt, by='trait')
label_dt <- cat_dt[, .N, by=category]
label_dt[, labs := paste0("n=", N)]

pv_dt <- cat_dt %>%
  wilcox_test(N ~ category) %>%
  add_significance("p")
  # filter(p.signif != 'ns')
  # add_y_position()
pv_dt

rare <- ggplot(cat_dt, aes(x=reorder(category, N, median), y=N)) +
  geom_boxplot(aes(fill=category), width=0.5, alpha=0.4) +
  geom_dotplot(aes(fill=category), binaxis='y', stackdir='center', dotsize=0.5, alpha=0.8) +
  geom_text(data=label_dt, aes(y=-38, label=labs), size=3, hjust=-0.0) +
  stat_pvalue_manual(pv_dt, label = "p.signif", tip.length = 0, y.position = c(175), coord.flip = TRUE) +
  scale_fill_manual(values = color_vec) +
  ylab('Gene associations') +
  xlab('') +
  # guides(x = guide_axis(angle = 45)) +
  theme_cowplot() +
  coord_flip() +
  theme_pub +
  theme(legend.position = 'none')
  
rare

# Number of GWAS index variants
gwas_dt <- read_parquet('/s/project/geno2pheno/data/enrichment_data/index_variants.pq') %>% as.data.table()
peak_cts <- gwas_dt[, .(peaks = .N), by=trait]
cat_dt <- merge(cat_dt, peak_cts, by='trait')

label_dt <- cat_dt[, .N, by=category]
label_dt[, labs := paste0("n=", N)]

pv_dt2 <- cat_dt %>%
  wilcox_test(peaks ~ category) %>%
  add_significance("p") %>%
  filter(p.signif != 'ns') 
pv_dt2

common <- ggplot() +
  # geom_boxplot() +
  geom_boxplot(data=cat_dt, aes(x=reorder(category, peaks, median), y=peaks, fill=category), width=0.5, alpha=0.4) +
  geom_dotplot(data=cat_dt, aes(x=reorder(category, peaks, median), y=peaks, fill=category), 
               binaxis='y', stackdir='center', dotsize=0.25, alpha=0.8) +
  geom_text(data=label_dt, aes(x=category, y=400, label=labs), size=3, hjust=-0.01) +
  ylab('GWAS index variants') +
  xlab('') +
  scale_y_log10() +
  scale_fill_manual(values = color_vec) +
  annotation_logticks(sides='b') +
  guides(x = guide_axis(angle = 45)) +
  stat_pvalue_manual(pv_dt2, label = "p.signif", tip.length = 0, y.position = c(4.5), coord.flip = TRUE) +
  theme_cowplot() +
  coord_flip() +
  theme_pub +
  theme(legend.position = 'none')

common

# ggarrange(herit, rare, common, nrow=1)

# ggarrange(lm_plot, ggarrange(herit, rare, common, nrow=1, labels = c("B", "C", "D")), ncol=1, labels = c("A", ""), heights = c(4,3))

ggsave('/s/project/geno2pheno/figures/resub_figures/figure3.png', fig3, width=18, height=13, units = "cm", bg = 'white')
# ggsave('/s/project/geno2pheno/figures/paper_figures/figure3.png', fig3, width=18, height=13, units = "cm", bg = 'white')

# Embedding comparison
om_pop <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_deepRVAT.csv")
om_pop[, emb:='Omics+PoPS']
enf <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_enformer_small_deepRVAT.csv")
enf[, emb:='Enformer\n(ns)']
g2v <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_gene2vec_deepRVAT.csv")
g2v[, emb:='gene2vec\n(ns)']
esm <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_esm2_deepRVAT.csv")
esm[, emb:='ESM2*']

plof <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_pLoF.csv")
plof[, `:=` (emb='Omics+PoPS*', genotype='pLoF')]
am_plof <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_am_plof.csv")
am_plof[, `:=` (emb='Omics+PoPS\n(+)', genotype='AlphaMissense\n+ pLoF')]

drvt <- rbindlist(list(om_pop, enf, esm, g2v), use.names = T)
drvt[, genotype := 'DeepRVAT']

comp_dt <- rbindlist(list(drvt, plof), use.names = T)
comp_dt[, .N, by=.(genotype, emb)]

plt <- comp_dt[, .N, by=.(genotype, emb, point_color)]
plt <- rbindlist(list(plt, list('DeepRVAT', 'Omics+PoPS', 'Significantly\nworse', 0)))

plt[['genotype']] <- factor(plt[['genotype']], levels = c('DeepRVAT', 'AlphaMissense\n+ pLoF', 'pLoF'))
plt[['point_color']] <- factor(plt[['point_color']], levels = c('Significantly\nworse', 'No significant\ndifference', 'Significantly\nbetter'))
plt[['emb']] <- factor(plt[['emb']], levels = c('Omics+PoPS', 'Enformer\n(ns)', 'gene2vec\n(ns)', 'ESM2*', 'Omics+PoPS\n(+)', 'Omics+PoPS*'))
plt

cols <- wes_palette("AsteroidCity1")
cols2 <- wes_palette("Zissou1", 10, type = "continuous")
fill_vec <- c('Significantly\nbetter' = cols[[4]], 'No significant\ndifference' = "gray", 'Significantly\nworse' = cols2[[7]])

emb_comp <- ggplot(plt, aes(x=emb, y=N, fill=point_color)) +
  geom_col(alpha=0.85, width=0.7, color='black', position = 'dodge') +
  xlab('\nEmbedding') +
  ylab('Traits') +
  scale_fill_manual(name='Phenotype\nprediction', values = fill_vec) +
  facet_grid(~genotype, scales = 'free_x', space='free') +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  theme_cowplot() +
  theme_pub +
  guides(fill = guide_legend(textsize=3, byrow = TRUE)) +  #, keyheight=1)) +
  theme(legend.position='right',
        panel.grid.major.y = element_line(colour="gray75", size=0.5),
        legend.spacing.y = unit(1, "cm"))
  
emb_comp

ggsave('/s/project/geno2pheno/figures/resub_figures/figure4_corr.png', emb_comp, width=15, height=10, units = "cm", bg = 'white')
# ggsave('/s/project/geno2pheno/figures/resub_figures/figure4.png', emb_comp, width=18, height=10, units = "cm", bg = 'white')


# fig3 <- ggarrange(lm_plot, ggarrange(rare, common, herit, ncol=1, labels = c("B", "C", "D"), font.label = list(size = 12)), ncol=2, labels = c("A", ""), widths = c(4,3), font.label = list(size = 12))
# fig3_ext <- ggarrange(fig3, emb_comp, nrow = 2, labels = c("", "E"), heights = c(1.75,1))
# ggsave('/s/project/geno2pheno/figures/resub_figures/figure3_ext.png', fig3_ext, width=18, height=19, units = "cm", bg = 'white')
