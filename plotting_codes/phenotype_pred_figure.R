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

model_name_vec <- c(
  "lm_sign_genes_deepRVAT" = "OLS",
  "omics_pops_bayesian_v107cov_flat_deepRVAT" = "Flat prior DeepRVAT",
  "omics_pops_bayesian_v108cov_deepRVAT" = "FuncRVAT"
)

# FuncRVAT predictions
dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/v108cov_deepRVAT_predictions_extended.pq') %>% as.data.table()
# dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v1_same_arch_deepRVAT_omics_pops_predictions_extended.pq') %>% as.data.table()


# LM model predictions
dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/lm_filteredv3_deepRVAT_predictions_extended.pq') %>% as.data.table()
dt_ols[, version := NA]
dt_ols[, .N, by=.(model, trait)]


dt <- rbindlist(list(dt, dt_ols), use.names=TRUE, fill = TRUE)


dt[, .N, by=.(model, trait)] 


model_list <- as.character(dt[, unique(model)])
baseline_models <- model_list[grepl('ols', model_list, ignore.case = TRUE, fixed = TRUE) | grepl('lm', model_list, ignore.case = TRUE, fixed = TRUE) | grepl('flat', model_list, ignore.case = TRUE, fixed = TRUE)]
# bayesian_models <- model_list[grepl('bayesian', model_list, ignore.case = TRUE, fixed = TRUE)]

dt[, version:=NULL]
dt[, genotype:=NULL]

# dt[, new_model_pred:=common_residual + pred]
# dt[, pred:=new_model_pred]     #NEW CHECK!!!!!!!!!!!!!!!!!
# dt[, new_model_pred:=NULL]

dt[, common_residual:=NULL]

# PREP data for bootstraps
dt_plot <- dcast(dt, ... ~ model, value.var = 'pred')
dt_plot <- melt(dt_plot, id.vars = c("trait", "phenocode", "individual", "trait_measurement", baseline_models), variable.name = "model", value.name = "model_pred")
dt_plot <- melt(dt_plot, id.vars = c("trait", "phenocode", "individual", "trait_measurement", "model", "model_pred"), variable.name = "baseline_model", value.name = "baseline_model_pred")
dt_plot[, unique(baseline_model)]

dt_plot <- dt_plot[baseline_model != 'lm_cov']
unique(dt_plot[, .(model, baseline_model)])


# Compare deepRVAT against deepRVAT and pLoF against pLoF
dt_plot <- dt_plot[(grepl('pLoF', model, ignore.case = TRUE, fixed = TRUE) & grepl('pLoF', baseline_model, ignore.case = TRUE, fixed = TRUE)) | (grepl('deepRVAT', model, ignore.case = TRUE, fixed = TRUE) & grepl('deepRVAT', baseline_model, ignore.case = TRUE, fixed = TRUE))]
unique(dt_plot[, .(model, baseline_model)])

dt_plot <- dcast(dt, ... ~ model, value.var = 'pred')
#
dt_plot[, `:=` (model = 'omics_pops_bayesian_v108cov_deepRVAT',
                model_pred = omics_pops_bayesian_v108cov_deepRVAT,
                baseline_model = 'lm_sign_genes_deepRVAT',
                baseline_model_pred = lm_sign_genes_deepRVAT,
                PRS_model = 'lm_cov',
                PRS_model_pred = lm_cov)]
#
dt_plot[, `:=` (omics_pops_bayesian_v108cov_deepRVAT = NULL,
                lm_sign_genes_deepRVAT = NULL,
                lm_cov = NULL)]

# # Compute R2 -------------------------------------------------
r2_general <-function(preds,actual){
  return(1- sum((preds - actual) ^ 2)/sum((actual - mean(actual))^2))
}

#
# # Bootstrap -------------------------------------------------
N_rep <- 2000
# list_rep <- lapply(1:N_rep, function(x) dt_plot[sample(1:.N, .N, replace=TRUE), .SD, by = .(model, trait, baseline_model)][,.(r2_diff = r2_general(model_pred, trait_measurement) - r2_general(baseline_model_pred, trait_measurement)), by=.(model, trait, baseline_model)])
list_rep <- lapply(1:N_rep, function(x) dt_plot[sample(1:.N, .N, replace=TRUE), .SD, by = .(model, trait, baseline_model)][,.(model_r2 = r2_general(model_pred, trait_measurement),
                                                                                                                              baseline_model_r2 = r2_general(baseline_model_pred, trait_measurement),
                                                                                                                              PRS_model_r2 = r2_general(PRS_model_pred, trait_measurement)), by=.(model, trait, baseline_model)])
names(list_rep) <- 1:length(list_rep)
#
bs_model_dt <- rbindlist(list_rep, idcol = "trial")
bs_model_dt[, `:=` (r2_diff = model_r2 - baseline_model_r2,
                    rel_delta_r2_model = model_r2 - PRS_model_r2,
                    rel_delta_r2_baseline = baseline_model_r2 - PRS_model_r2)]
#
std_err_dt <- bs_model_dt[,.(SE_rel_delta_r2_model = sd(rel_delta_r2_model),
               SE_rel_delta_r2_baseline = sd(rel_delta_r2_baseline)), by = .(model, baseline_model, trait)]
#
#
stats_plot <- bs_model_dt[, .(N_greater = sum(r2_diff > 0), N_lesser = sum(r2_diff < 0)), by=.(model, trait, baseline_model)]
stats_plot[, pval:= 2 * pmin((N_greater + 1)/(N_rep + 1), (N_lesser + 1)/(N_rep + 1))]
#
stats_plot <- merge(stats_plot, std_err_dt, by = c('model', 'baseline_model', 'trait'))
stats_plot
#
# # compute R2 for plot --------------------
a <- dcast(dt[, .(r2 = r2_general(pred, trait_measurement)), by = .(model, trait)], ...~model, value.var = 'r2')
a <- melt(a, id.vars = c('trait', 'lm_cov'), variable.name = "model", value.name = "r2")
a[, delta_r2 := r2 - lm_cov]
a[, rel_delta_r2 := delta_r2/lm_cov]

# # Relative Delta R2
dt_plot_r2_rel <- dcast(a[,.(trait, model, rel_delta_r2)], ... ~ model, value.var = 'rel_delta_r2')
dt_plot_r2_rel <- melt(dt_plot_r2_rel, id.vars = c("trait", baseline_models[baseline_models != 'lm_cov']), variable.name = "model", value.name = "model_r2")
dt_plot_r2_rel <- melt(dt_plot_r2_rel, id.vars = c("trait", "model", "model_r2"), variable.name = "baseline_model", value.name = "baseline_model_r2")
dt_plot_r2_rel  <- melt(dt_plot_r2_rel , id.vars = c("trait", 'lm_sign_genes_pLoF', 'lm_sign_genes_deepRVAT'), variable.name = "model", value.name = "model_r2")
dt_plot_r2_rel  <- melt(dt_plot_r2_rel , id.vars = c("trait", "model", "model_r2"), variable.name = "baseline_model", value.name = "baseline_model_r2")
#
dt_plot_r2_rel <- merge(dt_plot_r2_rel, stats_plot[baseline_model != 'lm_cov'], by = c("model", "trait", "baseline_model"))
#
dt_plot_r2_rel[, unique(model)]

dt_plot_r2_rel[model=='omics_pops_bayesian_v108cov_deepRVAT' & pval<0.05, .(trait, model_r2, baseline_model_r2, pval)]
#
dt_plot_r2_rel[, point_color := ifelse(pval<=0.05, ifelse(model_r2>=baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
dt_plot_r2_rel[point_color=='Significantly\nbetter' & baseline_model=='lm_sign_genes_deepRVAT', .N]
#
dt_plot_r2_rel[, `:=` (x_min = baseline_model_r2 - 1.96*SE_rel_delta_r2_baseline,
                       x_max = baseline_model_r2 + 1.96*SE_rel_delta_r2_baseline,
                       y_min = model_r2 - 1.96*SE_rel_delta_r2_model,
                       y_max = model_r2 + 1.96*SE_rel_delta_r2_model)]

# color_vec <- c('Significantly\nbetter' = "#33A02C", 'No significant\ndifference' = "gray", 'Significantly\nworse'="#FF5F30")
color_vec <- c('Significantly\nbetter' = "chartreuse4", 'No significant\ndifference' = "gray", 'Significantly\nworse'="firebrick3")

theme(axis.text=element_text(size=8),
      text = element_text(size=8))

# fwrite(dt_plot_r2_rel, "/s/project/geno2pheno/results/bootstrap_r2_results_with_SE.csv")
# dt_plot_r2_rel <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE.csv")

# 'omics_pops_bayesian_v1_same_arch_deepRVAT_omics_pops'
lm_plot <- ggplot(data = dt_plot_r2_rel[model=='omics_pops_bayesian_v108cov_deepRVAT' & baseline_model=='lm_sign_genes_deepRVAT'], aes(x=baseline_model_r2, y=model_r2, color=point_color)) +
  geom_abline(color="gray", linetype="dashed") +
  geom_point() +
  # geom_errorbar(aes(ymin=y_min, ymax=y_max, width = 0)) +
  # geom_pointrange(aes(ymin=y_min, ymax=y_max), fatten = 0.75) +
  # geom_errorbarh(aes(xmin=x_min, xmax=x_max, height = 0)) +
  geom_text_repel(data = dt_plot_r2_rel[baseline_model=='lm_sign_genes_deepRVAT' & pval<0.05],
                  aes(label = trait,
                      segment.size = 0.25,
                      segment.linetype=3),
                  # max.overlaps=0,
                  min.segment.length = 0,
                  box.padding = 1,
                  point.padding = 0,
                  force = 2,
                  size = 3,
                  show_guide  = FALSE,
                  arrow = arrow(length = unit(0.0015, "npc"))) + 
  xlab(TeX("Linear model on significant genes (Relative $\\Delta R^2$)")) + 
  ylab(TeX("FuncRVP (Relative $\\Delta R^2$)")) + 
  scale_y_sqrt(na.value=0) +
  scale_x_sqrt(na.value=0) +
  theme_cowplot() +
  scale_color_manual(name = "Significant (p < 0.05)",
                     values = color_vec,
                     labels = c('Significantly\nbetter' = paste0('Significantly\nbetter (', dt_plot_r2_rel[point_color=='Significantly\nbetter' & baseline_model=='lm_sign_genes_deepRVAT', .N],')'),
                                'No significant\ndifference' = paste0('No significant\ndifference (', dt_plot_r2_rel[point_color=='No significant\ndifference' & baseline_model=='lm_sign_genes_deepRVAT', .N],')'),
                                'Significantly\nworse' = paste0('Significantly\nworse (', dt_plot_r2_rel[point_color=='Significantly\nworse' & baseline_model=='lm_sign_genes_deepRVAT', .N],')'))) +
  guides(color = guide_legend(textsize=3, keyheight = 2.5)) +
  theme(legend.direction = "vertical",
        axis.text=element_text(size=9),
        text = element_text(size=9),
        # text = element_text(size=15),
        # legend.spacing.y = unit(2.0, 'cm'),
        legend.position = c(0.1, 0.8))

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

pv_herit <- asd %>%
  wilcox_test(aggregated_h2 ~ point_color, alternative = "less") %>%
  add_significance("p") %>%
  filter(p.signif != 'ns') %>%
  add_y_position()
pv_herit

label_dt <- asd[baseline_model=='lm_sign_genes_deepRVAT', .N, by=point_color]
label_dt[, labs := paste0("n=", N)]

herit <- ggplot(asd, aes(x=point_color, y=aggregated_h2)) +
  geom_boxplot(aes(fill=point_color), width=0.5, alpha=0.4) +
  geom_dotplot(aes(fill=point_color), binaxis='y', stackdir='center', dotsize=0.75, alpha=0.8) +
  geom_text(data=label_dt, aes(y=-0.005, label=labs), size=2.5, hjust=-0.01) +
  ylab('Burden heritability') +
  # scale_fill_manual(values = c('Significantly\nbetter' = "#33A02C", 'No significant\ndifference' = "gray")) +
  scale_fill_manual(values = color_vec) +
  theme_cowplot() +
  # guides(x = guide_axis(angle = 45)) +
  theme(legend.position = 'none',
        axis.text=element_text(size=9),
        text = element_text(size=9)) +
  stat_pvalue_manual(pv_herit, label = "p.signif", tip.length = 0, coord.flip = TRUE) +
  coord_flip() +
  xlab('')

herit

  
# Number of gene associations in three categories
cat_dt <- dt_plot_r2_rel[, list(trait, category = ifelse(pval<0.05, ifelse(model_r2>baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference'))]

genes_dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v108cov_deepRVAT_genes_extended.pq') %>% as.data.table()
genes_dt <- genes_dt[replicated==TRUE, .N, by=trait]

cat_dt <- merge(cat_dt, genes_dt, by='trait')
label_dt <- cat_dt[, .N, by=category]
label_dt[, labs := paste0("n=", N)]

pv_dt <- cat_dt %>%
  wilcox_test(N ~ category, alternative = "greater") %>%
  add_significance("p") %>%
  filter(p.signif != 'ns')
  # add_y_position()
pv_dt

rare <- ggplot(cat_dt, aes(x=reorder(category, N, median), y=N)) +
  geom_boxplot(aes(fill=category), width=0.5, alpha=0.4) +
  geom_dotplot(aes(fill=category), binaxis='y', stackdir='center', dotsize=0.75, alpha=0.8) +
  # geom_boxplot(width=0.5) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, alpha=0.6) +
  geom_text(data=label_dt, aes(y=-35, label=labs), size=2.5, hjust=-0.0) +
  stat_pvalue_manual(pv_dt, label = "p.signif", tip.length = 0, y.position = c(160, 175), coord.flip = TRUE) +
  # scale_fill_manual(values = c('Significantly\nbetter' = "#33A02C", 'No significant\ndifference' = "gray", 'Significantly\nworse'="#FF5F30")) +
  scale_fill_manual(values = color_vec) +
  ylab('Gene associations') +
  xlab('') +
  # guides(x = guide_axis(angle = 45)) +
  theme_cowplot() +
  coord_flip() +
  theme(legend.position = 'none',
        axis.text=element_text(size=9),
        text = element_text(size=9))
  
rare

# Number of GWAS index variants
gwas_dt <- read_parquet('/s/project/geno2pheno/data/enrichment_data/index_variants.pq') %>% as.data.table()
peak_cts <- gwas_dt[, .(peaks = .N), by=trait]
cat_dt <- merge(cat_dt, peak_cts, by='trait')

label_dt <- cat_dt[, .N, by=category]
label_dt[, labs := paste0("n=", N)]

pv_dt2 <- cat_dt %>%
  wilcox_test(peaks ~ category, alternative = "greater") %>%
  add_significance("p") %>%
  filter(p.signif != 'ns') 
pv_dt2

common <- ggplot() +
  # geom_boxplot() +
  geom_boxplot(data=cat_dt, aes(x=reorder(category, peaks, median), y=peaks, fill=category), width=0.5, alpha=0.4) +
  geom_dotplot(data=cat_dt, aes(x=reorder(category, peaks, median), y=peaks, fill=category), binaxis='y', stackdir='center', dotsize=0.35, alpha=0.8) +
  geom_text(data=label_dt, aes(x=category, y=400, label=labs), size=2.5, hjust=-0.01) +
  ylab('GWAS index variants') +
  xlab('') +
  scale_y_log10() +
  # scale_fill_manual(values = c('Significantly\nbetter' = "#33A02C", 'No significant\ndifference' = "gray", 'Significantly\nworse'="#FF5F30")) +
  scale_fill_manual(values = color_vec) +
  annotation_logticks(sides='b') +
  guides(x = guide_axis(angle = 45)) +
  stat_pvalue_manual(pv_dt2, label = "p.signif", tip.length = 0, y.position = c(4.25, 4.5), coord.flip = TRUE) +
  theme_cowplot() +
  coord_flip() +
  theme(legend.position = 'none',
        axis.text=element_text(size=9),
        text = element_text(size=9))

common

# ggarrange(herit, rare, common, nrow=1)

# ggarrange(lm_plot, ggarrange(herit, rare, common, nrow=1, labels = c("B", "C", "D")), ncol=1, labels = c("A", ""), heights = c(4,3))
fig3 <- ggarrange(lm_plot, ggarrange(herit, rare, common, ncol=1, labels = c("B", "C", "D"), font.label = list(size = 12)), ncol=2, labels = c("A", ""), widths = c(4,3), font.label = list(size = 12))
fig3

ggsave('/s/project/geno2pheno/figures/paper_figures/figure_3.png', fig3, width=18, height=13, units = "cm", bg = 'white')
