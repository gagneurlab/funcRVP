library(arrow)
library(magrittr)
library(data.table)
library(ggplot2)
library(ggpmisc)
library(tidyr)
library(rstatix)
library(cowplot)
library(cumstats)
library(timereg)
library(latex2exp)
library(ggpubr)


# FuncRVP
dt <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/v108cov_deepRVAT_genes_extended.pq') %>% as.data.table()
dt[, pd := 1-0.5*exp(-neglog_pval)]
dt[, significant := (pd >= 0.999)]       # CHANGING pd threshold here !!!!!!!!!!
signif_func <- dt[, .(trait, gene_id, significant)][, signif_func := significant][, significant:=NULL]


# OLS 200k
dt_ols <- read_parquet('/s/project/geno2pheno/predictions/bayesian/best_model_pred/ols_filteredv3_deepRVAT_genes_extended.pq') %>% as.data.table()
signif_ols <- dt_ols[, .(trait, gene_id, significant)][, signif_ols := significant][, significant:=NULL]


# OLS 300k
dt_test <- read_parquet('/s/project/geno2pheno/predictions/bayesian/ols_filteredv3_deepRVAT_300_genes_extended.pq') %>% as.data.table()
dt_test[, model := 'ols_deepRVAT_test']
dt_test[, beta_300k := beta]
dt_test[, pval := NULL]

signif_300k <- dt_test[, .(trait, gene_id, significant)][, signif_300k := significant][, significant:=NULL]

# Flat prior
dt_flat <- read_parquet('/s/project/geno2pheno/predictions/bayesian/v107cov_flat_deepRVAT_genes_extended.pq') %>% as.data.table()
dt_flat[, .N, by=.(model, trait)]
signif_flat <- dt_flat[, .(trait, gene_id, significant)][, signif_flat := significant][, significant:=NULL]

df <- rbindlist(list(dt[, .(model, trait, gene_id, neglog_pval, significant, beta)], dt_ols[, .(model, trait, gene_id, neglog_pval, significant, beta)]), use.names = TRUE)
df <- rbindlist(list(df, dt_flat[, .(model, trait, gene_id, neglog_pval, significant, beta)]), use.names = TRUE)
df <- merge(df, dt_test[, .(trait, gene_id, beta_300k)], by = c('trait','gene_id'))
df <- df[gene_id %in% dt[, unique(gene_id)]]

df[, beta_dir := beta*beta_300k >= 0]

# df <- df[trait %in% c('LDL_direct', 'Apolipoprotein_B', 'Mean_corpuscular_volume', 'IGF1', 'Creatinine', 'standing_height', 'Erythrocyte_count', 'Reticulocyte_count', 'Cholesterol', 'HDL_cholesterol')]

df[, rank := rank(-neglog_pval, na.last = TRUE, ties.method = 'random'), by=model]
df <- df[order(rank)]
df[, precision := cumsum(beta_dir)/rank, by=model]
df

# Rank - SSE plot
# df[, SE := ifelse(gene_id == 'ENSG00000162600' & trait=='Mean_sphered_cell_volume', 0, (beta - beta_300k)^2)]
df[, SE := (beta - beta_300k)^2]
df[, RSE := sqrt(SE)]
df[, SSE := cumsum(SE), by = model]
df[, RSSE := sqrt(SSE), by = model]
df[, MSE := cumsum(SE)/rank, by = model]
df[, RMSE := sqrt(MSE), by = model]
df[, AD := abs(beta - beta_300k)]
# df[, MAD := cummedian(AD), by = .(model)]


ols_signif <- unique(df[significant==TRUE & model=="ols_deepRVAT_cov", .(trait, gene_id)])
ols_signif[, ols_signif := TRUE]
df <- merge(df, ols_signif, on=c('trait', 'gene_id'), all.x = TRUE)
df[, ols_signif := ifelse(is.na(ols_signif), FALSE, ols_signif)]

func_signif <- unique(df[significant==TRUE & model=="omics_pops_bayesian_v108cov_deepRVAT", .(trait, gene_id)])
func_signif[, func_signif := TRUE]
df <- merge(df, func_signif, on=c('trait', 'gene_id'), all.x = TRUE)
df[, func_signif := ifelse(is.na(func_signif), FALSE, func_signif)]


# df[model=='omics_pops_bayesian_v107cov_flat_deepRVAT' & rank>=540 & rank<=550]
# df[model=='ols_deepRVAT_cov' & rank>=550 & rank<=570]
# df[trait=='Alkaline_phosphatase' & gene_id=='ENSG00000162600']
# df[rank<=1000 & SE>100]

model_name_vec <- c(
  "ols_deepRVAT_cov" = "Burden Test",
  "ols_pLoF_cov" = "OLS pLoF",
  "omics_pops_bayesian_v108cov_deepRVAT" = "FuncRVP",
  "omics_pops_bayesian_v107cov_flat_deepRVAT" = "Constant reg.",
  "omics_pops_bayesian_v98cov_deepRVAT" = "FuncRVP"
)

model_color_vec <- c(
  "ols_deepRVAT_cov" = "#FF7F00",
  "omics_pops_bayesian_v108cov_deepRVAT" = "#6A3D9A",
  "omics_pops_bayesian_v107cov_flat_deepRVAT" = "seagreen",
  "omics_pops_bayesian_v98cov_pLoF" = "#CAB2D6"
)

cust_theme <- list(
  theme(text = element_text(size=15))
)


small_cut <- df[significant==TRUE & model=="ols_deepRVAT_cov", max(rank)]
ols_c <- df[rank<=small_cut]
# ols_c[, cutoff := 'Burden Test significance cutoff']
ols_c[, cutoff := paste0('Top ', small_cut, ' associations')]

big_cut <- df[significant==TRUE & model=="omics_pops_bayesian_v108cov_deepRVAT", max(rank)]
func_c <- df[rank<=big_cut]
# func_c[, cutoff := 'FuncRVP discovery threshold']
func_c[, cutoff := paste0('Top ', big_cut, ' associations')]

facet_dt <- rbindlist(list(ols_c, func_c), use.names = TRUE)
# facet_dt$cutoff <- factor(facet_dt$cutoff, levels = c('Burden Test significance cutoff', 'FuncRVP discovery threshold'))

facet_dt[, beta_dir := ifelse(beta * beta_300k >= 0, TRUE, FALSE)]
facet_dt[, sign_precision := sum(beta_dir)/.N, by=.(model, cutoff)]
facet_dt[, RMSE := sqrt(sum((beta - beta_300k)^2)/.N), by=.(model, cutoff)]
facet_dt[, MAD := median(AD), by=.(model, cutoff)]

facet_metrics <- unique(facet_dt[, .(model, cutoff, sign_precision, RMSE, MAD)])
facet_metrics[, `:=` (sp_label = paste0("Sign Prec.: ", formatC(sign_precision, digits = 3)),
                      rmse_label = paste0("RMSE: ", formatC(RMSE, digits = 3)),
                      mad_label = paste0('MAD: ', formatC(MAD, digits = 3)))]


scatter_plot <- ggplot(facet_dt[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT'], aes(beta, beta_300k)) +
  geom_abline() +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_point(alpha=0.5, size=0.75) +
  geom_text(data = facet_metrics[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT'], aes(x=-4 , y=4 , label=mad_label), size=2.5) +
  # geom_text(data = facet_metrics[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT'], aes(x=-4.5 , y=3.2 , label=sp_label), size=4.5) +
  # geom_text(data = facet_metrics[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT'], aes(x=-4.5 , y=4.5 , label=rmse_label), size=4.5) +
  # xlab(TeX('$\\widehat{\\beta}$, discovery cohort')) +
  # ylab(TeX('$\\widehat{\\beta}$, evaluation cohort')) +
  xlab(TeX('Gene effects, discovery cohort')) +
  ylab(TeX('Gene effects, evaluation cohort')) +
  facet_grid(model~cutoff, scales='free_x', labeller=labeller(model=model_name_vec)) +
  theme_cowplot() +
  theme(axis.text=element_text(size=8),
        text = element_text(size=8))

scatter_plot



# Function returns result of jackknife confidence intervals of gene recall
jackknife <- function(results_dt){
  jk_list <- lapply(results_dt[, unique(trait)], FUN = function(Trait){
    results_dt <- results_dt[trait!=Trait][order(-neglog_pval)]
    
    results_dt[trait!=Trait, trait_rank:=1:.N, by = .(model)]
    results_dt[trait!=Trait, trait_rank_norm:=trait_rank/.N, by = .(model)]
    # results_dt[trait!=Trait, trait_RSSE:=sqrt(cumsum(SE)), by = .(model)]
    
    results_dt[trait!=Trait, trait_RMSE:=sqrt(cumsum(SE)/trait_rank), by = .(model)]
    results_dt[trait!=Trait, trait_SAD:=cumsum(AD), by = .(model)]
    results_dt[trait!=Trait, trait_MAD:=cummedian(AD), by = .(model)]
    
    results_dt[,.(model, gene_id, trait_rank, trait_rank_norm, trait_RMSE, trait_MAD, trait_SAD)]
  })
  
  names(jk_list) <- results_dt[, unique(trait)]
  jk_results_dt <- rbindlist(jk_list, idcol = "trait")
  
  jk_metrics <- jk_results_dt[, list(avg_MAD = mean(trait_MAD),
                                     med_MAD = median(trait_MAD),
                                     ci_low_MAD = quantile(trait_MAD, 0.025),
                                     ci_high_MAD = quantile(trait_MAD, 0.975),
                                     avg_RMSE = mean(trait_RMSE),
                                     med_RMSE = median(trait_RMSE),
                                     ci_low_RMSE = quantile(trait_RMSE, 0.025),
                                     ci_high_RMSE = quantile(trait_RMSE, 0.975)), by=.(model, trait_rank)]
  
  return(jk_metrics) 
}

cutoff = 1000
jk_metrics <- jackknife(df[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT' & rank<=cutoff])
jk_metrics[, trait_rank_rescaled := as.integer((trait_rank/.N)*df[, max(rank)]), by=model]
jk_metrics[model=='ols_deepRVAT_cov', significant := ifelse(trait_rank<=df[significant==TRUE & model=='ols_deepRVAT_cov', max(rank)], TRUE, FALSE)]
jk_metrics[model=='omics_pops_bayesian_v108cov_deepRVAT', significant := ifelse(trait_rank<=df[significant==TRUE & model=='omics_pops_bayesian_v108cov_deepRVAT', max(rank)], TRUE, FALSE)]
# jk_metrics[model=='omics_pops_bayesian_v107cov_flat_deepRVAT', significant := ifelse(trait_rank<=df[significant==TRUE & model=='omics_pops_bayesian_v107cov_flat_deepRVAT', max(rank)], TRUE, FALSE)]

label_jak <- df[significant==TRUE, .(n=max(rank), labs=paste0('n=',max(rank))), by=model]

rmse_jak <- ggplot() +
  # geom_line(data = jk_metrics[trait_rank_rescaled>25 & trait_rank_rescaled<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank_rescaled, y=avg_RMSE, color=model, linetype=significant)) +
  # geom_ribbon(data = jk_metrics[trait_rank_rescaled>25 & trait_rank_rescaled<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank_rescaled, ymin=ci_low_RMSE, ymax=ci_high_RMSE, fill=model), alpha=0.15) +
  geom_line(data = jk_metrics[trait_rank>25 & trait_rank<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank, y=avg_MAD, color=model, linetype=significant)) +
  geom_ribbon(data = jk_metrics[trait_rank>25 & trait_rank<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank, ymin=ci_low_MAD, ymax=ci_high_MAD, fill=model), alpha=0.15) +
  geom_text(data = label_jak[model != 'omics_pops_bayesian_v107cov_flat_deepRVAT'], aes(x=n , y=0.1 , label=labs, hjust=-0.12), size=2.5) +
  scale_color_manual(values=model_color_vec, labels=model_name_vec, name = "Model") +
  scale_fill_manual(values=model_color_vec, labels=model_name_vec, name = "Model") +
  scale_linetype_manual(values = c(`TRUE`='solid', `FALSE`='dashed'), name = "Significant") +
  # scale_alpha_manual(values = c(`TRUE`= 0.4, `FALSE`=0.1)) +
  # scale_x_continuous(labels = function(x) as.integer(x * df[, max(rank)])) +
  xlab('Genes ranked by confidence') +
  ylab('Median absolute\ndeviation') +
  theme_cowplot() +
  geom_vline(xintercept = df[model=='omics_pops_bayesian_v108cov_deepRVAT' & significant==TRUE, .N], linetype='dashed', color='gray') +
  geom_vline(xintercept = df[model=='ols_deepRVAT_cov' & significant==TRUE, .N], linetype='dashed', color='gray') +
  # scale_y_log10(limits=c(0.1,7.5)) +
  # annotation_logticks(sides = 'l') +
  # theme(legend.position = "bottom") +
  # guides(x = guide_axis(angle=60)) +
  theme(legend.position = c(0.02, 0.8),
        legend.box = "horizontal",
        legend.text = element_text(size=8),
        axis.text=element_text(size=8),
        text = element_text(size=8))
  
rmse_jak


library(png)
library(grid)

# panelA <- readPNG("/s/project/geno2pheno/figures/schema_beta.png")
# panelA <- readPNG("/s/project/geno2pheno/figures/beta_cartoon.png")
panelA <- readPNG("/s/project/geno2pheno/figures/beta_cartoon_new2.png")
panelA <- as.raster(panelA)
panelA <- rasterGrob(panelA, interpolate = TRUE)

# ggarrange(ggarrange(panelA, rmse_jak, labels = c("A", "C"), ncol=2), scatter_plot, labels = c("", "B"), ncol=1, heights = c(1,1.5))

# Sign precision plot
sp_df <- df[model=='omics_pops_bayesian_v108cov_deepRVAT'][order(rank)]
sp_df[, pd := 1-0.5*exp(-neglog_pval)]
sp_df[, pd_bin := qcut(pd, breaks=c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1))]

label_df <- sp_df[, .(num_assoc = .N), by=pd_bin]
label_df[, labs := paste0("associations\n", num_assoc)]
sp_df[, sp := sum(beta_dir)/.N, by=.(trait, pd_bin)]
sp_df <- unique(sp_df[, .(pd_bin, trait, sp)])

sp_plot <- ggplot(sp_df, aes(x=pd_bin, y=sp)) +
  geom_boxplot(fill="#6A3D9A", alpha=0.5) +
  # geom_jitter(fill="#6A3D9A", alpha=0.6) +
  # geom_dotplot(fill="#6A3D9A", binaxis='y', stackdir='center', dotsize=0.6, alpha=0.8) +
  # geom_text(data=label_df, aes(x=pd_bin, y=0.1, label=labs), size=4) +
  guides(x = guide_axis(angle=45)) +
  xlab('Estimated probability of direction') +
  ylab('Sign agreement\non evaluation cohort') +
  theme_cowplot() +
  theme(axis.text=element_text(size=8),
        text=element_text(size=8))

sp_plot


# LOEUF plot
gb_plof <- read_parquet('/s/project/geno2pheno/data/replication_sets/genebass_plof.pq') %>% as.data.table()
colnames(gb_plof) <- c("annotation", "gene_id", "phenocode", "GB_beta", "GB_beta_SE", "GB_pval", "trait")

signif <- df[, .(trait, gene_id, model, significant)]
signif <- dcast(signif, ...~model, value.var = "significant")

loeuf_dt <- merge(signif, gb_plof[, .(trait, gene_id, GB_beta)], by = c('trait', 'gene_id'))

# caf <- fread("/s/project/geno2pheno/data/plof_caf.tsv")
# dt <- merge(dt, caf, by = "gene_id")

gene_names <- fread('/s/project/geno2pheno/data/hgnc2ensg.tsv')
gene_names <- unique(gene_names[, .(gene_id = `Ensembl gene ID`, symbol = `Approved symbol`)])

loeuf_dt <- merge(loeuf_dt, gene_names, by = 'gene_id')
loeuf_dt[, both_signif := ols_deepRVAT_cov & omics_pops_bayesian_v108cov_deepRVAT]
loeuf_dt[, only_FuncRVAT_signif := !ols_deepRVAT_cov & omics_pops_bayesian_v108cov_deepRVAT]
loeuf_dt[, only_OLS_signif := ols_deepRVAT_cov & !omics_pops_bayesian_v108cov_deepRVAT]
loeuf_dt[, any_signif := ols_deepRVAT_cov | omics_pops_bayesian_v108cov_deepRVAT]

loeuf_dt[, discovery_method := ifelse(both_signif == TRUE, "Shared discoveries", 
                                      ifelse(only_FuncRVAT_signif == TRUE, "FuncRVP exclusive\ndiscoveries", 
                                             ifelse(only_OLS_signif == TRUE, "BD exclusive\ndiscoveries", "Not discovered")
                                             ))]

loeuf_dt[, label := paste(trait, symbol)]

gnomad <- fread('/s/project/geno2pheno/data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
all_loeuf <- gnomad[, .(gene_id, oe_lof_upper)]

loeuf <- merge(loeuf_dt[, .(gene_id, trait, discovery_method)], all_loeuf, by = 'gene_id')
loeuf[, .N, by=discovery_method]
loeuf <- loeuf[discovery_method %in% c('FuncRVP exclusive\ndiscoveries', 'BD exclusive\ndiscoveries', 'Shared discoveries')]

omim_all <- fread('/s/project/geno2pheno/data/ukb_omim_clean.tsv')
omim_all <- omim_all[trait %in% unique(dt$trait)]
omim <- copy(omim_all)
# omim <- unique(omim[, `:=` (`MIM Number` = NULL, trait = 'OMIM genes', discovery_method= 'OMIM genes')])
omim <- omim[, `:=` (`MIM Number` = NULL, trait = 'OMIM genes', discovery_method= 'OMIM genes')]
om_dt <- merge(all_loeuf, omim, by = 'gene_id')
all_loeuf[, `:=` (trait = 'All genes', discovery_method = 'All genes')]

loeuf <- rbindlist(list(all_loeuf, om_dt, loeuf), use.names = TRUE)
loeuf$discovery_method <- factor(loeuf$discovery_method, levels = c('FuncRVP exclusive\ndiscoveries', 
                                                                    'Shared discoveries', 
                                                                    'BD exclusive\ndiscoveries', 
                                                                    'All genes', 
                                                                    'OMIM genes'))


loeuf[, median(oe_lof_upper, na.rm = TRUE), by=discovery_method]

pv_l <- loeuf %>%
  wilcox_test(oe_lof_upper ~ discovery_method, alternative = 'less') %>%
  add_significance("p") %>%
  filter((group1 == 'FuncRVP exclusive\ndiscoveries') | (group2 == 'FuncRVP exclusive\ndiscoveries')) %>%
  filter(p.signif != 'ns') %>%
  add_y_position()
pv_l

label_l <- loeuf[, .N, by=discovery_method]
label_l[, labs := paste0("n=", N)]

loeuf[, med := median(oe_lof_upper, na.rm = TRUE), by = discovery_method]

loeuf_color_vec <- c(
  "BD exclusive\ndiscoveries" = "#FF7F00",
  "FuncRVP exclusive\ndiscoveries" = "#6A3D9A",
  "Shared discoveries" = "darkgreen",
  "All genes" = "ivory2", # "lightblue",
  "OMIM genes" = "ivory4" # 'maroon'
)

loeuf_plot <- ggplot() +
  geom_boxplot(data=loeuf, aes(x=discovery_method, y=oe_lof_upper, fill=discovery_method), alpha=0.6) +
  geom_text(data=label_l, aes(x=discovery_method, y=-0.15, label=labs), size=2.5) +
  geom_text(data=unique(loeuf[, .(discovery_method, med)]), aes(x=discovery_method, y=med, label = med), size = 2.5, vjust = -1) +
  scale_fill_manual(name="", values=loeuf_color_vec) +
  theme_cowplot() +
  # background_grid(major = 'x', size.major = 0.5) +
  scale_y_continuous(minor_breaks = seq(0, 2, 0.1)) +
  # coord_flip() +
  stat_pvalue_manual(pv_l, y.position = c(2.05, 2.15, 2.25), label = "p.signif", tip.length = 0, coord.flip = FALSE, size = 3) +
  guides(x = guide_axis(angle=45)) +
  xlab('') +
  ylab('LOEUF score') +
  theme(legend.position='none',
        panel.grid.minor.y = element_line(colour="gray85", size=0.25),
        panel.grid.major.y = element_line(colour="gray75", size=0.5),
        axis.text=element_text(size=8),
        text = element_text(size=8)) 

loeuf_plot

# ggarrange(ggarrange(sp_plot, panelA, labels = c("A", "B"), ncol=1, heights=c(1.5, 1)), ggarrange(scatter_plot, rmse_jak, labels = c("C", "D"), ncol=1, heights=c(2, 1)), ncol=2, widths=c(1,1.5))

# ggarrange(sp_plot, scatter_plot, panelA, rmse_jak, labels = c("A", "C", "B", "D"), ncol=2, nrow=2, heights=c(1.5,1), widths=c(1, 1.5))

# ggarrange(ggarrange(panelA, sp_plot, labels = c("A", "B"), ncol=2), ggarrange(scatter_plot, labels = c("C"), ncol=1), ggarrange(rmse_jak, loeuf_plot, labels = c("D", "E"), ncol=2), nrow=3, heights=c(1,1.5,1))
# ggarrange(ggarrange(panelA, sp_plot, labels = c("A", "B"), ncol=2), ggarrange(scatter_plot, labels = c("C"), ncol=1), ggarrange(rmse_jak, loeuf_plot, labels = c("D", "E"), ncol=2), nrow=3, heights=c(1,1.5,1))
# 
# ggarrange(ggarrange(panelA, sp_plot, labels = c("A", "B"), ncol=1, heights=c(1,2)), ggarrange(scatter_plot, labels = c("C"), ncol=1), ggarrange(rmse_jak, loeuf_plot, labels = c("D", "E"), ncol=1, heights = c(1, 2)), nrow=1, widths=c(1,1,1))
# 
# ggarrange(ggarrange(panelA, rmse_jak, sp_plot, loeuf_plot, labels = c("A", "D", "B", "E"), ncol=2, nrow=2, heights=c(1,2), widths=c(1, 1)), ggarrange(scatter_plot, labels = c("C"), ncol=1), ncol=1)

last_group <- ggarrange(rmse_jak, loeuf_plot, labels=c("D", "E"), ncol=1, heights=c(1,1.5), font.label = list(size=12))
# last_group
fig4 <- ggarrange(panelA, NULL, sp_plot, scatter_plot, NULL, last_group, nrow=2, ncol=3, labels=c("A","", "B", "C", "", ""), heights = c(1,2), widths = c(1.25,0.05,1), font.label = list(size=12))
fig4

ggsave('/s/project/geno2pheno/figures/paper_figures/figure_4.png', fig4, width=18, height=20, units = "cm", bg = 'white')
