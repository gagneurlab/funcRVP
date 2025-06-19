library(arrow)
library(scales)
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
library(diagis)
library(wesanderson)

theme_pub <- theme(legend.text = element_text(size=9),
                   axis.text=element_text(size=9),
                   text = element_text(size=9))

# Read the data
genotype = 'deepRVAT'

dt_func <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/v1NEWsplit_deepRVAT_testsplit0.25_omics_pops_genes_extended.pq') %>% as.data.table()
signif_func <- dt_func[, .(trait, gene_id, significant)][, signif_func := significant][, significant:=NULL]

dt_ols <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/at_filteredv3_deeprvat_0.25_genes_NEWsplit.pq') %>% as.data.table()
signif_ols <- dt_ols[, .(trait, gene_id, significant)][, signif_ols := significant][, significant:=NULL]

dt_test <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/at_filteredv3_deeprvat_0.25_test_genes_NEWsplit.pq') %>% as.data.table()
dt_test[, model := paste0('ols_',genotype,'_test')]
dt_test[, beta_test := beta]
signif_test <- dt_test[, .(trait, gene_id, significant)][, signif_test := significant][, significant:=NULL]

df <- rbindlist(list(dt_func[, .(model, trait, gene_id, neglog_pval, significant, beta)], dt_ols[, .(model, trait, gene_id, neglog_pval, significant, beta)]), use.names = TRUE)
df <- merge(df, dt_test[, .(trait, gene_id, beta_test)], by = c('trait','gene_id'))
df <- df[gene_id %in% dt[, unique(gene_id)]]

df[, beta_dir := beta*beta_test >= 0]
df[, rank := rank(-neglog_pval, na.last = TRUE, ties.method = 'random'), by=model]
df <- df[order(rank)]
df[, precision := cumsum(beta_dir)/rank, by=model]
df[, AD := abs(beta - beta_test)]
df[, SE := (beta - beta_test)^2]
df

# Sign precision plot
sp_df <- df[model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype)][order(rank)]
sp_df[, pd := 1-10^(-neglog_pval)]

# Define the custom transformation
linear_breaks_trans <- new_transform(
  name = "linear_breaks",
  transform = function(x) {
    sapply(x, function(val) {
      if (is.na(val)) {
        return(NA)
      } else if (val <= min(custom_breaks)) {
        return(0)
      } else if (val >= max(custom_breaks)) {
        return(1)
      } else {
        n_breaks <- length(custom_breaks)
        lower_idx <- max(which(custom_breaks <= val))
        upper_idx <- lower_idx + 1
        lower_break <- custom_breaks[lower_idx]
        upper_break <- custom_breaks[upper_idx]
        
        # Calculate the position between the lower and upper break (0 to 1)
        fraction <- (val - lower_break) / (upper_break - lower_break)
        
        # Map this fraction to the linear scale between the indices
        transformed_value <- (lower_idx - 1 + fraction) / (n_breaks - 1)
        return(transformed_value)
      }
    })
  },
  inverse = function(x) {
    if (any(is.na(x))) {
      return(rep(NA, length(x)))
    }
    if (any(x < 0, na.rm = TRUE) || any(x > 1, na.rm = TRUE)) {
      warning("Inverse of linear_breaks_trans is an approximation outside [0, 1]")
    }
    approx(x = seq(0, 1, length.out = length(custom_breaks)), y = custom_breaks, xout=x)$y
    # approx(x = seq(0, 1, length.out = length(custom_breaks)), y = custom_breaks, xout = x)$y
  },
  breaks = function(limits) {
    return(custom_breaks[custom_breaks >= limits[1] & custom_breaks <= limits[2]])
  },
  format = format_format()
)

# Your custom breaks that you want to appear linearly spaced
custom_breaks <- c(0.5, 0.75, 0.9, 0.99, 0.999, 1)

# Sign precision line plot
sp_plot <- ggplot(sp_df, aes(x=pd, y=precision)) +
  geom_line() +
  geom_hline(aes(yintercept = 0.95)) +
  geom_vline(aes(xintercept = 0.999), linetype='dashed', color='gray') +
  annotate("text", x=0.75, y=0.95, label='95% sign agreement\non test', size=3, hjust=0, parse = F) +
  xlab('Probability of direction') +
  ylab('Sign agreement\non test set') +
  theme_cowplot() +
  scale_x_continuous(
    transform = linear_breaks_trans,
    breaks = custom_breaks,
    labels = as.character(custom_breaks)
  ) +
  theme_pub + 
  theme(legend.position='none',
        panel.grid.minor.y = element_line(colour="gray85", size=0.25),
        panel.grid.major.y = element_line(colour="gray85", size=0.5))

sp_plot

model_name_vec <- c(
  "ols_deepRVAT_cov" = "Burden Test",
  "at_deepRVAT_cov" = "Burden Test",
  "ols_pLoF_cov" = "OLS pLoF",
  "omics_pops_bayesian_v108cov_deepRVAT" = "FuncRVP",
  "omics_pops_bayesian_v107cov_flat_deepRVAT" = "Constant reg.",
  "omics_pops_bayesian_v2cleansplitHO_deepRVAT" = "FuncRVP",
  "omics_pops_bayesian_v2NEWsplit_deepRVAT" = "FuncRVP",
  "omics_pops_bayesian_v1NEWsplit_deepRVAT" = "FuncRVP",
  "omics_pops_bayesian_v1NEWsplit_am_plof" = "FuncRVP"
)

model_color_vec <- c(
  "ols_deepRVAT_cov" = "#FF7F00",
  "ols_am_plof_cov" = "#FF7F00",
  "at_deepRVAT_cov" = "#FF7F00",
  "at_am_plof_cov" = "#FF7F00",
  "omics_pops_bayesian_v108cov_deepRVAT" = "#6A3D9A",
  "omics_pops_bayesian_v107cov_flat_deepRVAT" = "seagreen",
  "omics_pops_bayesian_v2cleansplitHO_deepRVAT" = "#6A3D9A",
  "omics_pops_bayesian_v2NEWsplit_deepRVAT" = "#6A3D9A",
  "omics_pops_bayesian_v1NEWsplit_deepRVAT" = "#6A3D9A",
  "omics_pops_bayesian_v1NEWsplit_am_plof" = "#6A3D9A"
)

cust_theme <- list(
  theme(text = element_text(size=15))
)

small_cut <- df[significant==TRUE & model==paste0("ols_", genotype,"_cov"), max(rank)]
ols_c <- df[rank<=small_cut]
# ols_c[, cutoff := 'Burden Test significance cutoff']
ols_c[, cutoff := paste0('Top ', small_cut, ' associations')]

big_cut <- df[significant==TRUE & model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype), max(rank)]
func_c <- df[rank<=big_cut]
# func_c[, cutoff := 'FuncRVP discovery threshold']
func_c[, cutoff := paste0('Top ', big_cut, ' associations')]

facet_dt <- rbindlist(list(ols_c, func_c), use.names = TRUE)
# facet_dt$cutoff <- factor(facet_dt$cutoff, levels = c('Burden Test significance cutoff', 'FuncRVP discovery threshold'))

facet_dt[, beta_dir := ifelse(beta * beta_test >= 0, TRUE, FALSE)]
facet_dt[, sign_precision := sum(beta_dir)/.N, by=.(model, cutoff)]
facet_dt[, RMSE := sqrt(sum((beta - beta_test)^2)/.N), by=.(model, cutoff)]
facet_dt[, MAD := median(AD), by=.(model, cutoff)]
facet_dt[['cutoff']] <- factor(facet_dt[['cutoff']], levels = c(paste0('Top ', small_cut, ' associations'),  paste0('Top ', big_cut, ' associations')))

facet_metrics <- unique(facet_dt[, .(model, cutoff, sign_precision, RMSE, MAD)])
facet_metrics[, `:=` (sp_label = paste0("Sign Prec.: ", formatC(sign_precision, digits = 3)),
                      rmse_label = paste0("RMSE: ", formatC(RMSE, digits = 3)),
                      mad_label = paste0('MAD: ', formatC(MAD, digits = 3)))]

scatter_plot <- ggplot(facet_dt, aes(beta, beta_test)) +
  geom_abline() +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_point(alpha=0.5, size=0.75) +
  geom_text(data = facet_metrics, aes(x=-3.25 , y=4 , label=mad_label), size=2.5) +
  xlab(TeX('Gene effects, Training samples')) +
  ylab(TeX('Gene effects, Test samples')) +
  facet_grid(model~cutoff, scales='free_x', labeller=labeller(model=model_name_vec)) +
  theme_cowplot() +
  theme_pub
  # theme(axis.text=element_text(size=8),
  #       text = element_text(size=8))

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

cutoff = 2000
jk_metrics <- jackknife(df[rank<=cutoff])
jk_metrics[, trait_rank_rescaled := as.integer((trait_rank/.N)*df[, max(rank)]), by=model]
jk_metrics[model==paste0("ols_", genotype,"_cov"), significant := ifelse(trait_rank<=df[significant==TRUE & model==paste0("ols_", genotype,"_cov"), max(rank)], TRUE, FALSE)]
jk_metrics[model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype), significant := ifelse(trait_rank<=df[significant==TRUE & model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype), max(rank)], TRUE, FALSE)]

label_jak <- df[significant==TRUE, .(n=max(rank), labs=paste0('n=',max(rank))), by=model]
label_jak[model==paste0("ols_", genotype,"_cov"), description:=paste0("FWER<5%\n", labs)]
label_jak[model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype), description:=paste0("pd>0.999\n", labs)]

rmse_jak <- ggplot() +
  # geom_line(data = jk_metrics[trait_rank_rescaled>25 & trait_rank_rescaled<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank_rescaled, y=avg_RMSE, color=model, linetype=significant)) +
  # geom_ribbon(data = jk_metrics[trait_rank_rescaled>25 & trait_rank_rescaled<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank_rescaled, ymin=ci_low_RMSE, ymax=ci_high_RMSE, fill=model), alpha=0.15) +
  geom_line(data = jk_metrics[trait_rank>25 & trait_rank<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank, y=avg_MAD, color=model, linetype=significant)) +
  geom_ribbon(data = jk_metrics[trait_rank>25 & trait_rank<=cutoff][sample(1:.N, .N*0.5)], aes(x=trait_rank, ymin=ci_low_MAD, ymax=ci_high_MAD, fill=model), alpha=0.15) +
  geom_vline(xintercept = df[model==paste0("omics_pops_bayesian_v1NEWsplit_", genotype) & significant==TRUE, .N], linetype='dashed', color='gray') +
  geom_vline(xintercept = df[model==paste0("ols_", genotype,"_cov") & significant==TRUE, .N], linetype='dashed', color='gray') +
  geom_text(data = label_jak, aes(x=c(250, 1115) , y=c(0.085, 0.085), label=description, hjust=-0.12), size=3, parse = F) +
  scale_color_manual(values=model_color_vec, labels=model_name_vec, name = "Model") +
  scale_fill_manual(values=model_color_vec, labels=model_name_vec, name = "Model") +
  scale_linetype_manual(values = c(`TRUE`='solid', `FALSE`='dashed'), name = "Significant") +
  xlab('Associations ranked by confidence') +
  ylab('Median Absolute Deviation') +
  theme_cowplot() +
  guides(linetype = "none") +
  theme_pub +
  theme(legend.position = c(0.04, 0.8),
        legend.box = "horizontal")
  
rmse_jak

# ggarrange(ggarrange(panelA, rmse_jak, labels = c("A", "C"), ncol=2), scatter_plot, labels = c("", "B"), ncol=1, heights = c(1,1.5))

# NEW LOEUF plot
dt_bayes <- dt_func[, .(model, trait, phenocode, gene_id, beta, pd, significant)]
dt_bayes[, abs_beta := abs(beta)]
dt_bayes <- dt_bayes[order(-abs_beta)]
dt_bayes <- dt_bayes[, rank := rank(-pd, na.last = T, ties.method = 'first'), by=model][order(rank)]
dt_bayes <- dt_bayes[significant==T]

dt_ols[, ols_significant := significant]
dt_test[, test_nom_pval := pval]
dt_test[, test_nom_significant := test_nom_pval<=0.05]

bayes_exc <- merge(dt_bayes, dt_ols[, .(trait, gene_id, pval, ols_significant)], by = c('trait', 'gene_id'), all.x = T)
bayes_exc <- merge(bayes_exc, dt_test[, .(trait, gene_id, test_nom_pval, test_nom_significant)], by = c('trait', 'gene_id'), all.x = T)
bayes_exc[, funcrvp_exclusive := (significant==T) & (ols_significant==F)]

gnomad <- fread('/s/project/geno2pheno/data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
avg_loeuf <- gnomad[transcript_type=='protein_coding', mean(oe_lof_upper, na.rm=T)]
bayes_exc <- merge(bayes_exc, gnomad[, .(gene_id, oe_lof_upper, classic_caf)], by = 'gene_id')

bayes_exc <- bayes_exc[order(rank)]
# low_cutoff <- bayes_exc[funcrvp_exclusive==T, min(rank)]
# up_cutoff <- bayes_exc[significant==T, max(rank)]

# Ranks to get Standard error of the mean
bayes_exc[, exc_count := rank(-pd, na.last=T, ties.method = 'first'), by = funcrvp_exclusive]
bayes_exc[!is.na(oe_lof_upper), mean_loeuf := cummean(oe_lof_upper), by = funcrvp_exclusive]
bayes_exc[!is.na(oe_lof_upper), se_loeuf := sqrt(cumvar(oe_lof_upper)), by = funcrvp_exclusive]
bayes_exc[!is.na(oe_lof_upper), sem_loeuf := se_loeuf/sqrt(exc_count)]
bayes_exc[!is.na(oe_lof_upper), up_ci_loeuf := mean_loeuf+1.96*sem_loeuf]
bayes_exc[!is.na(oe_lof_upper), low_ci_loeuf := mean_loeuf-1.96*sem_loeuf]

max_exc <- bayes_exc[funcrvp_exclusive==T, max(exc_count)]
max_non_exc <- bayes_exc[funcrvp_exclusive==F & significant==T, max(exc_count)]
bars <- bayes_exc[(funcrvp_exclusive==T & exc_count==max_exc) | (funcrvp_exclusive==F & exc_count==max_non_exc)]

num_lab_exc <- bars[funcrvp_exclusive==T, exc_count]
num_lab_sh <- bars[funcrvp_exclusive==F, exc_count]
bars[, plot_label:=ifelse(funcrvp_exclusive==T, paste0('FuncRVP\nexclusive\nassociations\n(',num_lab_exc,')'), paste0('Shared\nassociations\n(',num_lab_sh,')'))]

t_test <- scientific(t.test(bayes_exc[significant==T & funcrvp_exclusive==T, oe_lof_upper], bayes_exc[significant==T & funcrvp_exclusive==F, oe_lof_upper])$p.value, digits = 3)
plab <- paste0('two-sided t-test\n P=',t_test)

loeuf_bars <- ggplot(bars, aes(y=plot_label, x=mean_loeuf)) +
  geom_point(aes(fill=funcrvp_exclusive), alpha=0.85, size=2, color='black') +
  geom_errorbar(aes(xmin = low_ci_loeuf, xmax = up_ci_loeuf), width=0.25)+
  geom_vline(aes(xintercept = avg_loeuf), linetype='dashed', color='gray') +
  annotate("text", x=avg_loeuf-0.03, y=paste0('FuncRVP\nexclusive\nassociations\n(',num_lab_exc,')'), label='Genome-wide\naverage LOEUF', check_overlap=T, size=2.5, hjust=0) +
  annotate("text", x=0.81 ,y=paste0('Shared\nassociations\n(',num_lab_sh,')'), label=plab, size=2.5, hjust=0) +
  ylab("") +
  xlab("Mean LOEUF score") +
  theme_cowplot() +
  theme_pub + 
  theme(legend.position='none')
        
loeuf_bars

# Schematic -------------------------
library(png)
library(grid)
panelA <- readPNG("/s/project/geno2pheno/figures/resub_figures/beta_cartoon_new4.png")
panelA <- as.raster(panelA)
panelA <- rasterGrob(panelA, interpolate = TRUE)

last_group <- ggarrange(rmse_jak, loeuf_bars, labels=c("D", "E"), ncol=1, heights=c(1.25,1), font.label = list(size=12))
fig5 <- ggarrange(panelA, sp_plot, scatter_plot, last_group, nrow=2, ncol=2, labels=c("A", "B", "C", ""), heights = c(1,2), widths = c(1.25,1), font.label = list(size=12))
fig5

ggsave('/s/project/geno2pheno/figures/revision_figures/figure5_new.svg', fig5, width=19, height=17, units = "cm", bg = 'white')

# Upset plot of discoveries -------------------------
library(UpSetR)

assoc_dt <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/old_results/all_emb_associations.parquet') %>% as.data.table()
updt <- assoc_dt[beta_dir_cummean_all>=0.95, .(model, assoc)]

up_ls <- lapply(unique(updt$model), FUN = function(x){
  updt[model==x]$assoc
}
)
names(up_ls) <- unique(updt$model)

upplt <- upset(fromList(up_ls), order.by = "degree", point.size = 1.5, line.size = 1, text.scale = 1.25)
upplt

# Get upset plot that we can grob
uplt <- plot_grid(NULL, upplt$Main_bar, upplt$Sizes, upplt$Matrix, nrow=2, align='hv', rel_heights = c(2,1), rel_widths = c(1,2.5))
ggsave('/s/project/geno2pheno/figures/revision_figures/figure5_upset.svg', uplt, width=19, height=9, units = "cm", bg = 'white')
