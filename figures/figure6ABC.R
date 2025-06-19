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

theme_pub <- theme(axis.text=element_text(size=9),
                   text = element_text(size=9))

# bayes_dt <- read_parquet('/s/project/geno2pheno/novel_discoveries/repset_merged_GWAS_post_clumping_omics_pops_bayesian_v1NEWsplit_deepRVAT_plot_ready.pq') %>% as.data.table()
# bayes_dt <- fread('/s/project/geno2pheno/novel_discoveries/repset_merged_FuncRVP_associations_plot_ready.tsv')
bayes_dt <- fread('/s/project/geno2pheno/novel_discoveries/suppTable3.tsv')

bayes_dt[ols_significant==T, label:='Burden-test\n(DeepRVAT)'] 
bayes_dt[ols_significant==F & (Genebass_Backman_significant==T | DeepRVAT_significant==T), label:='Other methods'] 
bayes_dt[is.na(label), label:='Novel']

plot_dt <- bayes_dt[, .N, by=label]
plot_dt[['label']] <- factor(plot_dt[['label']], levels = c('DeepRVAT-based\nburden-test', 'Other methods', 'Novel'))

# cols <- brewer.pal(8,"Dark2")
cols <- wes_palette("Zissou1", 10, type = "continuous")

color_vec <- c(
  'DeepRVAT-based\nburden-test' = cols[[1]],
  'Other methods' = cols[[3]],
  "Dubious due to local\ncommon variant signal" = "gray", #cols[[4]],
  "Novel" = cols[[7]]
)

cats <- ggplot(plot_dt, aes(x=F, y=N, fill=label)) +
  geom_bar(position = position_stack(reverse = T), stat='identity', just=0.5, width = 0.5, color='black') +
  scale_fill_manual(name = 'FuncRVP\nassociaiton\noverlap', values = color_vec) +
  geom_text(aes(label = label, x=T), position = position_stack(reverse = T, vjust = 0.5), size=3) +
  # geom_text_repel(aes(label = label, x=T), position = position_stack(reverse = T, vjust = 0.5), size=4, point.padding = 0.25, box.padding = 1, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  # geom_label_repel(aes(label=label), position = position_stack(reverse = T, vjust = 0.5), alpha = 0.6, size=4.25, point.padding = 0.25, label.padding=1, seed = 1234, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  coord_cartesian(clip = "off") +
  geom_text(aes(label = N), position = position_stack(reverse = T, vjust = 0.5), show.legend = FALSE, size=3) + #, fontface = "bold") +
  xlab('') +
  ylab('FuncRVP associations') +
  # coord_flip() +
  theme_cowplot() +
  theme_pub +
  theme(legend.position = 'none',
        legend.justification = "center",
        # axis.line.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

cats

# LOEUF per category -------------------------------------
# gnomad <- fread('/s/project/geno2pheno/data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
# gw_avg_loeuf <- gnomad[transcript_type=='protein_coding', mean(oe_lof_upper, na.rm=T)]
# avg_loeuf <- merge(gnomad[, .(gene_id, oe_lof_upper, classic_caf)], bayes_dt, by='gene_id')
# avg_loeuf[, funcrvp_exc := ifelse(label == 'Novel', 'FuncRVP\nnovel\nassociations', 'Shared\nassociations')]
# avg_loeuf[!is.na(oe_lof_upper), mean_loeuf := mean(oe_lof_upper), by = funcrvp_exc]
# avg_loeuf[!is.na(oe_lof_upper), se_loeuf := sqrt(var(oe_lof_upper)), by = funcrvp_exc]
# avg_loeuf[!is.na(oe_lof_upper), sem_loeuf := se_loeuf/sqrt(.N), by=funcrvp_exc]
# avg_loeuf[!is.na(oe_lof_upper), up_ci_loeuf := mean_loeuf+1.96*sem_loeuf]
# avg_loeuf[!is.na(oe_lof_upper), low_ci_loeuf := mean_loeuf-1.96*sem_loeuf]
# 
# bars <- unique(avg_loeuf[, .(funcrvp_exc, mean_loeuf, up_ci_loeuf, low_ci_loeuf)])
# 
# t_test <- scientific(t.test(avg_loeuf[label=='Novel', oe_lof_upper], avg_loeuf[label!='Novel', oe_lof_upper])$p.value, digits = 3)
# plab <- paste0('two-sided t-test\n P=',t_test)
# 
# loeuf_bars <- ggplot(bars[!is.na(mean_loeuf)], aes(y=funcrvp_exc, x=mean_loeuf)) +
#   geom_point(alpha=0.85, size=2, color='black') +
#   geom_errorbar(aes(xmin = low_ci_loeuf, xmax = up_ci_loeuf), width=0.25)+
#   geom_vline(aes(xintercept = gw_avg_loeuf), linetype='dashed', color='gray') +
#   annotate("text", x=gw_avg_loeuf-0.04, y=paste0('FuncRVP\nnovel\nassociations'), label='Genome-wide\naverage LOEUF', check_overlap=T, size=3, hjust=0) +
#   annotate("text", x=0.78 ,y=paste0('Shared\nassociations'), label=plab, size=3, hjust=0) +
#   ylab("") +
#   xlab("Mean LOEUF score") +
#   theme_cowplot() +
#   theme_pub + 
#   theme(legend.position='none',
#         axis.text=element_text(size=9),
#         text = element_text(size=9))
# 
# loeuf_bars
# ggsave('/s/project/geno2pheno/figures/revision_figures/eshg_talk_loeuf_plot.png', loeuf_bars, width=12, height=6, units = "cm", bg = 'white')
# -------------------------------------

# post_clump <- bayes_dt[(ols_significant==F) & (Genebass_Backman_significant==F) & (ols_test_nom_significant==F | common_var_corrected_olsTest_significant==T)]
post_clump <- fread('/s/project/geno2pheno/novel_discoveries/suppTable4.tsv')
post_clump[, beta_ci_low := beta - 1.96*beta_se]
post_clump[, beta_ci_up := beta + 1.96*beta_se]
post_clump[, slope_ci_low := slope - 1.96*slope_SE]
post_clump[, slope_ci_up := slope + 1.96*slope_SE]
# post_clump <- bayes_dt[label=="Novel association"]
# post_clump[, label:=NULL]

post_clump[, beta_dir := slope*beta>=0]
post_clump[common_var_corrected_olsTest_significant==T & beta_dir==T, label_test:='Supports']
post_clump[common_var_corrected_olsTest_significant==F & (beta_ci_low<=slope_ci_up | beta_ci_up>=slope_ci_low), label_test:='Consistent']
post_clump[common_var_corrected_olsTest_significant==T & beta_dir==F, label_test:='Contradicts']
post_clump[label_test=='Contradicts']
# post_clump[common_var_corrected_olsTest_significant==F & (beta>=slope_ci_low & beta<=slope_ci_up)]
# post_clump[common_var_corrected_olsTest_significant==F & !(beta>=slope_ci_low & beta<=slope_ci_up)]
post_clump[['label_test']] <- factor(post_clump[['label_test']], levels = c('Contradicts', 'Consistent', 'Supports'))
test_text <- post_clump[, .N, by=label_test]

cols <- wes_palette("AsteroidCity1")
cols2 <- wes_palette("Chevalier1", 5, type = "continuous")
cols3 <- wes_palette("Zissou1", 10, type = "continuous")
col_vec <- c('Supports' = cols[[4]], 'Consistent' = cols2[[3]], 'Contradicts' = cols3[[7]])

test <- ggplot(post_clump) +
  geom_bar(aes(x=label_test, fill=label_test), stat='count', color='black', width=0.75) +
  scale_fill_manual(name = 'Associaiton\ncategory', values = col_vec) +
  geom_text(data = test_text, aes(x=label_test, y=N+4, label=N), size=3, hjust=0) +
  ylab('Novel associations') +
  xlab('Test-set gene\neffect estimate') +
  theme_cowplot() +
  theme_pub +
  coord_flip() +
  expand_limits(y = c(0,250)) +
  # guides(x = guide_axis(angle=45)) +
  theme(legend.position = 'none')

test

# post_clump[gwas_support==T & gwas_bootstrap_pval<0.05, label_gwas:='Yes']
# post_clump[is.na(label_gwas), label_gwas:='No']
# post_clump[['label_gwas']] <- factor(post_clump[['label_gwas']], levels = c('Yes', 'No'))
# gwas_text <- post_clump[, .N, by=label_gwas]



# GWAS enrichment
func_assoc <- 151
rand_assoc <- 53
total_assoc <- 241
gwas_text <- data.table(labs = c('FuncRVP\nnovel\nassociations', 'Genome-wide'),
                        N = c(func_assoc, rand_assoc),
                        y = c(func_assoc/total_assoc, rand_assoc/total_assoc),
                        y_min = c(binom.test(c(func_assoc,total_assoc-func_assoc), alternative=c("t"), conf.level=0.95)$conf.int[1], binom.test(c(rand_assoc,total_assoc-rand_assoc), alternative=c("t"), conf.level=0.95)$conf.int[1]),
                        y_max = c(binom.test(c(func_assoc,total_assoc-func_assoc), alternative=c("t"), conf.level=0.95)$conf.int[2], binom.test(c(rand_assoc,total_assoc-rand_assoc), alternative=c("t"), conf.level=0.95)$conf.int[2]))

gwas_text[['labs']] <- factor(gwas_text[['labs']], levels = c('Genome-wide', 'FuncRVP\nnovel\nassociations'))

col_vec <- c('FuncRVP\nnovel\nassociations' = cols[[4]], 'Genome-wide'='gray')

gwas <- ggplot(gwas_text, aes(x=labs, y=y)) +
  geom_bar(aes(fill=labs), stat='identity', color='black', width=0.75) +
  scale_fill_manual(values = col_vec) +
  geom_errorbar(aes(ymin=y_min, ymax=y_max), width=0.25) +
  geom_text(aes(y=y+0.1, label = N), size=3) +
  ylab('Percentage of GWAS association overlap') +
  xlab('') +
  theme_cowplot() +
  theme_pub +
  coord_flip() +
  theme(legend.position = 'none')

gwas


fig6_right <- ggarrange(test, gwas, ncol = 1, labels = c("B", "C"),heights = c(2.75,2), align='v')
fig6_right
fig6_top <- ggarrange(cats, fig6_right, nrow = 1, labels = c("A", ""), widths = c(1,2))
fig6_top

ggsave('/s/project/geno2pheno/figures/revision_figures/figure6ABC.png', fig6_top, width=18, height=10, units = "cm", bg = 'white')
# ggsave('/s/project/geno2pheno/figures/resub_figures/figure6ABC.png', fig6_top, width=18, height=10, units = "cm", bg = 'white')

# library(png)
# library(grid)
# panelD <- readPNG("/s/project/geno2pheno/figures/resub_figures/figure6D.png")
# panelD <- as.raster(panelD)
# panelD <- rasterGrob(panelD, interpolate = TRUE)
# 
# fig6 <- ggarrange(fig6_top, panelD, ncol = 1, labels = c("", "D"), heights = c(1,1))
# fig6
# 
# ggsave('/s/project/geno2pheno/figures/resub_figures/figure6_new.png', fig6, width=15, height=15, units = "cm", bg = 'white')

# fig6 <- ggarrange(cats, test, gwas, nrow = 1, labels = c("A", "B", "C"), widths = c(1,1.25,1), align = 'h')
# fig6
# ggsave('/s/project/geno2pheno/figures/resub_figures/figure6_new.svg', fig6, width=18, height=11, units = "cm", bg = 'white')
