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

# Embedding comparison
plt_df <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_models_phenopred_10k_bootstrap_more_stats.pq') %>% as.data.table()
fr <- plt_df[embedding %in% c('Omics+PoPS', 'Enformer', 'gene2vec', 'ESM2')]
fr[, genotype:='DeepRVAT']
fr[embedding=='Omics+PoPS', emb:='Omics+PoPS']
fr[embedding=='Enformer', emb:='Enformer\n(ns)']
fr[embedding=='gene2vec', emb:='gene2vec\n(ns)']
fr[embedding=='ESM2', emb:='ESM2*']
fr[, point_color := ifelse(pval<=0.05, ifelse(N_greater>N_lesser, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
fr[, point_color_fdr := ifelse(pval_fdr<=0.05, ifelse(N_greater>N_lesser, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
fr[['emb']] <- factor(fr[['emb']], levels = c('Omics+PoPS', 'Enformer\n(ns)', 'gene2vec\n(ns)', 'ESM2*', 'Omics+PoPS\n(+)', 'Omics+PoPS*'))
fr

plof <- fread("/s/project/geno2pheno/results/bootstrap_r2_results_with_SE_v1NEWsplit_omics_pops_pLoF.csv")
plof[, `:=` (emb='Omics+PoPS*', genotype='pLoF')]
plof[, pval_fdr := p.adjust(plof$pval, method = 'BY')]
plof[, point_color := ifelse(pval<=0.05, ifelse(model_r2>=baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
plof[, point_color_fdr := ifelse(pval_fdr<=0.05, ifelse(model_r2>=baseline_model_r2, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
plof[, .N, by=point_color]
plof[, .N, by=point_color_fdr]

comp_dt <- rbindlist(list(fr[, .(genotype, emb, pval, pval_fdr, point_color, point_color_fdr)], plof[, .(genotype, emb, pval, pval_fdr, point_color, point_color_fdr)]), use.names = T)
comp_dt[, .N, by=.(genotype, emb)]

plt <- comp_dt[, .N, by=.(genotype, emb, point_color_fdr)]

em_ls <- data.table(
  emb = c('Omics+PoPS', 'gene2vec\n(ns)', 'Omics+PoPS*', 'Omics+PoPS*'),
  point_color_fdr = c('Significantly\nworse', 'Significantly\nworse', 'Significantly\nworse', 'Significantly\nbetter'),
  genotype = c('DeepRVAT', 'DeepRVAT', 'pLoF', 'pLoF'),
  N = 0
)

plt <- rbindlist(list(plt, em_ls), use.names = T)
# plt <- rbindlist(list(plt, list('DeepRVAT', 'Omics+PoPS', 'Significantly\nworse', 0)))

fill_vec <- c(
  'Significantly\nbetter'= '#6C8645',
  'Significantly\nworse'= '#E3B710',
  'No significant\ndifference'= 'gray'
)

plt[['point_color_fdr']] <- factor(plt[['point_color_fdr']], levels = c('Significantly\nworse', 'No significant\ndifference', 'Significantly\nbetter'))

emb_comp <- ggplot(plt, aes(x=emb, y=N, fill=point_color_fdr)) +
  geom_col(alpha=0.85, width=0.7, color='black', position = 'dodge') +
  xlab('\nEmbedding') +
  ylab('Traits') +
  scale_fill_manual(name='Phenotype prediction\nFDR-adjusted p-value', values = fill_vec) +
  facet_grid(~genotype, scales = 'free_x', space='free') +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  theme_cowplot() +
  guides(fill = guide_legend(textsize=3, byrow = TRUE)) +  #, keyheight=1)) +
  theme(legend.position='right',
        panel.grid.major.y = element_line(colour="gray75", size=0.5),
        legend.spacing.y = unit(5, "cm"))

emb_comp


## Upset plot -----------------------
library(UpSetR)

listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5, 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
updt <- fr[point_color_fdr=='Significantly\nbetter', .(embedding, trait)]

up_ls <- lapply(unique(updt$embedding), FUN = function(x){
  updt[embedding==x]$trait
  }
)
names(up_ls) <- unique(updt$embedding)

upplt <- upset(fromList(up_ls), order.by = "freq", point.size = 3.5, line.size = 1.5, text.scale = 1.25)
upplt
 
# Get upset plot that we can grob
uplt <- plot_grid(NULL, upplt$Main_bar, upplt$Sizes, upplt$Matrix, nrow=2, align='hv', rel_heights = c(2,1), rel_widths = c(1,2.5))


fig4 <- plot_grid(emb_comp, uplt, nrow=2, align='hv', rel_heights = c(1,1), labels = c("A", "B"))
fig4

ggsave('/s/project/geno2pheno/figures/revision_figures/figure4_new.svg', fig4, width=18, height=20, units = "cm", bg = 'white')
