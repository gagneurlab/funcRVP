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
fr <- read_parquet('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/all_models_phenopred_10k_bootstrap_more_stats.pq') %>% as.data.table()
# fr <- plt_df[embedding %in% c('Omics+PoPS', 'Enformer', 'gene2vec', 'ESM2')]
fr[, genotype:='DeepRVAT']
# fr[embedding=='Omics+PoPS', emb:='Omics+PoPS']
# fr[embedding=='Enformer', emb:='Enformer\n(ns)']
# fr[embedding=='gene2vec', emb:='gene2vec\n(ns)']
# fr[embedding=='ESM2', emb:='ESM2*']
fr[, point_color := ifelse(pval<=0.05, ifelse(N_greater>N_lesser, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
fr[, point_color_fdr := ifelse(pval_fdr<=0.05, ifelse(N_greater>N_lesser, 'Significantly\nbetter', 'Significantly\nworse'), 'No significant\ndifference')]
unique(fr$embedding)

comp_dt <- fr[, .(genotype, embedding, pval, pval_fdr, point_color, point_color_fdr)]
comp_dt[, .N, by=.(genotype, embedding)]

plt <- comp_dt[, .N, by=.(genotype, embedding, point_color_fdr)]

em_ls <- data.table(
  embedding = c('Omics+PoPS', 'Omics', 'PoPS', 'gene2vec', 'Codons'),
  point_color_fdr = c('Significantly\nworse', 'Significantly\nworse', 'Significantly\nworse', 'Significantly\nworse', 'Significantly\nbetter'),
  genotype = 'DeepRVAT',
  N = 0
)

plt <- rbindlist(list(plt, em_ls), use.names = T)

fill_vec <- c(
  'Significantly\nbetter'= '#6C8645',
  'Significantly\nworse'= '#E3B710',
  'No significant\ndifference'= 'gray'
)

plt[['point_color_fdr']] <- factor(plt[['point_color_fdr']], levels = c('Significantly\nworse', 'No significant\ndifference', 'Significantly\nbetter'))

emb_comp <- ggplot(plt, aes(x=embedding, y=N, fill=point_color_fdr)) +
  geom_col(alpha=0.85, width=0.7, color='black', position = 'dodge') +
  xlab('\nEmbedding') +
  ylab('Traits') +
  scale_fill_manual(name='Phenotype prediction\nFDR-adjusted p-value', values = fill_vec) +
  facet_wrap(~genotype) +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  theme_cowplot() +
  guides(fill = guide_legend(textsize=3, byrow = TRUE)) +  #, keyheight=1)) +
  theme(legend.position='right',
        panel.grid.major.y = element_line(colour="gray75", size=0.5),
        legend.spacing.y = unit(5, "cm"))

emb_comp
