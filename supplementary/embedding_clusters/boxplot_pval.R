library(data.table)
library(magrittr)
library(rstatix)

a <- fread('/s/project/geno2pheno/funcrvp/paper_revisions/predictions/embedding_nn_dist_stats.tsv')

lolim <- 10
uplim <- 50
lab_dt <- unique(a[num_pairs>=lolim & num_pairs<=uplim, .(trait, comparison_bonf)])[, .N, by=comparison_bonf]
lab_dt

cols <- palette.colors(8, palette = "Dark2")
color_vec <- c('significantly better' = cols[[1]], 'no significant difference' = "gray")

pv_dt <- a[num_pairs>=lolim & num_pairs<=uplim] %>%
  wilcox_test(min_distance ~ comparison_bonf) %>%
  add_significance("p") %>%
  add_y_position()
pv_dt

rare <- ggplot(a[num_pairs>=lolim & num_pairs<=uplim], aes(x=comparison_bonf, y=min_distance)) +
  geom_boxplot(aes(fill=comparison_bonf), width=0.5, alpha=0.5) +
  geom_text(data=lab_dt, aes(y=-38, label=N), size=4, hjust=-0.0) +
  stat_pvalue_manual(pv_dt, label = "p.signif", tip.length = 0, y.position = c(400), coord.flip = TRUE) +
  scale_fill_manual(values = color_vec) +
  ylab('Nearest neighbor Distance for Associated Genes') +
  xlab('') +
  guides(color = guide_legend(textsize=9, keyheight = 2.5)) +
  theme_cowplot() +
  coord_flip() +
  theme(legend.position = 'none')

rare
