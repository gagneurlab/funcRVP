library(ggplot2)
library(data.table)
library(extrafont)
library(xkcd)
# vignette("xkcd-intro")


# x <- seq(0,1,0.01)
x <- rnorm(250, sd=0.5)
eps <- rnorm(length(x), sd = 0.05)

# m=2

# y_func <- 2.1*x 
# y_ols <- 1.9*x

disco_dt <- data.table(
  'gene impairment' = x,
  'phenotype' = 0.15*x + rnorm(length(x), sd = 0.1),
  'beta FuncRVAT' = 0.17*x,
  'beta OLS' = 0.13*x 
)

ggplot() +
  geom_point(data=disco_dt, aes(x=`gene impairment`, y=phenotype), alpha=0.5) +
  geom_line(data=disco_dt, aes(x=`gene impairment`, y=`beta OLS`, color="beta OLS"), size = 1.5) +
  geom_line(data=disco_dt, aes(x=`gene impairment`, y=`beta FuncRVAT`,color="beta FuncRVAT"), size = 1.5) +
  # ggtitle('Discovery cohort') +
  theme_cowplot() +
  scale_color_manual(values = c("beta OLS" = "#FF7F00", "beta FuncRVAT" = "#6A3D9A"),
                     labels = c("beta OLS" = expression(bold(beta['OLS'])), "beta FuncRVAT" = expression(bold(beta['FuncRVAT'])))) +
  labs(color = "") +
  theme(legend.position = c(0.6, 0.2),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=20, face="bold"))

  # annotate('text', x=0.68*plot_dt[, max(`gene impairment`)], y=0.75*plot_dt[, max(`beta FuncRVAT`)], label = "beta FuncRVAT", size = unit(5, "pt")) +
  # annotate('text', x=0.65*plot_dt[, max(`gene impairment`)], y=0.68*plot_dt[, max(`beta OLS`)], label = "beta OLS", size = unit(5, "pt")) +

  

val_dt <- data.table(
  'gene impairment' = x,
  'phenotype' = 0.15*x + rnorm(length(x), sd = 0.1),
  # 'beta FuncRVAT' = 2.2*x,
  'beta OLS' = 0.14*x 
)

ggplot(val_dt, aes(x=`gene impairment`)) +
  geom_point(aes(y=phenotype), alpha=0.5) +
  # geom_line(aes(y=`beta OLS`), color="#00AFFF", size = 1.25) +
  geom_line(data=disco_dt, aes(x=`gene impairment`, y=`beta OLS`, color="beta OLS"), size = 1.25) +
  ggtitle('Replication cohort') +
  theme_cowplot() +
  scale_color_manual(values = c("beta OLS" = "#00AFFF"),
                     labels = c("beta OLS" = expression(bold(beta['OLS'])))) +
  labs(color = "") +
  theme(legend.position = c(0.6, 0.2),
        axis.text = element_text(size=15),
        axis.title = element_text(size=18, face="bold"),
        legend.text = element_text(size=18, face="bold"))
  # annotate('text', x=0.5*plot_dt[, max(`gene impairment`)], y=0.85*plot_dt[, max(`beta OLS`)], label = "beta OLS", size = unit(5, "pt")) +
  


b = rnorm(250, sd=0.5)
comp_dt <- data.table(
  'beta discovery' = b,
  'beta validation' = b + rnorm(length(b), sd = 0.25) 
)

ggplot(comp_dt, aes(x=`beta validation`, y=`beta discovery`)) +
  geom_point(alpha=0.5) +
  geom_abline() +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_vline(xintercept=0, linetype='dashed') +
  xlab(TeX('$\\beta$, discovery cohort')) +
  ylab(TeX('$\\beta$, replication cohort')) +
  theme_cowplot() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18, face="bold"))

