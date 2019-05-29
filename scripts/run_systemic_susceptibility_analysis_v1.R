#! /usr/bin/Rscript

# script to examine shared susceptibility within plant patches
# goal is to examine correlation between damaged and un-damaged samples within the same subplot.
# expectation is that there will be positive correlation if:
  # a.) plants have variation in base-line susceptibility or exposure rates to infection by bacteria
  # b.) herbivore damage induces susceptibility systemically, so un-damaged leaves will also show increased infection

# load header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))

# load data
NPD <- read.csv(here("data/NPD_final.csv"))

NPD[,'ln_bac_host'] <- log(NPD[,'bac']/NPD[,'host'])
NPD$plot_num <- factor(NPD$plot_num)
NPD$sp_id <- factor(paste0(NPD$plot_num,"_",NPD$sub_plot_tx))

NP_sp_cor <- NPD[,names(NPD) %in% c('sp_id','herb_dmg','ln_bac_host','plot_num')]
NP_sp_cor2 <- reshape2::dcast(NP_sp_cor, formula = plot_num + sp_id ~ herb_dmg, value.var = "ln_bac_host", fun.aggregate = function(x) mean(x))

# grab only those sub_plots with a mark in both columns:
NP_sp_cor2 <- NP_sp_cor2[complete.cases(NP_sp_cor2),]
names(NP_sp_cor2) <- c("plot_num","sp_id","none","dmg")
ggplot(NP_sp_cor2, aes(x = dmg, y = none, col = plot_num, group = plot_num)) +
  geom_point(alpha = 0.4) +
  geom_line() +
  theme_phy1() #+ ggtitle("NP")

# run linear regression
sys_cor_m2 <- brm(bf(none ~ dmg), family = gaussian, data = NP_sp_cor2, cores = 4)
#sys_cor_m3 <- brm(bf(none ~ dmg + (1|plot_num)), family = gaussian, data = NP_sp_cor2, cores = 4)
saveRDS(sys_cor_m2, file = here("models/sys_cor_mod1.rds"))

me_sys1 <- marginal_effects(sys_cor_m2, "dmg", points = FALSE, plot = FALSE) # %>% ggsave(filename = here("figs/supp_fig_systemic_bacteria.pdf"), width = 2.5, height = 2)

# now re-construct plot with points:
sys_plot1 <- ggplot() +
  geom_hline(yintercept = 0, lty = 3, col = "gray80") +
  geom_vline(xintercept = 0, lty = 3, col = "gray80") +
  geom_line(data = me_sys1[[1]], aes(x = dmg, y = estimate__), col = "gray40", lwd = 1) +
  geom_line(data = me_sys1[[1]], aes(x = dmg, y = lower__), col = "gray80") +
  geom_line(data = me_sys1[[1]], aes(x = dmg, y = upper__), col = "gray80") +
  geom_point(data = NP_sp_cor2, aes(x = dmg, y = none), col = "gray20", alpha = 0.5) +
  scale_y_continuous(limits = c(-5.25,8), breaks = seq(-4,8,2)) +
  scale_x_continuous(limits = c(-5.25,8), breaks = seq(-4,8,2)) +
  xlab("ln-ratio with herbivory") +
  ylab("ln-ratio without herbivory")
ggsave(sys_plot1, file = here("figs/systemic_plot1.pdf"), width = 2.5, height = 2.25)

# spit out table of model coefficients:
sys_coefs <- tidy(sys_cor_m2)[1:3,]
sys_coefs2 <- cbind(term = sys_coefs[,1], round(sys_coefs[,-c(1)], 2)) # round and create output data.frame
sys_coefs2$coef <- with(sys_coefs2, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string
sys_coefs3 <- sys_coefs2[,c('term','coef')]
sys_coefs3$term <- c('Intercept','Slope','$\\sigma_{\\varepsilon}$')
names(sys_coefs3) <- c('Term','Estimate')

# write to .tex file:
con <- file(here("tables/Sys_cor_model_table.tex"))
open(con, 'wr')
writeLines(kable(sys_coefs3, "latex", caption = "Model coefficients systemic correlation in bacterial abundance", booktabs = T, escape = FALSE) %>%
             kable_styling(position = "center"), con = con, sep = ''
)
close(con)
