# script to run model for systemic susceptibility
# last updated 2018-SEP-20 by PTH

#header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))


# this model will ipmly that we should be including subplot as random intercept effect in the bacterial models. In the bASV-level models, this term gets neglected; instead, any subplot-level intercept variation gets stuck into the residual.
# These are paired samples after all.

NPD$sp_id <- paste0(NPD$plot_num,"_",NPD$sub_plot_tx)
NP_sp_cor <- NPD[,names(NPD) %in% c('sp_id','herb_dmg','ln_bac_host')]
NP_sp_cor2 <- dcast(NP_sp_cor, formula = sp_id ~ herb_dmg, value.var = "ln_bac_host", fun.aggregate = function(x) mean(x))

# grab only those sub_plots with a mark in both columns:
NP_sp_cor2 <- NP_sp_cor2[complete.cases(NP_sp_cor2),]
names(NP_sp_cor2) <- c("sp_id","none","dmg")
ggplot(NP_sp_cor2, aes(x = dmg, y = none)) + geom_point(alpha = 0.4) + theme_phy1() #+ ggtitle("NP")

# run model:
sys_cor_m2 <- brm(bf(none ~ dmg), family = gaussian, data = NP_sp_cor2)

# save model:

# generate table-form output for supplemental:

