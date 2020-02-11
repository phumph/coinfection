#!/usr/bin/env Rscript

# run_CFU_mod_EL.R
# script to run model examining log10 CFU in calibration set from EL (n=101 samples)

source("./phy_functions.R")
source("./phy_header.R")

# read data
ELD <- read.csv(file = file.path("../data/ELD_final.csv"))

# grab subset of rows where logLeafCFU>0:
ELc <- ELD[complete.cases(ELD[,'logLeafCFU']),]

cfu1 <- ggplot(ELc, aes(y = logLeafCFU, x = factor(herb_dmg), col = factor(herb_dmg),alpha = 0.2)) +
  geom_jitter(width = 0.15) + xlab("") +
  # stat_summary(geom = "crossbar",
  # fun.ymin = function(z) { quantile(z,0.25) },
  # fun.ymax = function(z) { quantile(z,0.75) },
  # fun.y = median,
  # size = 0.2) +
  geom_boxplot(width = 0.8) +
  ylab("log10 CFU g leaf") +
  scale_y_continuous(limits = c(0.5,10), breaks = seq(2,10,2)) +
  scale_color_manual(values = c("gray40","steelblue")) +
  theme_phy1() + theme(legend.position = "none")

# run quick t-test to report in SI:
t.test(ELc$logLeafCFU ~ ELc$herb_dmg)
