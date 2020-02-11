#!/usr/bin/env Rscript

## run_brmsfit_bASVs.R
## last updated: 27-JUN-18 by PTH

# script runs Bayesian estimation of model coefficients for Family-level bASV abundance across given dataset.
# outputs are model .RDS objects fitted with params set in this script
# these .RDS objects get summarized and all plots outputted in subsequent script.

# inputs:
  # DAT = input data.frame containing all relevant information (full path)
  # FAMS = input .csv with target families listed, one per row (full path)
  # OUTDIR = output directory in which to store .rds model files

# outputs:
  # in OUTDIR directory, outputs .RDS for each model (n=6 per family) for each family listed in FAMS

suppressPackageStartupMessages(library(brms, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(parallel, quietly = TRUE))

# need to grab command line arguments (there are two)
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please supply three input files: <1 = DATA file> <2 = metadata file> <3 = output directory (defaults to ./)\n", call.=FALSE)
}

DAT <- read.csv(args[1], header = T)
FAMS <- read.csv(args[2], header = F)
FAMS <- paste0(FAMS$V1)
OUTDIR <- args[3]

# ensure OUTDIR has trailing '/'
OUTDIR <- ifelse((substr(OUTDIR, nchar(OUTDIR),nchar(OUTDIR))=='/'),
                 OUTDIR,
                 paste0(OUTDIR,'/'))

## start running models:

# for debugging
# DAT <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv", header = T)
# FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1)

DAT <- DAT[DAT$bASV_count>0,]
DAT$log_ratio <- log(DAT$bASV_count/DAT$host)

# next steps:
  # 1. filter DAT
  # 2. Run models, saving each one

# set number of cores:
options(mc.cores = parallel::detectCores())
mod.names <- c('ga0','ga1','ga2','skn0','skn1','skn2','skn3','skn4')

# run models for each family in FAMS
for(s in 1:length(FAMS)){

  # create tmp dataframe for focal FAM
  dd <- dplyr::filter(DAT, Family == FAMS[s])
  dd$herb_dmg <- factor(dd$herb_dmg)

  #for debugging:
  cat(sprintf("\nWorking on Family %s...(%s of %s) \n", FAMS[s], s, length(FAMS)))

  #fit gaussian models (effectively log-linear models)
  cat(sprintf("\nWorking on model 1 of %s...", length(mod.names)))
  ga0 <- brm(bf(log_ratio ~ (1|bASV)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("2 of %s...", length(mod.names)))
  ga1 <- brm(bf(log_ratio ~ herb_dmg + (1|bASV)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("3 of %s...", length(mod.names)))
  ga2 <- brm(bf(log_ratio ~ herb_dmg + (1|bASV), sigma = ~ herb_dmg), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  #ga3 <- brm(bf(log_ratio ~ herb_dmg + (1 + herb_dmg|bASV)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))
  #ga4 <- brm(bf(log_ratio ~ herb_dmg + (1 + herb_dmg|bASV), sigma = ~ herb_dmg), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  # fit skew-normal models
  skn0 <- brm(bf(log_ratio ~ (1|bASV)), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("4 of %s...", length(mod.names)))
  skn1 <- brm(bf(log_ratio ~ herb_dmg + (1|bASV)), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("5 of %s...", length(mod.names)))
  skn2 <- brm(bf(log_ratio ~ herb_dmg + (1|bASV), sigma ~ herb_dmg), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("6 of %s...", length(mod.names)))
  skn3 <- brm(bf(log_ratio ~ herb_dmg + (1 + herb_dmg|bASV)), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("7 of %s...", length(mod.names)))
  skn4 <- brm(bf(log_ratio ~ herb_dmg + (1 + herb_dmg|bASV), sigma ~ herb_dmg), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99))
  cat(sprintf("8 of %s...", length(mod.names)))

  # save models
  mods <- list(ga0,ga1,ga2,skn0,skn1,skn2,skn3,skn4)
  names(mods) <- mod.names

  for(m in 1:length(mods)){
    cat(sprintf("Saving model '%s' for Family %s\n", names(mods)[m], FAMS[s]))
    saveRDS(mods[[m]], file = paste0(OUTDIR,FAMS[s],'_mod_',names(mods)[m],'.rds'))
  }
  cat("done!\n")
}

cat("\n\n(now we're really done.)")

# execute script via the terminal
# Rscript run_brmsfit_bASVs.R /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/ELD_long_v1.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/
