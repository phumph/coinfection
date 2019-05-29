#!/usr/bin/env Rscript

## run_brmsfit_bASVs-v-hormones-models.R
## last updated: 04-APR-19 by PTH

# script runs Bayesian estimation of model coefficients for Family-level bASV abundance across given dataset.
# outputs are model .RDS objects fitted with params set globally in this script
# these .RDS objects get summarized and all plots outputted in subsequent script.
# this is just the model caller!

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
FAMS <- FAMS[!FAMS %in% c('Enterococcaceae','Rhizobiales_Incertae_Sedis')]

OUTDIR <- args[3]

# for debugging:
#DAT <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv", header = T)
#FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1)


#OUTDIR <- here('models/bact_v_hormones/NP/')

# ensure OUTDIR has trailing '/'
OUTDIR <- ifelse((substr(OUTDIR, nchar(OUTDIR),nchar(OUTDIR))=='/'),
                 OUTDIR,
                 paste0(OUTDIR,'/'))

## start running models:

# for debugging
# DAT <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv", header = T)
# FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1)

DAT <- DAT[DAT$bASV_count>0,]
DAT$log_ratio <- log(DAT$bASV_count / DAT$host)

# next steps:
# 1. filter DAT
# 2. Run models, saving each one

# set number of cores:
options(mc.cores = parallel::detectCores())
#mod.names <- c('ga0','ga1','ga2','skn0','skn1','skn2','skn3','skn4')

mod.names <- c('h0d','h1d','h0d2','h1d2','h0u','h1u','h0u2','h1u2')

# run models for each family in FAMS
for(s in 1:length(FAMS)){
  # create tmp dataframe for focal FAM
  dd0 <- dplyr::filter(DAT, Family == FAMS[s],
                      herb_dmg == 0)
  dd0$plot_num <- factor(dd0$plot_num)
  dd0$stem_tx <- factor(dd0$stem_tx, levels = c('MC','JATX','SATX'))
  dd0$stem_tx <- relevel(dd0$stem_tx, ref = 'MC')

  dd1 <- dplyr::filter(DAT, Family == FAMS[s],
                       herb_dmg == 1)
  dd1$plot_num <- factor(dd1$plot_num)
  dd1$stem_tx <- factor(dd1$stem_tx, levels = c('MC','JATX','SATX'))
  dd1$stem_tx <- relevel(dd1$stem_tx, ref = 'MC')

  #for debugging:
  cat(sprintf("\nWorking on Family %s...(%s of %s) \n", FAMS[s], s, length(FAMS)))

  #fit gaussian models (effectively log-linear models)
  cat(sprintf("\nWorking on model 1 of %s...", length(mod.names)))
  h0d <- brm(bf(log_ratio ~ (1|bASV)), data = dd1, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("2 of %s...", length(mod.names)))
  h1d <- brm(bf(log_ratio ~ stem_tx + (1|bASV)), data = dd1, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("3 of %s...", length(mod.names)))
  # add plot num:
  h0d2 <- brm(bf(log_ratio ~ (1|bASV) + (1|plot_num)), data = dd1, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("4 of %s...", length(mod.names)))
  h1d2 <- brm(bf(log_ratio ~ stem_tx + (1|bASV) + (1|plot_num)), data = dd1, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  # models on undamaged leaf sets
  cat(sprintf("5 of %s...", length(mod.names)))
  h0u <- brm(bf(log_ratio ~ (1|bASV)), data = dd0, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("6 of %s...", length(mod.names)))
  h1u <- brm(bf(log_ratio ~ stem_tx + (1|bASV)), data = dd0, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("7 of %s...", length(mod.names)))
  # add plot num:
  h0u2 <- brm(bf(log_ratio ~ (1|bASV) + (1|plot_num)), data = dd0, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  cat(sprintf("8 of %s...", length(mod.names)))
  h1u2 <- brm(bf(log_ratio ~ stem_tx + (1|bASV) + (1|plot_num)), data = dd0, family = gaussian, inits = "0", control = list(adapt_delta = 0.99))

  # save models
  mods <- list(h0d,h1d,h0d2,h1d2,h0u,h1u,h0u2,h1u2)
  names(mods) <- mod.names

  for(m in 1:length(mods)){
    cat(sprintf("Saving model '%s' for Family %s\n", names(mods)[m], FAMS[s]))
    saveRDS(mods[[m]], file = paste0(OUTDIR,FAMS[s],'_mod_',names(mods)[m],'.rds'))
  }
  cat("done!\n")
}
cat("\n\n(now we're really done.)")

# calling this script:
# Rscript --vanilla run_brmsfit_bASVs-v-hormones-models.R /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams_set1.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/bact_v_hormones/NP/
# Rscript --vanilla run_brmsfit_bASVs-v-hormones-models.R /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams_set2.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/bact_v_hormones/NP/
# Rscript --vanilla run_brmsfit_bASVs-v-hormones-models.R /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams_set3.csv /Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/bact_v_hormones/NP/
