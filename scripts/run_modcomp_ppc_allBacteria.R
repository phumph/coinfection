# script to compare model performance of all bacteria models for EL and NP
# last updated 2018-SEP-20 by PTH

# header
library(knitr)
library(markdown)
library(rmarkdown)
library(brms)
library(here)

source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

OUT_DIR <- here("models/")
REP_TEMPLATE <- here("scripts/brmsfit_modrep_template_allBact.Rmd")
CFU_MOD_FILE <- here("models/cfu_m1.rds")
cfu_mod <- readRDS(file = CFU_MOD_FILE) # this part is hard-coded for now
