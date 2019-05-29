#!/usr/bin/env Rscript

## run_brmsfit_allBact-v-hormones-models.R
## last updated: 04-APR-19 by PTH

# header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))

library(parallel)
options(mc.cores = parallel::detectCores())

# load data
ASVs <- read.csv(file = here("data/ASV_table_26-JUN-2018.csv"), row.names = 1)
bTAX <- read.csv(file = here("data/bTAX_table_26-JUN-2018.csv"), row.names = 1)
#mEL  <- read.csv(file = here("data/EL_sample_data_final.csv"), row.names = 1)
mNP  <- read.csv(file = here("data/NP_sample_data_final.csv"), row.names = 1)
mNP$sp_id <- with(mNP,paste0(plot_num,'_',sub_plot_tx))

# clean up counts and re-total
fams_to_remove <- c('Enterococcaceae','Rhizobiales_Incertae_Sedis')
bASVs_to_remove <- paste0(bTAX$unique.id[bTAX$Family %in% fams_to_remove])

# add column that sums up all reads as well a host-derived reads:
ASVs <- ASVs[,!names(ASVs) %in% bASVs_to_remove]
ASVs[,'lib']  <- rowSums(ASVs)
ASVs[,'host'] <- rowSums(ASVs[,names(ASVs) %in% c("mt","cp")])
ASVs[,'bac']  <- rowSums(ASVs[,1:(length(ASVs[1,])-5)])

# sanity check these sum up to 'lib'
#ASVs[,'bac'] + ASVs[,'fungi'] + ASVs[,'host'] == ASVs[,'lib'] # yup. checks out. I can add.

#ELD <- merge(ASVs, mEL, by = "row.names")
NPD <- merge(ASVs, mNP, by = "row.names")

#ELD[,'log_ratio'] <- log(ELD[,'bac']/ELD[,'host'])
NPD[,'log_ratio'] <- log(NPD[,'bac']/NPD[,'host'])
NPD$plot_num <- factor(NPD$plot_num)
# assign order to stem_tx factor:
NPD$stem_tx <- factor(NPD$stem_tx, levels = c('MC','JATX','SATX'))
NPD$stem_tx <- relevel(NPD$stem_tx, ref = 'MC')

# construct models
# do herb_dmg==0 and ==1 separately
run_tot_mods <- function(dd, site, ...){

  #mod.names <- c('ga0','ga0b','ga1','ga2','skn0','skn0b','skn1','skn2')
  mod.names <- c('h0a','h1a','h2a','h3a')
  cat(sprintf("Running models for site %s using all bacteria as response..\n\n", site))

  # fit gaussian models
  cat(sprintf("1 of %s...", length(mod.names)))
  h0a <- brm(bf(log_ratio ~ 1), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 8000, thin = 4)

  cat(sprintf("2 of %s...", length(mod.names)))
  h1a <- brm(bf(log_ratio ~ stem_tx), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 8000, thin = 4)

  cat(sprintf("3 of %s...", length(mod.names)))
  # add plot num:
  h2a <- brm(bf(log_ratio ~ (1|plot_num)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 8000, thin = 4)

  cat(sprintf("4 of %s...", length(mod.names)))
  h3a <- brm(bf(log_ratio ~ stem_tx + (1|plot_num)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 8000, thin = 4)

  mods <- list(h0a,h1a,h2a,h3a)

  names(mods) <- mod.names

  for(m in 1:length(mods)){
    #cat(sprintf("Saving model '%s' for Family %s\n", names(mods)[m], FAMS[s]))
    saveRDS(mods[[m]], file = paste0(here("models/bact_v_hormones/allBact/"),'all-bact-mod_',site,"_",names(mods)[m],'.rds'))
  }
  return(mods)
}

# run models, saving along the way:
#EL_mods <- run_tot_mods(dd = ELD, site = 'EL')
NP_d_mods <- run_tot_mods(dd = droplevels(NPD[NPD$herb_dmg==1,]), site = 'NP')
NP_u_mods <- run_tot_mods(dd = droplevels(NPD[NPD$herb_dmg==0,]), site = 'NP')

# # for debugging:
# the_mods <- NP_d_mods

# function to run comparisons:
compare_tot_mods <- function(the_mods, DAT, SET){

  #DAT_FILE <- here("data/ELD_long_v1.csv")
  #DAT <- read.csv(DAT_FILE, header = T)
  # SET <- unlist(strsplit(DAT_FILE,"/"))
  # SET <- unlist(strsplit(SET[length(SET)],'_'))[1]
  # if (SET=='NP'){
  #   names(DAT)[1] <- 'samples'
  #   DAT$sp_id <- with(DAT, paste0(plot_num,'_',sub_plot_tx))
  # }

  #REP_TEMPLATE <- here("scripts/brmsfit_modrep_template_allBact-v-hormones.Rmd")

  # LOAD CFU MODEL
  #CFU_MOD_FILE <- here("models/cfu_m1.rds")
  #cfu_mod <- readRDS(file = CFU_MOD_FILE) # this part is hard-coded for now

  mod_names <- c('h0a','h1a','h2a','h3a')

  # first add lOO-IC to each of the models:
  for(i in seq_along(the_mods)){
    the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'loo', reloo = TRUE) # assign LOOIC to each model object in modlist
    #the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'waic') # assign wAIC to each
    assign(mod_names[i], the_mods[[i]]) # turn objects named in the_mods into objects with those names
  }

  # generate LOO objects for comparisons, below
  mods <- list(h0a,h1a,h2a,h3a)
  l1 <- loo(h0a)
  l2 <- loo(h1a)
  l3 <- loo(h2a)
  l4 <- loo(h3a)

  loo_list <- list(l1,l2,l3,l4) # put them all in one place

  loo1 <- do.call(rbind, sapply(loo_list, function(x) with(x, data.frame(model_name, looic, se_looic)), simplify = F))
  loo1$loo_id <- c('l1','l2','l3','l4') # assign IDs for referenceing LOO results

  # need to assign loo objects as call-able by name:
  for(l in seq_along(loo_list)){
    assign(loo1$loo_id[l], loo_list[[l]])
  }

  loo1 <- loo1[order(loo1$looic),] # sort this data.frame by the lowest LOOIC value
  best_mod <- get(paste0(loo1$model_name[1])) # define best model based on lowest LOOIC

  # grab top five models; compare
  loo_comp1 <- brms::compare_ic(get(paste0(loo1[1,'loo_id'])),
                                get(paste0(loo1[2,'loo_id'])),
                                get(paste0(loo1[3,'loo_id'])),
                                get(paste0(loo1[4,'loo_id']))
                                )

  # grab model coefficients from hormone models:
  coef1 <- rbind(
    tidy(h3a,prob=0.95,par_type = "non-varying"),
    tidy(h1a,prob=0.95,par_type = "non-varying")
  )
  coef1[,-1] <- round(coef1[,-1],3)
  coef1$model <- c(rep('h3a',3),rep('h1a',3))

  # write this as table:
  con <- file(here("tables/allBact-v-hormones_v1.tex"))
  open(con, 'wr')
  writeLines(kable(coef1, "latex", caption = "Model coefficients for JA and SA effects on total bacteria (site NP)", booktabs = T, escape = FALSE) %>%
               kable_styling(position = "center"), con = con, sep = ''
  )
  close(con)

  # now write table with loo-ic comparisons:
  write.csv(loo1,here("models/bact_v_hormones/allBact/allBact-v-hormones_v1_LOO.txt"),quote=F,row.names=F)
  write.csv(loo_comp1$ic_diffs__,here("models/bact_v_hormones/allBact/allBact-v-hormones_v1_LOO_comp.txt"),quote=F,row.names=T)

  # export coefficients for plotting supplemental figure:
  coef2 <- rbind(
    tidy(h3a,prob=0.5,par_type = "non-varying"),
    tidy(h1a,prob=0.5,par_type = "non-varying")
  )
  coef2[,-1] <- round(coef2[,-1],3)
  coef2$model <- c(rep('h3a',3),rep('h1a',3))

  c_for_p <- cbind(coef1[grep('stem_',coef1$term),c('term','lower')],
                   coef2[grep('stem_',coef2$term),c('lower','estimate','upper')],
                   coef1[grep('stem_',coef1$term),c('upper','model')]
  )
  c_for_p$taxon <- 'allBact'
  names(c_for_p) <- c('term','q0.025','q0.50','estimate','q0.75','q0.975','model','taxon')

  # write this to file:
  write.csv(c_for_p,file = here("models/bact_v_hormones/hormone_coefficients_for_plot.txt"),quote=F,row.names=F)

  # #### GENERATE REPORT ####
  # rmarkdown::render(input = paste0(REP_TEMPLATE),
  #                   output_format = "pdf_document",
  #                   output_file = paste0("modrep_",SET,"_",Sys.Date(),".pdf"),
  #                   output_dir = here("models/all-bacteria/"))
}

compare_tot_mods(NP_d_mods, DAT = NPD, SET = 'NP_d')
compare_tot_mods(NP_u_mods, DAT = NPD, SET = 'NP_u')

# now run script 'prep_ppcfu_rep-level.R` to generate figure output for log median across sample types as well as
# posterior difference in number of doublings between sample types
