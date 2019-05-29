# Total herbivory analysis for EL and NP
# last updated 2018-SEP-20 by PTH

# header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))

library(parallel)
options(mc.cores = parallel::detectCores())

# load data
ASVs <- read.csv(file = here("data/ASV_table_26-JUN-2018.csv"), row.names = 1)
bTAX <- read.csv(file = here("data/bTAX_table_26-JUN-2018.csv"), row.names = 1)
mEL  <- read.csv(file = here("data/EL_sample_data_final.csv"), row.names = 1)
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

ELD <- merge(ASVs, mEL, by = "row.names")
NPD <- merge(ASVs, mNP, by = "row.names")

ELD[,'log_ratio'] <- log(ELD[,'bac']/ELD[,'host'])
NPD[,'log_ratio'] <- log(NPD[,'bac']/NPD[,'host'])

# construct models;
run_tot_mods <- function(dd, site, ...){

  mod.names <- c('ga0','ga0b','ga1','ga2','skn0','skn0b','skn1','skn2')

  # fit gaussian models
  ga0 <- brm(bf(log_ratio ~ (1|sp_id)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  ga0b <- brm(bf(log_ratio ~ (1|sp_id), sigma = ~ herb_dmg), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  ga1 <- brm(bf(log_ratio ~ herb_dmg + (1|sp_id)), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  ga2 <- brm(bf(log_ratio ~ herb_dmg + (1|sp_id), sigma = ~ herb_dmg), data = dd, family = gaussian, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  # fit skew-normal models
  skn0 <- brm(bf(log_ratio ~ (1|sp_id)), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  skn0b <- brm(bf(log_ratio ~ (1|sp_id), sigma = ~ herb_dmg), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  skn1 <- brm(bf(log_ratio ~ herb_dmg + (1|sp_id)), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  skn2 <- brm(bf(log_ratio ~ herb_dmg + (1|sp_id), sigma = ~ herb_dmg), data = dd, family = skew_normal, inits = "0", control = list(adapt_delta = 0.99), iter = 10000)

  mods <- list(ga0,ga0b,ga1,ga2,skn0,skn0b,skn1,skn2)

  names(mods) <- mod.names

  for(m in 1:length(mods)){
    #cat(sprintf("Saving model '%s' for Family %s\n", names(mods)[m], FAMS[s]))
    saveRDS(mods[[m]], file = paste0(here("models/all-bacteria/"),'all-bact-mod_',site,"_",names(mods)[m],'.rds'))
  }
  return(mods)
}

# run models, saving along the way:
EL_mods <- run_tot_mods(dd = ELD, site = 'EL')
NP_mods <- run_tot_mods(dd = NPD, site = 'NP')

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

  REP_TEMPLATE <- here("scripts/brmsfit_modrep_template_allBact.Rmd")

  # LOAD CFU MODEL
  CFU_MOD_FILE <- here("models/cfu_m1.rds")
  cfu_mod <- readRDS(file = CFU_MOD_FILE) # this part is hard-coded for now

  mod_names <- c('ga0','ga0b','ga1','ga2','skn0','skn0b','skn1','skn2')

  # first add lOO-IC to each of the models:
  for(i in seq_along(the_mods)){
    the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'loo', reloo = TRUE) # assign LOOIC to each model object in modlist
    #the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'waic') # assign wAIC to each
    assign(mod_names[i], the_mods[[i]]) # turn objects named in the_mods into objects with those names
  }

  # generate LOO objects for comparisons, below
  mods <- list(ga0,ga0b,ga1,ga2,skn0,skn0b,skn1,skn2)
  l1 <- loo(ga0)
  l2 <- loo(ga0b)
  l3 <- loo(ga1)
  l4 <- loo(ga2)
  l5 <- loo(skn0)
  l6 <- loo(skn0b)
  l7 <- loo(skn1)
  l8 <- loo(skn2)

  loo_list <- list(l1,l2,l3,l4,l5,l6,l7,l8) # put them all in one place

  loo1 <- do.call(rbind, sapply(loo_list, function(x) with(x, data.frame(model_name, looic, se_looic)), simplify = F))
  loo1$loo_id <- c('l1','l2','l3','l4','l5','l6','l7','l8') # assign IDs for referenceing LOO results

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
                                get(paste0(loo1[4,'loo_id'])),
                                get(paste0(loo1[5,'loo_id'])))

  #h1 <- table(unique(DAT[,c('Row.names','herb_dmg')])[,2])[[2]] # counts number of damaged leaves in dataset
  #h0 <- table(unique(DAT[,c('Row.names','herb_dmg')])[,2])[[1]] # counts number of un-damaged leaves in dataset
  h1 <- sum(DAT$herb_dmg)
  h0 <- length(DAT$herb_dmg) - h1

  # now, need to generate pp samples using best_mod and input data:
  sim_dat <- DAT[,c('herb_dmg','Row.names','sp_id')] # creates newdata df contaiing each bASV across both herb_dmg leaf types.
  names(sim_dat) <- c('herb_dmg','samples','sp_id') # assign col headers to newdata
  #sim_dat$samp_id <- c(1:(h0+h1)) # add arbitrary sample IDs to rows of newdata, where each sample is uniquely identified and contains all bASVs.
  #sim_dat$herb_dmg <- factor(sim_dat$herb_dmg)
  sim_dat$sp_id <- factor(sim_dat$sp_id)
  sim_dat_l <- split(sim_dat, sim_dat$samples) # split this up into a list based on samp_id and run each through our posterior predictive simulation, below.

  ## Simualte R y_rep datasets for each leaf (i.e., one 'draw' per joint posterior sample) and concatenate them all together
  # use sapply() and return a matrix. Need to supply a list as input.
  R <- 200 # first define number of posterior simulated datasets to draw

  pplr_by_samp_id <- lapply(sim_dat_l, FUN = leaf_level_pplr, the_mod = best_mod, nsamples = R) # this returns an (R*length(bASVs)) by length(samp_id) matrix; takes 1.5 minutes on MacBook Pro 2012

  # this results object needs to be added back to the original data.frame elements of input list
  # do this with a loop since the results pplr_by_samp_id isn't itself a list
  for(i in seq_along(sim_dat_l)){
    sim_dat_l[[i]] <- cbind(sim_dat_l[[i]], t(pplr_by_samp_id[[i]])) # need the transpose of the pplr element, since the brms::posterior_predict returns an R by n matrix
  }

  # now flatten this results list and melt it into long-format for summaries and plots:
  sim_dat_l2 <- do.call(rbind, sim_dat_l)
  sim_dat_l2 <- reshape2::melt(sim_dat_l2, id.vars = c('sp_id','herb_dmg', 'samples'),
                               variable.name = "rep",
                               value.name = 'ln_bac_host')

  sim_cfu_1 <- pp_cfu_per_leaf_all(sim_dat_l2, cfu_mod, samp_col = 'samples', the_fact = 'sp_id') # simulate posterior predicted CFU and log CFU values per sample, where each leaf per rep gets its own joint posterior draw from cfu_mod; takes 30 sec.

  # export for book-keeping:
  write.csv(sim_cfu_1, file = paste0(here("models/all-bacteria/"),"ppcfu_",SET,'.csv')) # export all simulation results for book-keeping; large file!! Better ways to do all this, for sure.

  # generate summary plots to put into report:
  sum_leaf_rep <- dplyr::group_by(sim_cfu_1, samples, herb_dmg, rep) %>% summarise(tot_leaf_pp_cfu = sum(pp_cfu))

  # next, need to calculate averages across damage types for each rep:
  sum_rep <- dplyr::group_by(sum_leaf_rep, herb_dmg, rep) %>% summarise(mu_tot_leaf_pp_cfu = mean(tot_leaf_pp_cfu),
                                                                        med_tot_leaf_pp_cfu = median(tot_leaf_pp_cfu),
                                                                        log_mu_tot_leaf_pp_cfu = log(mu_tot_leaf_pp_cfu,10),
                                                                        log_med_tot_leaf_pp_cfu = log(med_tot_leaf_pp_cfu,10))

  # save rep-level summaries for book-keeping and later plotting:
  write.csv(sum_rep, file = paste0(here("models/all-bacteria/"),"sum-all-rep_",SET,'.csv'))

  all_plot1 <- ggplot(sum_rep, aes(x = factor(herb_dmg), y = log_mu_tot_leaf_pp_cfu, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    theme_phy1() + ggtitle(paste0(SET)) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  all_plot1_med <- ggplot(sum_rep, aes(x = factor(herb_dmg), y = log_med_tot_leaf_pp_cfu, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted median log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    theme_phy1() + ggtitle(paste0(SET)) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # plot random rep to inspect:
  the_rep <- sample(unique(sim_cfu_1$rep),1)
  raw1 <- dplyr::filter(sim_cfu_1, rep == the_rep) %>% group_by(samples, herb_dmg) %>% summarise(tot_leaf_pp_cfu = sum(pp_cfu))
  plot_raw1 <- ggplot(raw1, aes(x = factor(herb_dmg), y = log(tot_leaf_pp_cfu,10), col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted log CFU per g leaf") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    theme_phy1() + ggtitle(paste0(SET)) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # second plot: average across the log-ratio values:
  lr_avg <- dplyr::filter(sim_cfu_1, rep == the_rep) %>% group_by(samples, herb_dmg) %>% summarise(mu_lr_leaf_rep = mean(ln_bac_host))
  lr_raw <- ggplot(lr_avg, aes(x = factor(herb_dmg), y = mu_lr_leaf_rep, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted log-ratio of 16S counts") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    theme_phy1() + ggtitle(paste0(SET)) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # plot distribution of median for all bacteria:
  all_sum <- dplyr::group_by(sum_rep, herb_dmg) %>% summarise(log_mu_mu = mean(log_mu_tot_leaf_pp_cfu),
                                                                  log_mu_0.025 = quantile(log_mu_tot_leaf_pp_cfu, 0.025),
                                                                  log_mu_0.975 = quantile(log_mu_tot_leaf_pp_cfu, 0.975),
                                                                  log_mu_0.50  = quantile(log_mu_tot_leaf_pp_cfu, 0.5),
                                                                  log_mu_0.75  = quantile(log_mu_tot_leaf_pp_cfu, 0.75),
                                                                  log_mu_0.25  = quantile(log_mu_tot_leaf_pp_cfu, 0.25),
                                                                  log_med_0.025 = quantile(log_med_tot_leaf_pp_cfu, 0.025),
                                                                  log_med_0.975 = quantile(log_med_tot_leaf_pp_cfu, 0.975),
                                                                  log_med_0.50  = quantile(log_med_tot_leaf_pp_cfu, 0.5),
                                                                  log_med_0.75  = quantile(log_med_tot_leaf_pp_cfu, 0.75),
                                                                  log_med_0.25  = quantile(log_med_tot_leaf_pp_cfu, 0.25))
  # save full summaries for book-keeping and later plotting:
  write.csv(all_sum, file = paste0(here("models/all-bacteria/"),"sum-all_",SET,'.csv'))

  all_plot2 <- ggplot(all_sum) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.75, ymax = log_mu_0.25, col = factor(herb_dmg)), size = 2) +
    geom_point(aes(x = factor(herb_dmg), y = log_mu_mu), col = "white", shape = 95, size = 8) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    #scale_y_continuous(limits = c(6,9), breaks = c(6,7,8,9)) +
    theme_phy1() + ggtitle(paste0(SET)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  all_plot2_med <- ggplot(all_sum) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_med_0.025, ymax = log_med_0.975, col = factor(herb_dmg)), size = 0.5) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_med_0.75, ymax = log_med_0.25, col = factor(herb_dmg)), size = 2) +
    geom_point(aes(x = factor(herb_dmg), y = log_med_0.50), col = "white", shape = 95, size = 8) +
    ylab("predicted median log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','gray40')) +
    #scale_y_continuous(limits = c(6,9), breaks = c(6,7,8,9)) +
    theme_phy1() + ggtitle(paste0(SET)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  #### GENERATE REPORT ####
  rmarkdown::render(input = paste0(REP_TEMPLATE),
                    output_format = "pdf_document",
                    output_file = paste0("modrep_",SET,"_",Sys.Date(),".pdf"),
                    output_dir = here("models/all-bacteria/"))
}

compare_tot_mods(EL_mods, DAT = ELD, SET = 'EL')
compare_tot_mods(NP_mods, DAT = NPD, SET = 'NP')

# now run script 'prep_ppcfu_rep-level.R` to generate figure output for log median across sample types as well as
# posterior difference in number of doublings between sample types
