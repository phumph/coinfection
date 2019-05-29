#!/usr/bin/env Rscript
# R script to automate report generation from brmsfit object, stored as .rds object.
# Last updated: PTH 18-JUL-18

# load required libraries
library(knitr)
library(markdown)
library(rmarkdown)
library(brms)
library(here)
library(latex2exp)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# args = commandArgs(trailingOnly=TRUE)
#
# if (length(args)<5) {
#   stop("Positional arguments are: [1] input dir where models are stored; [2] output dir where reports will be put; [3] Rmarkdown template file that generates the report; [4] leaf dataset; and [5] 1-column .csv containing Families to parse.", call.=FALSE)
# }

#IN_DIR <- args[1]
#IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/EL/'
# IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/NP/' # positional arg 1
IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/EL/' # positional arg 1

#OUT_DIR <- args[2]
# OUT_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/NP/reports/' # positional arg 2
OUT_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/EL/reports/' # positional arg 2

#REP_TEMPLATE <- args[3]
REP_TEMPLATE <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/scripts/brmsfit_modrep_template.Rmd' # positional arg 3

#DAT_FILE <- args[4]
# DAT_FILE <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv" # positional arg 4
DAT_FILE <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/ELD_long_v1.csv" # positional arg 4
DAT <- read.csv(DAT_FILE, header = T)
SET <- unlist(strsplit(DAT_FILE,"/"))
SET <- unlist(strsplit(SET[length(SET)],'_'))[1]

#FAMS <- args[5]
FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1) # positional arg 5

# edit this with latest version of CFU mod
CFU_MOD_FILE <- here("models/cfu_m1.rds")

# LOAD CFU MODEL
cfu_mod <- readRDS(file = CFU_MOD_FILE) # this part is hard-coded for now

#### FUNCTION DEFINITIONS ####
# all of these have been moved to phy_functions.R by PTH 2018-SEP-11

#### MAIN CALLS ####

# start by looping through each family in FAMS:
for (m in seq_along(FAMS)){
  the_mods <- grab_rds(FAMS[m], IN_DIR) # pull up .rds model results for a given family from IN_DIR
  mod_names <- paste0(sapply(fam_files, grab_mod_names)) # grab model names for easy assignment:

  # add LOOIC and WAIC to models
  for(i in seq_along(the_mods)){
    the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'loo') # assign LOOIC to each model object in modlist
    #the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'waic') # assign wAIC to each
    assign(mod_names[i], the_mods[[i]]) # turn objects named in the_mods into objects with those names
  }

  # generate LOO objects for comparisons, below
  l1 <- loo(ga0)
  l2 <- loo(ga1)
  l3 <- loo(ga2)
  l4 <- loo(skn0)
  l5 <- loo(skn1)
  l6 <- loo(skn2)
  l7 <- loo(skn3)
  l8 <- loo(skn4)

  loo_list <- list(l1,l2,l3,l4,l5,l6,l7,l8) # put them all in one place

  # capture relevant attributes of loo objects and put them into a data.frame
  # NOTE: use of looic and se_looic as attributes will be depreciated soon!! If this breaks in the future, this is likely why...
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

  #### ADDRESSING THE ZERO PROBLEM ####
  # Before we generate simulations at Family-level using best model for posterior predicted log CFU per leaf class, we have to decide whether missing values are true sampled zeros
  # or whether missing observations are probably there but below our limit of detection for a given bASV in a given sample.

  # We can take the original data, find the zeros, compute the detection limit, and score each (for each bASV) as to whether they should be included in the predicted distribution or not.
  # The idea is that if the model predicts log ratio values below this detection limit (of count == 1), then we might as well include the smaple since we have no information
  # that justifies its exclusion (i.e., the bASV is probably there but at too-low-to-detect densities).

  # In contrast, if a given sample contains fewer than a reasonable random draw from the log-ratio predictive distribution (in this case, below 0.01 of the predicted simulations),
  # then there's a good chance this leaf did not actually contain this bASV. This is a "true" zero, and we assume that the bASV wasn't in the leaf in the first place.

  # Of course, only on-average high abundance bASVs will offer enough resolution to score missing observations as true zeros,
  # so we should expect to see values below cutoff only for high abundance bASVs.
  # This is a built-in observational bias common to all taxon-counting problems across ecology
  # In this sense, this approach allows low-abundance bASVs to nonetheless be simulated in leaves for which the predicted counts overlap well with the observation limit.
  # These counts will not affect the overall counts much for the Family-level abundance, since the predicted counts are low in the first place (hence why they overlap ~1 count in the sample).
  # Instead, this approach "conservatively" reduces the imputed presence of high-abundance bASVs whose observational biases would have a larger effect on Family-level abundance predictions.

  dd <- dplyr::filter(DAT, Family == FAMS[m]) # grab focal family
  ddz <- dd[dd$bASV_count<1,] # find rows with zero-count
  ddz$zc <- log(1/ddz$host) # define cutoff for these samples. This cutoff is equivalent to sampling 1 bASV sequence count in the focal library.

  # determine which best_mod was fit, then grab relevant newdat parameters for posterior predictions:
  if(1 %in% grep("0", loo1$model_name)){
    newdat <- data.frame(unique(best_mod$data$bASV))
    names(newdat) <- c('bASV')
  } else {
    newdat <- expand.grid(unique(best_mod$data$bASV), unique(best_mod$data$herb_dmg))
    names(newdat) <- c('bASV','herb_dmg')
  }

  zsamps <- 1000 # define number of posterior samples from which to compute the cutoff

  # draw posterior samples from best_mod for determining per-sample cutoff
  # Done this way in one go, each row of newdata per nsamples shares parameter draws from the joint posterior
  # but this is appropriate for now since we're looking at the distribution of simulated values per sample.
  post_samps1 <- brms::posterior_predict(best_mod, newdata = newdat, nsamples = zsamps)
  post_df <- data.frame(cbind(newdat,t(post_samps1))) # annotate posterior draws from best_mod for indexing

  # calculate cutoffs for posterior draws
  # given level is 99%
  # again, do this differently depending on whether best_mod included herb_dmg as term or not:
  if(1 %in% grep("0", loo1$model_name)){
    post_df2 <- data.frame(bASV = post_df$bASV,
                           pcz = apply(post_df[,-c(1:2)], 1, function(x) quantile(x, probs = c(0.01))))
  } else {
    post_df2 <- data.frame(bASV = post_df$bASV,
                           herb_dmg = post_df$herb_dmg,
                           pcz = apply(post_df[,-c(1:2)], 1, function(x) quantile(x, probs = c(0.01))))
  }

  # run through each bASV PRESENT IN best_mod. First, grab the levels:
  bASVs <- unique(best_mod$data$bASV)

  # calculate number of leaves in each insect damage class:
  h1 <- table(unique(DAT[,c('Row.names','herb_dmg')])[,2])[[2]] # counts number of damaged leaves in dataset
  h0 <- table(unique(DAT[,c('Row.names','herb_dmg')])[,2])[[1]] # counts number of un-damaged leaves in dataset
  bASV.res <- data.frame() # initialize results df for zero determination

  # run through bASVs and calculates sample prevalence and corrected predicted prevalence
  for(i in 1:length(bASVs)){
    bz <- dplyr::filter(ddz, bASV == paste0(bASVs[i])) # subset zero-count df for target bASV

    # if this bASV in fact has no zeros at all:
    if(dim(bz)[1]==0){
      a <- 0.5 # define Beta distribution parameter alpha according to Jeffrey's prior (info below)
      b <- 0.5 # define Beta distribution parameter beta according to Jeffrey's prior (info below)
      df1 <- data.frame(bASV = bASVs[i],
                        h0, h1, n0 = h0, n1 = h1, p0 = h0, p1 = h1, pp0_a = h0+a, pp0_b = b, pp1_a = h1+a, pp1_b = b) # put this all in one place
      bASV.res <- rbind(bASV.res, df1) # add to growing results df
      next # move to next bASV; nothing else to calculate.
    } else {
      n0 <- h0 - length(bz$herb_dmg[bz$herb_dmg == '0']) # calculate sample prevalence in un-damaged leaves
      n1 <- h1 - length(bz$herb_dmg[bz$herb_dmg == '1']) # calculate sample prevalence in damaged leaves
      post1 <- dplyr::filter(post_df2, bASV == paste0(bASVs[i])) # sub-set posterior data-frame for target bASV

      # if best_mod did not include herb_dmg, just evaluate cutoff as though leaves of each type were interchangeable
      if(1 %in% grep("0", loo1$model_name)){
        zc_res <- table(bz$herb_dmg, bz$zc <= post1$pcz) # evaluate zeros against posterior 1% cutoffs and define results
      } else {
        zc_res <- table(bz$herb_dmg, bz$zc <= post1$pcz[match(bz$herb_dmg,post1$herb_dmg)]) # evaluate zeros against posterior 1% cutoffs and define results
      }

      # here's how to interpret this table, zc_res:
      # rows are herb_dmg == 0,1
      # cols reflect whether sample cutoff is below 1% of the posterior predictive distribution for leaves in a given class
      # TRUE means need to keep as "true" zero, since cutoff is below detection limit.
      # FALSE means these zeros are probably a detection limit thing and should be included in sums of Family-level abundance (i.e., sums of bASV predicted abundance per leaf).
      # so, for each leaf class, we sum up the number of FALSES and add this quantity to the sample prevalence estimate to get predicted prevalence:
      # this "adds back" false zeros.
      p0 <- n0 + sum(zc_res[1,1])
      p1 <- n1 + sum(zc_res[2,1])

      # let's use these sampling-corrected prevalence estimates to generate posterior estimates of Beta distribution parameters alpha and beta for posterior prevalence.
      # We'll use Jeffrey's (J's) prior as the conjugate prior, Beta(0.5,0.5), for estimating posterior prevalence in each leaf damage class.
      # This choice is largely for simplicity's sake. Choices here won't impact much, as other sources of variation will swamp this (i.e., in the log_ratio --> logCFU predictions).
      # One could try Bayes' (Beta(1,1)) prior, or generate empirical priors based on the overall relationship between average abundance and (zero-corrected) prevalence, after iterating from a non-informative prior.
      # That is, models here could be iterated, using Jeffrey's prior to arrive at a "naive" posterior, then adjusting this posterior based on the posterior distribution of prevalences given data on abundance.
      # But this seems overly complex and not necessarily more justifiable.
      # Moving on...
      a <- 0.5 # J's prior for alpha of the prior Beta distribution
      b <- 0.5 # J's prior for beta of the same
      pp0_a <- a + p0 # number of "successes" in un-damaged class, for posterior sampling
      pp0_b <- b + (h0 - p0) # number of "failures" in un-damaged class
      pp1_a <- a + p1 # number of "successes" in damaged class
      pp1_b <- b + (h1 - p1) # number of "failures" in damaged class
      # now we can use these parameters to draw from Beta(a+s, b+f), i.e., the posterior predicted prevalence estimate.
      # These draws will be used for each damage class when constructing simulated data as input to the predicted log CFU model.
      # Basically, this spread out the zeros across the simulations and allows for uncertainty in presence/absence of bASVs from leaf samples when simulating abundances.

      # export this data for focal bASV:
      df1 <- data.frame(bASV = bASVs[i],
                        h0,h1,n0,n1,p0,p1,pp0_a,pp0_b,pp1_a,pp1_b)
      bASV.res <- rbind(bASV.res, df1) # add to growing results list
    }
  }

  write.csv(bASV.res, file = paste0(IN_DIR,"ppcfu_data/prev-ests_",FAMS[m],'.csv')) # export prevalence estimates for safe-keeping

  #### POSTERIOR PREDICTIONS ####
  ## construct bASV-wise posterior predictions for CFU.
  # This should be straight-forward: basically, for predicted infection intensities, we don't need to care about prevalence.
  # This is because these predicted distributions are implicitly conditional on having sampled the bASV in a given leaf type (damages or un-damaged).
  # We just need to grab distributions of estimates from both damage classes to calculate means (or medians, and relevant quantiles).
  # Let's take R draws from the posterior of the log_ratio model (best_mod) and use this as input (i.e., simulated predictor values) to calculate predicted log CFU.
  # Each leaf sampled will contain all bASVs for the Family and will get a unique draw from the joint posterior of the CFU mod, per rep in R reps.
  # That is, each leaf needs to represent a unique "sample" from the posterior, and each simulated "dataset" (of which there are R) consists of h0 un-damaged leaves and h1 damaged leaves.
  # Each leaf in h0,h1 gets its own posterior draw PER REP. Meaning, this process generates (h0+h1)*R draws from the joint posterior of the best_mod as well as the cfu_mod.

  # Once this is done, we'll have a population of estimates for each bASV and we can generate summary statistics for plotting and comparisons.
  # Following this, we will draw again but estimate a presence/absence for each simulated leaf, for determining how to sum predicted abundances across bASV for the Family-wide estimates.
  # Basicaly, this weights the average contribution of a bASV to the Family-level total by the probability of its occurence within a leaf of a given type (i.e., group-level prevanece estimates).

  sim_dat <- expand.grid(c(rep(0,h0),rep(1,h1)), bASVs) # creates newdata df contaiing each bASV across both herb_dmg leaf types.
  names(sim_dat) <- c('herb_dmg', 'bASV') # assign col headers to newdata
  sim_dat$samp_id <- rep(c(1:(h0+h1)),length(bASVs)) # add arbitrary sample IDs to rows of newdata, where each sample is uniquely identified and contains all bASVs.
  sim_dat_l <- split(sim_dat, sim_dat$samp_id) # split this up into a list based on samp_id and run each through our posterior predictive simulation, below.

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
  sim_dat_l2 <- reshape2::melt(sim_dat_l2, id.vars = c('bASV','herb_dmg', 'samp_id'),
                               variable.name = "rep",
                               value.name = 'ln_bac_host')

  sim_cfu_1 <- pp_cfu_per_leaf(sim_dat_l2, cfu_mod) # simulate posterior predicted CFU and log CFU values per sample, where each leaf per rep gets its own joint posterior draw from cfu_mod; takes 30 sec.

  #### ADDING PREVALENCE ESTIMATES ####
  # Next step is to weight each bASV simulated count by prevalence. For each rep of the simulated data above (sim_cfu_1), draw a 1/0 Bernoulli random variable.
  # In this case, the probability p of success for the Bernoulli draw will be drawn from the posterior Beta distribution given, for each leaf class, for each bASV, as stored in bASV.res
  # We will generate a vector of bernoulli random variables based on posterior draws from this prevalence Beta disributions.
  # This draws one trial per leaf per bASV per sim dataset per rep.
  # first we need to re-order sim_cfu_1 so we can simply cbind the output of sim_bernoulli.
  # output format is: bASV., samp_id, rep, Currently, order of sim_cfu_1 is samp_id,rep,bASV.
  sim_cfu_1 <- dplyr::arrange(sim_cfu_1, bASV, rep, samp_id)
  # also arrange bASV.res by bASV to create same alphabetical arbitrary order as above:
  bASV.res <- dplyr::arrange(bASV.res, bASV)

  # finally, do the bernoulli draws per sample/rep for each bASV, and append results to sim_cfu_1
  # this set of functions applies sim_bernoulli to each row of bASV.res, taking as four positional arguments the Beta distribution parameters
  # and returns a vectorized output the same length of nrows(sim_cfu_1).
  sim_cfu_1$hits <- c(with(bASV.res, mapply(sim_bernoulli, pp0_a, pp0_b, pp1_a, pp1_b)))

  # # CHECKING CODE
  # # Let's see if average prevalence reflects, more or less, the predicted prevalence:
  # prev_check <- dplyr::group_by(sim_cfu_1, bASV, herb_dmg) %>% summarise(prev = sum(hits)/length(hits))
  # prev_check <- reshape2::dcast(prev_check, bASV ~ herb_dmg)
  # # add posterior mean prevalence to bASV.res for comparison to simulated draws
  # bASV.res <- dplyr::mutate(bASV.res, pmu_0 = pp0_a / (pp0_a+pp0_b), pmu_1 = pp1_a / (pp1_a+pp1_b))
  # prev_check2 <- merge(prev_check, bASV.res[,c('bASV','pmu_0','pmu_1')], by = 'bASV', sort = F)
  # names(prev_check2)[c(2,3)] <- c('p0','p1')
  # # ggplot(prev_check2, aes(x = pmu_1, y = p1)) + geom_point() + geom_abline(intercept = 0, slope = 1) # plot to examine correspondence.
  # # OK this checks out. The average prevalences are quite close to the predicted means of the respective posterior beta distributions.

  # export sim_cfu_1 for completeness:
  write.csv(sim_cfu_1, file = paste0(IN_DIR,"ppcfu_data/ppcfu_",FAMS[m],'.csv')) # export all simulation results for book-keeping; large file!! Better ways to do all this, for sure.


  #### GENERATING PLOTS ####
  ## PART 1: Family-level analysis, summing counts across all bASVs for the focal Family.

  # sum predicted CFU counts across bASVs per leaf sample, and rep; chuck out leaves per rep with hits==0 (i.e., predicted prevalence was 0 for this leaf-rep combo.).
  fam_sum_leaf_rep <- dplyr::filter(sim_cfu_1, hits == 1) %>% group_by(samp_id, herb_dmg, rep) %>% summarise(tot_leaf_pp_cfu = sum(pp_cfu))

  # next, need to calculate averages across damage types for each rep:
  fam_sum_rep <- dplyr::group_by(fam_sum_leaf_rep, herb_dmg, rep) %>% summarise(mu_tot_leaf_pp_cfu = mean(tot_leaf_pp_cfu),
                                                                                med_tot_leaf_pp_cfu = median(tot_leaf_pp_cfu),
                                                                                log_mu_tot_leaf_pp_cfu = log(mu_tot_leaf_pp_cfu,10),
                                                                                log_med_tot_leaf_pp_cfu = log(med_tot_leaf_pp_cfu,10),
                                                                                prev = length(tot_leaf_pp_cfu))

  # we will output these data so that we may summarise them for each Family all at once
  fam_sum_rep$Family <- FAMS[m] # add Family name column to results

  # output
  write.csv(fam_sum_rep, file = paste0(IN_DIR,"ppcfu_data/fam-sum-rep_",FAMS[m],'.csv'))

  # we can make a quick plot here to check out distributions are sensible
  fam_plot1 <- ggplot(fam_sum_rep, aes(x = factor(herb_dmg), y = log_mu_tot_leaf_pp_cfu, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # same plot but with medians as summary statistic to plot
  fam_plot1_med <- ggplot(fam_sum_rep, aes(x = factor(herb_dmg), y = log_med_tot_leaf_pp_cfu, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted median log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # finally, plot the distribution of the raw simualted cfu data (not summarised) as boxplots from a random rep:
  # first generate subsetted dataframe containing 1 random rep and plot the predicted/simulated cfu data from all leaves in each damage class.
  the_rep <- sample(unique(sim_cfu_1$rep),1)
  fam_raw1 <- dplyr::filter(sim_cfu_1, hits == 1, rep == the_rep) %>% group_by(samp_id, herb_dmg) %>% summarise(tot_leaf_pp_cfu = sum(pp_cfu))
  fam_plot_raw1 <- ggplot(fam_raw1, aes(x = factor(herb_dmg), y = log(tot_leaf_pp_cfu,10), col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted log CFU per leaf") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # for completeness, let's generate plot of the simulted log-ratio values, averaged across bASVs, across damage types:
  fam_lr_avg <- dplyr::filter(sim_cfu_1, hits == 1, rep == the_rep) %>% group_by(samp_id, herb_dmg) %>% summarise(mu_lr_leaf_rep = mean(ln_bac_host))
  fam_lr_raw <- ggplot(fam_lr_avg, aes(x = factor(herb_dmg), y = mu_lr_leaf_rep, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted log-ratio of 16S counts") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # Next we generate Family-level summary statistics on the distribution of log(E[cfu]) per damage class.
  # To be clear, here we are calculating the mean of the log of the predicted mean CFU counts and displaying credible posterior predicted intervals around these metrics.
  # This is for outputting to table (via xtable) in the model report and for generating all-Family plot panel in Fig. 2.
  fam_sum <- dplyr::group_by(fam_sum_rep, herb_dmg) %>% summarise(log_mu_mu = mean(log_mu_tot_leaf_pp_cfu),
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

  fam_sum$Family <- FAMS[m] # add Family tag to output

  # export these summary statistics so we can collate them all later and plot each Family all together in Fig. 2.
  write.csv(fam_sum, file = paste0(IN_DIR,"ppcfu_data/fam-sum-all_",FAMS[m],'.csv'))

  # plot another way:
  fam_plot2 <- ggplot(fam_sum) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.75, ymax = log_mu_0.25, col = factor(herb_dmg)), size = 2) +
    geom_point(aes(x = factor(herb_dmg), y = log_mu_mu), col = "white", shape = 95, size = 8) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    #scale_y_continuous(limits = c(6,9), breaks = c(6,7,8,9)) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  ## PART 2: bASV-level plotting and summary-statistics
  # make rep-level summary of bASV log(E[cfu])
  bASV_rep <- dplyr::filter(sim_cfu_1, hits == 1) %>% group_by(bASV, herb_dmg, rep) %>% summarise(mu_pp_cfu = mean(pp_cfu),
                                                                                                  med_pp_cfu = median(pp_cfu),
                                                                                                  log_mu_pp_cfu = log(mu_pp_cfu,10),
                                                                                                  log_med_pp_cfu = log(med_pp_cfu,10),
                                                                                                  prev = length(pp_cfu)) %>% arrange(desc(log_mu_pp_cfu))

  # now generate summary statistics for log(E[cfu]) for each bASV for each leaf damage class:
  bASV_sum <- dplyr::group_by(bASV_rep, bASV, herb_dmg) %>% summarise(log_mu_mu = mean(log_mu_pp_cfu),
                                                                      log_mu_0.025 = quantile(log_mu_pp_cfu, 0.025),
                                                                      log_mu_0.975 = quantile(log_mu_pp_cfu, 0.975),
                                                                      log_mu_0.50  = quantile(log_mu_pp_cfu, 0.5),
                                                                      log_mu_0.75  = quantile(log_mu_pp_cfu, 0.75),
                                                                      log_mu_0.25  = quantile(log_mu_pp_cfu, 0.25),
                                                                      log_med_0.025 = quantile(log_med_pp_cfu, 0.025),
                                                                      log_med_0.975 = quantile(log_med_pp_cfu, 0.975),
                                                                      log_med_0.50  = quantile(log_med_pp_cfu, 0.5),
                                                                      log_med_0.75  = quantile(log_med_pp_cfu, 0.75),
                                                                      log_med_0.25  = quantile(log_med_pp_cfu, 0.25)) %>% arrange(desc(log_mu_mu))

  # export both files
  write.csv(bASV_rep, file = paste0(IN_DIR,"ppcfu_data/bASV-rep_",FAMS[m],'.csv'))
  write.csv(bASV_sum, file = paste0(IN_DIR,"ppcfu_data/bASV-all_",FAMS[m],'.csv'))

  # re-level rep-level bASV for plotting from highest to lowest predicted CFU
  bASV_rep$bASV <- factor(bASV_rep$bASV, levels = paste0(unique(bASV_sum$bASV)))
  bASV_sum$bASV <- factor(bASV_sum$bASV, levels = paste0(unique(bASV_sum$bASV)))

  # plot distribution of means per bASV across leaf types:
  bASV_plot1 <- ggplot(bASV_rep, aes(x = bASV, y = log_mu_pp_cfu, col = factor(herb_dmg))) +
    #geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    #scale_y_continuous(limits = c(4.5,10), breaks = c(5,6,7,8,9,10)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("log mean pp_cfu")

  # now plot in our fancy minimalist line style:
  bASV_plot2 <- ggplot(bASV_sum) +
    geom_linerange(aes(x = bASV, ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5, position = position_dodge(1)) +
    geom_linerange(aes(x = bASV, ymin = log_mu_0.25, ymax = log_mu_0.75, col = factor(herb_dmg)), size = 1.25, position = position_dodge(1)) +
    geom_point(aes(x = bASV, y = log_mu_mu, group = factor(herb_dmg)), size = 6, shape = 95, col = "white", position = position_dodge(1)) +
    theme_phy1() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "top") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    ylab("log mean pp_cfu")

  #### GENERATE REPORT ####
  rmarkdown::render(input = paste0(REP_TEMPLATE),
                    output_format = "pdf_document",
                    output_file = paste0("modrep_",FAMS[m],"_",Sys.Date(), ".pdf"),
                    output_dir = paste0(OUT_DIR))
}
