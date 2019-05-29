# #!/usr/bin/env Rscript


# load required libraries
library(knitr)
library(markdown)
library(rmarkdown)
library(brms)
library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# args = commandArgs(trailingOnly=TRUE)
#
# if (length(args)<5) {
#   stop("Positional arguments are: [1] input dir where models are stored; [2] output dir where reports will be put; [3] Rmarkdown template file that generates the report; [4] prevalence/abundance stats file' [5] OTU table", call.=FALSE)
# }

## R script to automate report generation from brmsfit object, stored as .rds object.

#IN_DIR <- args[1]
#IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/EL/'
IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/NP/'
#OUT_DIR <- args[2]
OUT_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/NP/reports/'
#REP_TEMPLATE <- args[3]
REP_TEMPLATE <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/scripts/brmsfit_modrep_template.Rmd'
DAT_FILE <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv"
DAT <- read.csv(DAT_FILE, header = T)
SET <- unlist(strsplit(DAT_FILE,"/"))
SET <- unlist(strsplit(SET[length(SET)],'_'))[1]
FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1)

# edit this with latest version of CFU mod
CFU_MOD_FILE <- here("models/cfu_m1.rds")

# LOAD CFU MODEL
cfu_mod <- readRDS(file = CFU_MOD_FILE)

#### FUNCTION DEFINITIONS ####
grab_rds <- function(...){
  fam_files <- Sys.glob(paste0(IN_DIR, FAMS[m],'*'))
  fam_mods <- list()
  for (f in 1:length(fam_files)){
    fam_mods[[f]] <- readRDS(file = fam_files[f])
  }
  return(fam_mods)
}

grab_mod_names <- function(x){
  #gsub('.rds','',unlist(strsplit(gsub(IN_DIR,'',x),'_'))[3])
  gsub('.rds','',grep('.rds',unlist(strsplit(gsub(IN_DIR,'',x),'_')), value = T))
}

bernoulli_trials<-function(p)
{
  U<-runif(length(p),0,1) # draw p numbers from uniform distribution
  outcomes <- U<p # determine if draw is less than p, taken as a vector (i.e., element-wise comparison between elements of U and p). If the uniform number is less than p, score as success.
  return(as.numeric(outcomes)) # return binary vector of same length as p.
}

sim_bernoulli <- function(a0,b0,a1,b1,...){
  # generate vector of probabilities to use in bernoulli trials
  prob_vec <- rep(c(rbeta(n = h0, shape1 = a0, shape2 = b0),rbeta(n = h1, shape1 = a1, shape2 = b1)),R)
  # now conduct bernoulli trials, one for each element of prob_vec
  berns <- bernoulli_trials(prob_vec)
  return(berns)
}

# function def to generate predicted CFUs from draws from joint posterior of cfu_mod; returns predicted CFU on linear (i.e. count) scale.
cfu_post_pred.old <- function(x, IDs, log=FALSE, ...){
  l <- length(unique(IDs))
  jp <- data.frame(cfu_mod)
  cfu_jpd <- jp[sample(1:length(jp[,1]),l,replace=T),] # draw from joint posterior for each row of the
  cfu_jpd2 <- cfu_jpd[IDs,] # replicate each draw from joint posterior across all sub-items for each leaf sample
  cfu <- cfu_jpd[,1] + cfu_jpd[,2] * x + rnorm(n = 1, mean = 0, sd = cfu_jpd[,3])

  ifelse(log==TRUE,return(10^cfu),return(cfu))
}




# let's look at the joint posterior:
# cfu_jpd <- data.frame(t(as.matrix(cfu_mod)[sample(1:4000,1000),]))
# ggplot(data.frame(t(cfu_jpd[,1:1000])), aes(x = b_Intercept, y = b_ln_bac_host)) + geom_point(alpha = 0.1) + geom_density_2d()

#### MAIN CALLS ####
# start by generating a loop for each family:
for (m in 1:length(FAMS)){
#for (m in 7:length(FAMS)){
  the_mods <- grab_rds(FAMS[m], IN_DIR) # pull up .rds model results for a given family from IN_DIR
  mod_names <- paste0(sapply(fam_files, grab_mod_names)) # grab model names for easy assignment:

  # add LOOIC and WAIC to models
  for(i in 1:length(the_mods)){
  #for(i in 1:3){
    the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'loo') # assign LOOIC to each model object in modlist
    #the_mods[[i]] <- add_ic(the_mods[[i]], ic = 'waic')
    assign(mod_names[i], the_mods[[i]]) # turn objects named in the_mods into objects with those names
  }

  # capture objects themselves in named list
  #models <- sapply(mod_names, get, simplify = FALSE)

  l1 <- loo(ga0)
  l2 <- loo(ga1)
  l3 <- loo(ga2)
  l4 <- loo(skn0)
  l5 <- loo(skn1)
  l6 <- loo(skn2)
  l7 <- loo(skn3)
  l8 <- loo(skn4)

  loo_list <- list(l1,l2,l3,l4,l5,l6,l7,l8)

  # loo_comps <- compare_ic(l1,l2,l3,l4,l5,l6,l7,l8)
  # loo_comps2 <- as.data.frame(do.call(rbind,loo_comps))

  # use of looic and se_looic as attributes will be depreciated soon!! If this breaks in the future, this is likely why...
  loo1 <- do.call(rbind, sapply(loo_list, function(x) with(x, data.frame(model_name, looic, se_looic)), simplify = F))
  loo1$loo_id <- c('l1','l2','l3','l4','l5','l6','l7','l8')

  # need to assign looic objects:
  for(l in 1:length(loo_list)){
    assign(loo1$loo_id[l], loo_list[[l]])
  }

  loo1 <- loo1[order(loo1$looic),]

  # define best model based on lowest LOOIC
  best_mod <- get(paste0(loo1$model_name[1]))

  # grab top three; compare
  loo_comp1 <- compare_ic(get(paste0(loo1[1,'loo_id'])),
                          get(paste0(loo1[2,'loo_id'])),
                          get(paste0(loo1[3,'loo_id'])),
                          get(paste0(loo1[4,'loo_id'])),
                          get(paste0(loo1[5,'loo_id'])))

  ## Now generate simulations at Family-level using best model for posterior predicted log CFU per leaf class:
  # first let's look at the zero problem.
  # We can take the original data, find the zeros, compute the detection limit, and score each (for each bASV) as to whether they should be included in the predicted distribution or not.
  # The idea is that if the model predicts log ratio values below this detection limit (of count == 1), then we might as well include the smaple since we have no information
  # that justifies its exclusion (i.e., the bASV is probably there but at too-low-to-detect densities).
  # In contrast, if a given sample contains fewer than a reasonable random draw from the log-ratio predictive distribution (in this case, below 0.01 of the predicted simulations),
  # then there's a good chance this leaf did not actually contain this bASV.
  # In this sense, this approach uses well-sampled bASVs to exclude counts from leaves they likely were absent from, while allowing low-abundance bASVs to nonetheless be simulated
  # in leaves for which the predicted counts overlap well with the observation limit. These counts will not affect the overall counts much for the Family-level abundance.

  dd <- dplyr::filter(DAT, Family == FAMS[m]) # grab focal family
  ddz <- dd[dd$bASV_count<1,] # find rows with zero-count
  ddz$zc <- log(1/ddz$host) # define cutoff for these samples. This cutoff is equivalent to sampling 1 bASV sequence count in the library.

  # posterior-predict distribution of log ratios for each bASV for both leaf types form best_mod
  # determine which best_mod was fit, then grab relevant newdat parameters for posterior predictions:

  if(1 %in% grep("0", loo1$model_name)){
    newdat <- data.frame(unique(best_mod$data$bASV))
    names(newdat) <- c('bASV')
  } else {
    newdat <- expand.grid(unique(best_mod$data$bASV), unique(best_mod$data$herb_dmg))
    names(newdat) <- c('bASV','herb_dmg')
  }
  zsamps <- 1000
  # newdat2 <- newdat[rep(seq.int(1,nrow(newdat)),reps),]

  post_samps1 <- brms::posterior_predict(best_mod, newdata = newdat, nsamples = zsamps)

  # annotate posterior draws from best_mod for indexing
  post_df <- data.frame(cbind(newdat,t(post_samps1)))

  # generate plots of distributions of posterior predicted log-ratios (from random sim); aveaged across all bASVs. This is basically the marginal effect of herb_dmg on Y
  post_df_long <- reshape2::melt(post_df, id.vars = c('bASV','herb_dmg'), value.name = 'log_ratio', variable.name = 'rep')
  post_df_long$rep <- sapply(post_df_long$rep, function(x) gsub('X','',x))
  post_df_long$rep <- as.numeric(post_df_long$rep)
  post_df_long2 <- group_by(post_df_long, herb_dmg, rep) %>% summarise(mu_lr = mean(log_ratio))
  the_reps <- sample(unique(post_df_long$rep),100)

  lr_p1 <- ggplot() + geom_jitter(data = post_df_long2[post_df_long2$rep %in% the_reps,], aes(x = factor(herb_dmg), y = mu_lr, col = factor(herb_dmg)), width = 0.15, alpha = 0.3) +
    geom_boxplot(data = post_df_long2, aes(x = factor(herb_dmg), y = mu_lr, col = factor(herb_dmg)), alpha = 0.3) + xlab("DMG") + ylab("E[log-ratio]") +
    scale_color_manual(values = c("gray40","dodgerblue")) +
    theme_phy1() + theme(legend.position = "none")
  #lr_p1

  # now plot the distribution of means per rep, for each bASV:
  post_df_long_bASV <- group_by(post_df_long, bASV, herb_dmg, rep) %>% summarise(mu_lr = mean(log_ratio)) %>% group_by(bASV, herb_dmg) %>%
    summarise(mu_mu_lr = mean(mu_lr),
              mu_mu_lr_0.025 = quantile(mu_lr, 0.025),
              mu_mu_lr_0.975 = quantile(mu_lr, 0.975),
              mu_mu_lr_0.50  = quantile(mu_lr, 0.5),
              mu_mu_lr_0.75  = quantile(mu_lr, 0.75),
              mu_mu_lr_0.25  = quantile(mu_lr, 0.25)) %>% arrange(desc(mu_mu_lr))

  post_df_long_bASV$bASV <- factor(post_df_long_bASV$bASV, levels = paste0(unique(post_df_long_bASV$bASV)))

  bASV_plot_lr <- ggplot(post_df_long_bASV) +
    geom_linerange(aes(x = bASV, ymin = mu_mu_lr_0.025, ymax = mu_mu_lr_0.975, col = factor(herb_dmg)), size = 0.5, position = position_dodge(1)) +
    geom_linerange(aes(x = bASV, ymin = mu_mu_lr_0.25, ymax = mu_mu_lr_0.75, col = factor(herb_dmg)), size = 1.25, position = position_dodge(1)) +
    geom_point(aes(x = bASV, y = mu_mu_lr, col = factor(herb_dmg)), size = 6, shape = 95, col = "white", position = position_dodge(1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_color_manual(values = c('gray40','dodgerblue')) +
    theme(legend.position = "top")
  #bASV_plot_lr

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
  bASV.res <- data.frame() # initialize results df

  # run through bASVs and calculates sample prevalence, predicted prevalence
  for(i in 1:length(bASVs)){
    bz <- dplyr::filter(ddz, bASV == paste0(bASVs[i])) # subset zero-count df for target bASV

    # if this bASV in fact has no zeros:
    if(dim(bz)[1]==0){
      a <- 0.6
      b <- 0.5
      df1 <- data.frame(bASV = bASVs[i],
                        h0, h1, n0 = h0, n1 = h1, p0 = h0, p1 = h1, pp0_a = h0+a, pp0_b = h0+b, pp1_a = h1+a, pp1_b = h1+b)
      bASV.res <- rbind(bASV.res, df1)
      next
    } else {

      n0 <- h0 - length(bz$herb_dmg[bz$herb_dmg == '0']) # calculate sample prevalence in un-damaged leaves
      n1 <- h1 - length(bz$herb_dmg[bz$herb_dmg == '1']) # calculate sample prevalence in damaged leaves
      post1 <- dplyr::filter(post_df2, bASV == paste0(bASVs[i])) # sub-set posterior data-frame for target bASV

      # if best_mod did not include herb_dmg, just evaluate cutoff as though leaves were intechangeable
      if(1 %in% grep("0", loo1$model_name)){
        zc_res <- table(bz$herb_dmg, bz$zc <= post1$pcz) # evaluate zeros against posterior 1% cutoffs and define results
      } else {
        zc_res <- table(bz$herb_dmg, bz$zc <= post1$pcz[match(bz$herb_dmg,post1$herb_dmg)]) # evaluate zeros against posterior 1% cutoffs and define results
      }

      # here's how to interpret this table, zc_res:
      # rows are herb_dmg == 0,1
      # cols are sample cutoff is below 1% of the posterior predictive distribution for that leaf class
      # TRUE means need to keep as "true" zero
      # FALSE means these zeros are probably a detection limit thing and should be included in sums of Family-level abundance (i.e., sums of bASV predicted abundance per leaf).
      # so, for each leaf class, we sum up the number of FALSES and add this quantity to the sample prevalence estimate to get predicted prevalence:
      p0 <- n0 + sum(zc_res[1,1])
      p1 <- n1 + sum(zc_res[2,1])
      # let's use these sampling-corrected prevalence estimates to generate posterior estimates of Beta distribution parameters alpha and beta for posterior prevalence
      # We'll use Jeffrey's (J's) prior as the conjugate Beta(0.5,0.5) for estimating posterior prevalence in each leaf damage class.
      # This choice is largely for simplicity's sake. Choices here won't impact much, as other sources of variation will swamp this (i.e., in the log_ratio --> logCFU predictions).
      # One could try Bayes' (Beta(1,1)) prior, or generate empirical priors based on the overall relationship between average abundance and (zero-corrected) prevalence.
      # Models here could be iterated, using Jeffrey's prior to arrive at a "naive" posterior, then adjusting this posterior based on the posterior distribution of prevalences given data on abundance.
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
      # export this data for focal bASV:
      df1 <- data.frame(bASV = bASVs[i],
                        h0,h1,n0,n1,p0,p1,pp0_a,pp0_b,pp1_a,pp1_b)
      bASV.res <- rbind(bASV.res, df1) # add to growing results list
    }
  }

  # export prevalence estimates
  write.csv(bASV.res, file = paste0(IN_DIR,"ppcfu_data/prev-ests_",FAMS[m],'.csv'))

  #### POSTERIOR PREDICTIONS ####
  ## construct bASV-wise posterior predictions for CFU
    # This should be straight-forward: basically, for predicted infection intensities, we don't need to care about prevalence.
    # This is because these predicted distributions are implicitly conditional on having sampled the bASV in a given leaf type (damages or un-damaged).
    # We just need to grab distribution of estimates from both damage classes to calculate means (or medians, and relevant quantiles).
    # Let's take S draws from the posterior of the log_ratio model (cfu_mod) and calculate predicted log CFU based on draw from a single row of the joint posterior of the CFU mod.
    # Then we'll have a population of estimates for each bASV and we can generate summary statistics for plotting from this.
    # Following this, we will draw again but estimate a presence/absence for each simulated leaf, for summing across bASV for the Family-wide estimates.
    # Basicaly, this weights the average contribution of a bASV to the Family-level total by the probability of its occurence within a leaf of a given type (i.e., group-level prevanece estimates).

  # sim_dat <- expand.grid(c(0,1), bASVs) # replicate the input data in terms of number of leaves of each class
  sim_dat <- expand.grid(c(rep(0,h0),rep(1,h1)), bASVs) # replicate this an arbitrary number of times to calculate means per class per rep dataset
  names(sim_dat) <- c('herb_dmg', 'bASV') # assign col headers to newdata
  sim_dat$samp_id <- rep(c(1:(h0+h1)),length(bASVs)) # add arbitrary sample IDs to rows of newdata

  # split this up into a list based on samp_id and run each through our posterior predictive simulation:
  sim_dat_l <- split(sim_dat, sim_dat$samp_id)

  ## Simualte R y_rep datasets for each leaf (i.e., one 'draw' per joint posterior sample) and concatenate them all together
  # use sapply() and return a matrix. Need to supply a list as input
  R <- 200 # first define number of posterior simulated datasets to draw

  # draw log_ratio posterior samples from best_mod
  leaf_level_pplr <- function(x, ...){
    y_rep <- brms::posterior_predict(best_mod, newdata = sim_dat_l[[1]], nsamples = R)
    return(y_rep)
  }

  pplr_by_samp_id <- lapply(sim_dat_l, FUN = leaf_level_pplr, nsamples = R) # this returns an (R*length(bASVs)) by length(samp_id) matrix
  # this results object needs to be added back to the original data.frame elements of input list
  # do this with a loop since the results pplr_by_samp_id isn't itself a list
  for(i in seq_along(sim_dat_l)){
    sim_dat_l[[i]] <- cbind(sim_dat_l[[i]], t(pplr_by_samp_id[[i]])) # need the transpose of the pplr element, since the brms::posterior_predict returns an R by n matrix
  }

  # now flatten this results list and melt it into long-format for summaries and plots:
  sim_dat_l2 <- do.call(rbind, sim_dat_l)
  sim_dat_l3 <- reshape2::melt(sim_dat_l2, id.vars = c('bASV','herb_dmg', 'samp_id'),
                               variable.name = "rep",
                               value.name = 'ln_bac_host')

  # adding pp_cfu using tapply:
  cfu_post_pred_lowest <- function(x, log=FALSE, ...){
    #l <- length(unique(IDs))
    jp <- data.frame(cfu_mod)
    cfu_jpd <- jp[sample(1:length(jp[,1]),1,replace=T),] # draw from joint posterior
    deviates <- sapply(cfu_jpd[,3], function(s) rnorm(n = length(x), mean = 0, sd = s))
    cfu <- cfu_jpd[,1] + cfu_jpd[,2] * x + deviates

    ifelse(log==FALSE,return(10^cfu),return(cfu))
  }


  # function to generate joint posterior draws per samp_id for each rep, with observation-level normal deviates per sigma draw per sample
  # this matrix will be of the same dimension as input matrix
  pp_cfu_per_leaf <- function(x, ...){
    # for each samp_id and rep, generate
    # make samp_id draws
    jp <- data.frame(cfu_mod)
    dat1 <- expand.grid(samp_id = c(1:max(as.numeric(x$samp_id))),
                        rep = c(1:max(as.numeric(x$rep))))
    n <- dim(dat1)[1]
    cfu_jpd <- cbind(dat1, jp[sample(1:length(jp[,1]),n,replace=T),]) # draw from joint posterior
    # now add to original data.frame
    the_rows <- rep(c(1:dim(cfu_jpd)[1]), each = length(unique(x$bASV)))
    x2 <- cbind(x, cfu_jpd[the_rows,-c(1,2)])
    # now draw normal deviates for each row:
    x2$deviate <- sapply(x2$sigma, function(x) rnorm(n = 1, mean = 0, sd = x))
    # now calculate pp_cfu
    x2$log_pp_cfu <- x2$b_Intercept + x2$b_ln_bac_host * x2$ln_bac_host + x2$deviate
    x2$pp_cfu <- 10^x2$log_pp_cfu
    return(x2)
  }

  sim_cfu_1 <- pp_cfu_per_leaf(sim_dat_l3, cfu_mod)

  # let's quickly look at summaries
  # let's do Family-level sums
  sim_cfu_fam <- dplyr::group_by(sim_cfu_1, samp_id, rep, herb_dmg) %>% summarise(tot = sum(pp_cfu)) %>% group_by(herb_dmg, rep) %>% summarise(mu = mean(tot))

  ggplot(sim_cfu_fam, aes(x = factor(herb_dmg), y = log(mu,10))) + geom_jitter(width = 0.15, alpha = 0.2) + geom_violin(alpha = 0.2)
  # THAT's more like what I'm looking for...



  # # look at how these simulated values fall out:
  # sim_dat0 <- reshape2::melt(cbind(sim_dat, t(best_mod_post)),
  #                            id.vars = c('bASV','herb_dmg', 'samp_id'),
  #                            variable.name = "rep",
  #                            value.name = 'ln_bac_host')
  #
  # # average across samp_id for each level of herb_dmg:
  # sim_dat0b <- dplyr::group_by(sim_dat0, bASV, herb_dmg, rep) %>% summarise(mu_lr = mean(ln_bac_host))
  # sdt1 <- ggplot(sim_dat0b[sim_dat0b$bASV == 'Pseudomonas_3',], aes(x = factor(herb_dmg), y = mu_lr, group = rep)) + geom_point(alpha = 0.3) + geom_line(alpha = 0.3) +
  #   scale_y_continuous(limits = c(-5,3))
  #
  # # the question is whether each newdata instance gets single draw. I suspect this is the case. let's average across samples and then plot these to find out:
  # sim_dat0c <- dplyr::group_by(sim_dat0, bASV, samp_id, herb_dmg) %>% summarise(mu_lr = mean(ln_bac_host))
  # sdt2 <- ggplot(sim_dat0c[sim_dat0c$bASV == 'Pseudomonas_3',], aes(x = factor(herb_dmg), y = mu_lr, group = samp_id)) + geom_point(alpha = 0.3) + geom_line(alpha = 0.3) +
  #   scale_y_continuous(limits = c(-5,3))
  # ggarrange(plotlist = list(sdt1, sdt2), ncol = 2)
  #
  # # the difference is super clear. The plot on the left represents the range of the DATA simulated, while the right plot shows the distribution of the summary statistic (i.e., what we are interested in plottig).
  # # So, we need to draw once from the jpd for each simulated sample.
  #
  # # now let's make draws from cfu_mod joint posterior on sim_dat0 to simulate log_cfu
  # xxx <- brms::posterior_predict(cfu_mod, newdata = data.frame(ln_bac_host = sim_dat0$ln_bac_host))

  # we have to make sure that we are in fact drawing from the posterior in the way that we expect.
  # main question is whether nsamples corresponds to a single draw from the jpd for each row of newdata, or whether each datapoint gets its own draw.


  # send these samples through cfu_post_pred
  # have to do this PER SAMPLE, not per rep:
  sim_dat2 <- reshape2::melt(cbind(sim_dat, apply(best_mod_post, 1, cfu_post_pred, IDs = sim_dat$samp_id)),
                             id.vars = c('bASV','herb_dmg', 'samp_id'),
                             variable.name = "rep",
                             value.name = 'pp_log_cfu')

  # # OK this isn't producing the result I expect. Try by hand a few rows:
  # test.df2 <- cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[1,], IDs = sim_dat$samp_id))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[2,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[3,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[4,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[5,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[6,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[7,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[8,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[9,], IDs = sim_dat$samp_id)))
  # test.df2 <- rbind(test.df2,
  #                   cbind(sim_dat,pp_cfu = cfu_post_pred(x = best_mod_post[10,], IDs = sim_dat$samp_id)))
  # test.df2$rep <- rep(c(1:10),each = length(best_mod_post[1,]))

  # # now look at means..
  # test.df3 <- dplyr::group_by(test.df2, herb_dmg, rep) %>% summarise(tot = sum(pp_cfu)) %>% group_by(herb_dmg, rep) %>% summarise(mu_tot = mean(tot))
  # ggplot(test.df3, aes(x = factor(herb_dmg), y = log(mu_tot,10), group = rep))+geom_point()+geom_line()
  #

  # export for reproducibility
  write.csv(sim_dat2, file = paste0(IN_DIR,"/ppcfu_data/ppcfu_bASV-level_",FAMS[m],'.csv'))
#
#   # try to figure out whether this is actually working
#   # this averages across samples (samp_id) for each rep for each damage class. The distribution of these values should represent the variation in the summary statistic for this bASV
#   testdat1 <- dplyr::filter(sim_dat2, bASV == 'Pseudomonas_3') %>% group_by(herb_dmg, rep) %>% summarise(mu_log = mean(pp_log_cfu))
#   ggplot(testdat1, aes(x = factor(herb_dmg), y = mu_log, group = rep)) + geom_point() + geom_line()
#
#
#   # do the same with medians, which is less susceptible to Jensen's inequality than means:
#   testdat2 <- dplyr::filter(sim_dat2, bASV == 'Pseudomonas_3') %>% group_by(herb_dmg, rep) %>% summarise(med = median(pp_cfu))
#   ggplot(testdat2, aes(x = factor(herb_dmg), y = log(med,10), group = rep)) + geom_point() + geom_line()


  # generate plot of distribution of CFU per leaf class (not means):
  # lr_plot1 <- ggplot()
  # choose a random rep
  # the_rep <- as.numeric(sample(unique(sim_dat2$rep),1))
  # fam_plot1_raw <- ggplot(sim_dat2[sim_dat2$rep==the_rep,], aes(x = factor(herb_dmg), y = pp_cfu)) + geom_jitter(width = 0.15) + geom_boxplot() + theme_phy1()

  # test_dat1 <- dplyr::filter(sim_dat2, rep == the_rep, bASV == 'Pseudomonas_3')
  #
  # # look at rep 5
  # test_dat1 <- dplyr::filter(sim_dat2, rep == 5, bASV == 'Pseudomonas_3')
  #
  # tp1 <- ggplot(test_dat1, aes(x = factor(herb_dmg), y = log(pp_cfu,10))) + geom_jitter(width = 0.15) + geom_boxplot() + theme_phy1()
  #
  # # now group by rep and take the mean predicted CFU count per damage class
  # test_dat2 <- dplyr::filter(sim_dat2, bASV == 'Pseudomonas_3') %>% group_by(herb_dmg, rep) %>% summarise(tot = mean(pp_cfu))
  #
  # # test_dat2 <- dplyr::group_by(test_dat1, herb_dmg) %>% summarise(tot = mean(pp_cfu))
  #
  # tp2 <- ggplot(test_dat2, aes(x = factor(herb_dmg), y = log(tot,10))) + geom_jitter(width = 0.15) + geom_boxplot() + theme_phy1()
  # tp3 <- ggplot(test_dat2, aes(x = factor(herb_dmg), y = log(tot,10))) + geom_point() + theme_phy1() + geom_line(aes(group = rep), alpha = 0.1) + geom_boxplot()
  #
  # # need to plot the distribution of the differences, on the log scale. This will show the distribution of effect sizes of herbivore damage


  # calculate mean and median in each sample class per rep
  sim_dat3 <- dplyr::group_by(sim_dat2, bASV, herb_dmg, rep) %>% summarise(mu_pp_cfu = mean(pp_cfu),
                                                                           med_pp_cfu = median(pp_cfu))

  # now calculate mean among reps of replicated sample means, as well as 95% posterior predictive credible intervals
  sim_dat4 <- dplyr::group_by(sim_dat3, bASV, herb_dmg) %>% summarise(log_mu_pp_cfu = mean(log(mu_pp_cfu,10)),
                                                                      log_mu_0.025 = quantile(log(mu_pp_cfu,10), 0.025),
                                                                      log_mu_0.975 = quantile(log(mu_pp_cfu,10), 0.975),
                                                                      log_mu_0.50  = quantile(log(mu_pp_cfu,10), 0.5),
                                                                      log_mu_0.75  = quantile(log(mu_pp_cfu,10), 0.75),
                                                                      log_mu_0.25  = quantile(log(mu_pp_cfu,10), 0.25),
                                                                      log_med_pp_cfu = median(log(med_pp_cfu,10)),
                                                                      log_med_0.025 = quantile(log(med_pp_cfu,10), 0.025),
                                                                      log_med_0.975 = quantile(log(med_pp_cfu,10), 0.975),
                                                                      log_med_0.50  = quantile(log(med_pp_cfu,10), 0.5),
                                                                      log_med_0.75  = quantile(log(med_pp_cfu,10), 0.75),
                                                                      log_med_0.25  = quantile(log(med_pp_cfu,10), 0.25)) %>%
    arrange(desc(log_mu_pp_cfu))

  # Next step is to weight each bASV simulated count by prevalence. For each rep of the simulated data above (sim_dat2), draw a 1/0 Bernoulli random variable.
  # In this case, the probability p of success for the Bernoulli draw will be drawn from the posterior Beta distribution given, for each leaf class, for each bASV, as stored in bASV.res

  # generate vector of bernoulli random variables based on posterior draws from prevalence Beta disributions.
  # This draws one trial per leaf per bASV per sim dataset per rep. That's a lot of rows.
  sim_dat2$hits <- c(apply(X = bASV.res[,-c(1)], MARGIN = 1, FUN = sim_bernoulli, a0 = pp0_a, b0 = pp0_b, a1 = pp1_a, b1 = pp1_b))

  # take this data frame and capture only the simulated prevalence hits;
  # sum predicted CFU counts across bASVs per leaf sample, and rep;
  fam_sum_leaf_rep <- dplyr::filter(sim_dat2, hits == 1) %>% group_by(samp_id, herb_dmg, rep) %>% summarise(tot_leaf_pp_cfu = sum(pp_cfu))

  ############# ADDRESS THIS POINT  #############
  # decision time. Do these averages sum across all leaves "sampled" or is this basically just infection intensity?
  # another option is to simply assign predicted cfu to zero and include these in the mean calculation among leaves. This seems to make less sense to me.
  #############

  # proceeding after having removed bASV-samp_id rows with predicted hit of 0.
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
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # same plot but with medians as summary statistic to plot
  fam_plot1_med <- ggplot(fam_sum_rep, aes(x = factor(herb_dmg), y = log_med_tot_leaf_pp_cfu, col = factor(herb_dmg))) +
    geom_jitter(width = 0.15, alpha = 0.2) +
    geom_boxplot(alpha = 0.2) +
    ylab("predicted median log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # we are calculating the mean of the log of the predicted mean CFU counts and displaying credible posterior predicted intervals around these metrics:
  # this is for outputting to table (via xtable) in the model report:
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

  # plot another way:
  fam_plot2 <- ggplot(fam_sum) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5) +
    geom_linerange(aes(x = factor(herb_dmg), ymin = log_mu_0.75, ymax = log_mu_0.25, col = factor(herb_dmg)), size = 1.25) +
    geom_point(aes(x = factor(herb_dmg), y = log_mu_mu), col = "white", shape = 95, size = 8) +
    ylab("predicted mean log CFU") + xlab("herbivory") +
    scale_color_manual(values = c('gray40','dodgerblue')) +
    theme_phy1() + ggtitle(paste0(FAMS[m])) +
    theme(legend.position = "none",
          plot.title = element_text(size = 6))

  # let's try with position=position_dodge(1)
  # fam_plot3 <- ggplot(fam_sum) +
  #   geom_linerange(aes(x = 1, ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5, position=position_dodge(0.1)) +
  #   geom_linerange(aes(x = 1, ymin = log_mu_0.75, ymax = log_mu_0.25, col = factor(herb_dmg)), size = 1.25, position=position_dodge(0.1)) +
  #   geom_point(aes(x = 1, y = log_mu_mu, fill = factor(herb_dmg)), col = "white", shape = '|', size = 6, position=position_dodge(0.1)) +
  #   ylab("predicted mean log CFU") + xlab("herbivory") +
  #   scale_color_manual(values = c('gray40','dodgerblue')) +
  #   theme_phy1() + ggtitle(paste0(FAMS[m])) +
  #   scale_x_continuous(limits = c(0.90,1.4)) +
  #   theme(legend.position = "none",
  #         plot.title = element_text(size = 6),
  #         axis.text.y = element_blank(),
  #         axis.line.y = element_blank()) + coord_flip()
  # fam_plot3
  # Now we need to construct newdata for posterior predictions. This is a three-step process:
    # 1. Using posterior prevalence estimates for each bASV, create a df where each is represented in each damage class according to the draw from the prevalence posteriors
    # 2. Use this df as input to posterior_predict from brms to generate predicted y_rep values for log_ratio
    # 3. Use this simulated log_ratio as input to posterior_predict for the CFU model.
    # 4. Finally, with these data, sum up the exponentiated predicted y_rep_cfu values per leaf by Family, and also by bASV individally. Store all of this info in .csv for later plotting

  # # generate newdat for model checking:
  # nsamp <- 100
  # newdat <- brms::posterior_predict(best_mod, nsamples = nsamp)
  # #this file has lots of rows.. some for each bASV. Need to compile this into a reasonable format
  #
  # # load up function to perform posterior simulations for Family-level analysis
  # # for debugging:
  # mod <- cfu_mod
  # S <- 10
  # nd <- nd1
  # ppsim_fam_cfu <- function(mod, nd, ppstat = "mean", ppfact = "herb_dmg", S){
  #   require(dplyr)
  #   # define factor levels
  #   ppfact_levels <- levels(nd[,paste0(ppfact)])
  #   # create results matrix of correct dimensions
  #   Sres <- data.frame()
  #   # perform simulations and summary statisics
  #   for (s in 1:S){
  #     yrep <- posterior_predict(mod, newdat = nd, nsamples = 1) # generate samples
  #     Sres <- rbind(Sres,data.frame(dplyr::group_by(data.frame(yrep = t(yrep), ppfact = nd[,ppfact]), ppfact) %>% summarise(yrep = mean(yrep)), s))
  #   }
  #
  #   yrep <- data.frame(ppfact = nd[,ppfact],yrep = t(yrep), s = s, grp = 'sim')
  #   Sres[,'grp'] <- 'means'
  #
  #   return(list(Y0 = yrep, Y1 = Sres))
  #   #return(rbind(Sres,yrep))
  # }


  # re-level bASV for plotting along x-axis
  sim_dat4$bASV <- factor(sim_dat4$bASV, levels = paste0(unique(sim_dat4$bASV)))

  # generate plot of individual bASV predicted log mean CFUs
  bASV_plot1 <- ggplot(sim_dat4) +
    geom_linerange(aes(x = bASV, ymin = log_mu_0.025, ymax = log_mu_0.975, col = factor(herb_dmg)), size = 0.5, position = position_dodge(1)) +
    geom_linerange(aes(x = bASV, ymin = log_mu_0.25, ymax = log_mu_0.75, col = factor(herb_dmg)), size = 1.25, position = position_dodge(1)) +
    geom_point(aes(x = bASV, y = log_mu_pp_cfu, col = factor(herb_dmg)), size = 6, shape = 95, col = "white", position = position_dodge(1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_color_manual(values = c('gray40','dodgerblue')) +
    theme(legend.position = "top")


  #### GENERATE REPORT ####
  rmarkdown::render(input = paste0(REP_TEMPLATE),
                    output_format = "pdf_document",
                    output_file = paste0("modrep_",FAMS[m],"_",Sys.Date(), ".pdf"),
                    output_dir = paste0(OUT_DIR))
}
