#! /usr/bin/Rscript
# Rscript to generate sub-plot level predicted bacterial abundance
# last updated 2018-SEP-18 by PTH

# for EL, since we don't have even sampling across subplots, we will have to compute an average and go with that. Calculate it per-stem, and report this as density of bacteria per stem at the patch level.
# This will have to be compared to herbivore density as well, also calculated as number of mined leaves per stem. Clearly, patches with more leaves have the potential to be larger reservoirs for bacteria. So this aspect is not a problem.

library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# load dataset to supply the newdata for predictive simulations
FOCAL.DAT.NAME <- here("data/EL.field.hormone.exp.master.txt")
focal.dat <- read.table(FOCAL.DAT.NAME,T,"\t")
names(focal.dat)[3] <- 'condition'

focal.dat$subplot <- with(focal.dat, paste0(plot,'_',tx))
focal.dat$plant <- with(focal.dat, paste0(sub.plot, stem.num))
focal.dat$condition <- relevel(focal.dat$condition, ref = 'CTR')

#ggplot(subplot_data, aes(x = n.mined.leaves.sum)) + geom_density()

# load model to sample coefficient posteriors
FOCAL.MOD.NAME <- here("models/EL/Pseudomonadaceae_mod_skn2.rds") # best model for EL
focal.mod <- readRDS(FOCAL.MOD.NAME)

#### MAIN CALLS ####

# 1. generate newdat template
nd1 <- focal.dat[,c('plot','subplot','plant','n.leaves','n.mined.leaves')]
nd1$bASV <- factor('Pseudomonas_3')
nd1 <- nd1[complete.cases(nd1),]

# keep tabs on distribution of counts:
# cpre <- table(table(nd1$subplot))

# brute-force it since sapply won't handle data.frame col names very well...:
res <- data.frame()
for (k in 1:nrow(nd1)){
  res <- rbind(res,row_rep(nd1[k,]))
  #print(k)
}
res$bASV <- factor(res$bASV)

nd2_l <- split(res, res$sample_id) # transform into list

# extract posterior distributions of coefficients from focal.mod
coef_post <- data.frame(as.matrix(focal.mod))

# this model for EL does not have different slope for Pseudomonas_3
# grab slope coefficients as vectors
a0 <- coef_post[,'b_Intercept'] # overall intercept of herb_dmg==0
b0 <- coef_post[,'b_herb_dmg1'] # overall slope effect of herb_dmg1==1
alpha <- coef_post[,'alpha'] # defines skew-normal shape parameter for estimating residual error
sigma0 <- coef_post[,'b_sigma_Intercept'] # defines posterior of standard deviation of residual error for herm_dmg==0
sigma1 <- coef_post[,'b_sigma_herb_dmg1'] # defines posterior of standard deviation of residual error for herm_dmg==1
a1 <- coef_post[,grep('r_bASV.Pseudomonas_3.Intercept.',names(coef_post),value = T)] # focal bASV-specific additive intercept
#b1 <- coef_post[,grep(paste0(nd2_l[[1]]$bASV,'\\.'),names(coef_post),value = T)[2]] # focal bASV-specific additive slope

# calculate vector of posterior predicted log-ratio when herb_dmg==0
non_dmg_pplr <- (a0+a1) + rskew_normal(n = length(sigma0), mu = 0, sigma = sigma0, alpha = alpha)
yes_dmg_pplr <- (a0+a1) + b0 + rskew_normal(n = length(sigma1), mu = 0, sigma = sigma1, alpha = alpha)

# now sample these simulated posterior log-ratio values depending on whether leaf is damaged or not:
leaf_level_pplr_precomputed <- function(x, no, yes, ...){
  i <- sample(1:length(no), size = R, replace = TRUE)
  if(x$herb_dmg == 0) {
    y_rep <- no[i]
  } else {
    y_rep <- yes[i]
  }
  return(y_rep)
}

# now add log_ratio factor to nd2:
R <- 200
nd2_l2 <- lapply(nd2_l, function(x) leaf_level_pplr_precomputed(x, no = non_dmg_pplr, yes = yes_dmg_pplr, R))

# combine with original data.frame
for(i in seq_along(nd2_l)){
  nd2_l[[i]] <- cbind(nd2_l[[i]], t(nd2_l2[[i]])) # need the transpose of the pplr element, since the brms::posterior_predict returns an R by n matrix
}

# now flatten this results list and melt it into long-format for summaries and plots:
nd3 <- do.call(rbind, nd2_l) # this take a while.. not sure why.
nd3 <- reshape2::melt(nd3, id.vars = c('bASV','herb_dmg', 'sample_id','plot','subplot','plant','n.mined.leaves','leaf_id','plant_id'),
                      variable.name = "rep",
                      value.name = 'ln_bac_host')
write.csv(nd3, file = here("models/EL/patch_level_Pseudomonas_pplr.csv"), quote = F)
rm(nd2_l)

# next step is to transform pplr into log cfu for each leaf:
cfu_mod <- readRDS(here("models/cfu_m1.rds")) # load relevant CFU model
jp <- data.frame(cfu_mod)
n <- dim(nd3)[1]
nd4 <- cbind(nd3, jp[sample(1:length(jp[,1]), n, replace=T),]) # draw from joint posterior
nd4$deviate <- sapply(nd4$sigma, function(x) rnorm(n = 1, mean = 0, sd = x)) # draw normal deviates for each row:
nd4$log_pp_cfu <- nd4$b_Intercept + nd4$b_ln_bac_host * nd4$ln_bac_host + nd4$deviate # now calculate pp_cfu from model coefficients
nd4$pp_cfu <- 10^nd4$log_pp_cfu # exponentiate y_rep values
rm(nd3)

### generate prevalence estimates for each leaf.. and keep this the same for all reps so it introduced noise at the sample level not rep level
## go through each sample_id and draw, according to herb_dmg status, from the posterior beta distribution given in the prevalence estimates file:
# load prevalence estimates:
prev_ests <- read.csv(here("models/EL/ppcfu_data/prev-ests_Pseudomonadaceae.csv"))
prev_ests <- prev_ests[prev_ests$bASV == 'Pseudomonas_3',]

# use beta distribution to draw p for drawing subsequent Bernoulli random variable for each of the leaves given by sample_id
the_samples <- unique(nd4[,c('sample_id','herb_dmg')])

# break into damaged and un-damaged; rbind; make the calls; bind back to nd4:
the_samples_no  <- the_samples[the_samples$herb_dmg==0,]
the_samples_yes <- the_samples[the_samples$herb_dmg==1,]
the_samples_no$pr  <- rbeta(n = length(the_samples_no[,1]), shape1 = prev_ests$pp0_a, shape2 = prev_ests$pp0_b)
the_samples_yes$pr <- rbeta(n = length(the_samples_yes[,1]), shape1 = prev_ests$pp1_a, shape2 = prev_ests$pp1_b)
the_samples_all <- rbind(the_samples_no,the_samples_yes)

# now generate bernoulli draws for each element of pr:
the_samples_all$hit <- bernoulli_trials(the_samples_all$pr)

# add back to nd4:
nd5 <- merge(nd4, the_samples_all[,c('sample_id','hit')], by = 'sample_id', sort = F)
nd5$pp_cfu <- nd5$pp_cfu * nd5$hit
rm(nd4)

# save to disc
write.csv(nd5, file = here("models/EL/patch_level_Pseudomonas_ppcfu.csv"), quote = F)

#### Summaries and plots ####

# read in data, if not already loaded:
nd5 <- read.csv(file = here("models/EL/patch_level_Pseudomonas_ppcfu.csv"))

# first, sum up per subplot per rep:
nd5_sp_rep_plant <- dplyr::group_by(nd5, subplot, plant, rep) %>% summarise(sp_sum_pp_cfu = sum(pp_cfu),
                                                                            log_sp_sum_pp_cfu = log(sp_sum_pp_cfu, 10),
                                                                            n.leaves.tot = length(sample_id),
                                                                            herb_dmg.sum = sum(herb_dmg))

# at this point, we have R rows of each plant from each subplot. The number of plants varies across subplots with an average of 14.5.
# we need to sum up across stems and compute an average for each subplot. We can then also average the number of mined leaves.
nd5_sp_rep <- dplyr::group_by(nd5_sp_rep_plant, subplot, rep) %>% summarise(sp_mean_pp_cfu = mean(sp_sum_pp_cfu),
                                                                            log_sp_mean_pp_cfu = log(sp_mean_pp_cfu, 10),
                                                                            n.leaves.sp = length(n.leaves.tot),
                                                                            herb_dmg.sp = mean(herb_dmg.sum),
                                                                            n.stems = length(n.leaves.tot)) %>% arrange(herb_dmg.sp)


# now we will perform an extrapolation to bring the n.stems up to n=16 per subplot. This is so we can plot the data
# on the same scale as the NP dataset. All we do is multiple sp_mean_pp_cfu and herb_dmg.sp by a constant (16) to get the new columns.
nstems <- 16
nd5_sp_rep <- dplyr::mutate(nd5_sp_rep,
                            sp_mean_pp_cfu.ext = round(nstems*sp_mean_pp_cfu),
                            herb_dmg.sp.ext = round(nstems*herb_dmg.sp),
                            log_sp_mean_pp_cfu.ext = log(sp_mean_pp_cfu.ext,10))

# plot this to check what's up:
nd5_sp_rep$subplot <- factor(nd5_sp_rep$subplot, levels = paste0(unique(nd5_sp_rep$subplot)))

# ggplot(nd5_sp_rep, aes(x = herb_dmg.sp, y = log_sp_mean_pp_cfu)) + geom_point()

# now capture quantiles of this quantity among the reps:
nd5_sp_summary <- dplyr::group_by(nd5_sp_rep, subplot) %>% summarise(log_med_0.025 = quantile(log_sp_mean_pp_cfu.ext, 0.025),
                                                                     log_med_0.975 = quantile(log_sp_mean_pp_cfu.ext, 0.975),
                                                                     log_med_0.50  = quantile(log_sp_mean_pp_cfu.ext, 0.5),
                                                                     log_med_0.75  = quantile(log_sp_mean_pp_cfu.ext, 0.75),
                                                                     log_med_0.25  = quantile(log_sp_mean_pp_cfu.ext, 0.25),
                                                                     n.leaves.sp = unique(n.leaves.sp),
                                                                     herb_dmg.sp.ext = unique(herb_dmg.sp.ext),
                                                                     n.stems = unique(n.stems)) %>% arrange(herb_dmg.sp.ext)


# merge back with other info from subplots:
# sp_dat <- unique(nd5_sp_rep[,c('plot','subplot','n.mined.leaves')]) %>% group_by(subplot) %>%
#   summarise(sum.n.mined.leaves = sum(n.mined.leaves))

# clean up workspace
rm(nd5)

# nd5_sp_summary_all <- merge(nd5_sp_summary, sp_dat, by = 'subplot', sort = F) %>% arrange(desc(sum.n.mined.leaves))

#### plots ####
nd5_sp_summary$subplot <- factor(nd5_sp_summary$subplot, levels = paste0(nd5_sp_summary$subplot))

co_infection_p1 <- ggplot(nd5_sp_summary, aes(x = herb_dmg.sp.ext, y = log_med_0.50)) +
  geom_linerange(aes(ymin = log_med_0.025, ymax = log_med_0.975), size = 0.33, col = "gray80") +
  geom_linerange(aes(ymin = log_med_0.25, ymax = log_med_0.75), size = 0.75, col = "gray70") +
  geom_point(alpha = 0.5, col = "gray20") +
  xlab("mined leaves per patch") +
  ylab("pp. bacterial load per patch") +
  scale_y_continuous(limits = c(6.4,8.5), breaks = seq(6.5,8.5,0.5)) +
  scale_x_continuous(limits = c(0,110), breaks = seq(0,100,20)) #+
  #annotation_logticks(sides = 'l')
#co_infection_p1

# distribution of herbivory per sub-plot
dp1 <- ggplot(nd5_sp_summary, aes(x = herb_dmg.sp.ext)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  scale_x_continuous(limits = c(0,110), breaks = seq(0,100,20))

# distribution of median ppcfu per patch:
dp3 <- ggplot(nd5_sp_summary, aes(x = log_med_0.50)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  scale_x_continuous(limits = c(6.4,8.5), breaks = seq(6.5,8.5,0.5)) +
  coord_flip()

ggarrange(plotlist = list(dp1, dp1, #placeholder
                          co_infection_p1, dp3),
          widths = c(1,0.5),
          heights = c(0.5, 1),
          align = 'hv',
          ncol = 2, nrow = 2) %>% ggsave(filename = here("figs/Fig3_EL_v2.pdf"), width = 3, height = 2.5)

# now make plot of quantile of herbivore damage versus percent contribution to bacterial propagules:
herb_rank <- sort(nd5_sp_summary$herb_dmg.sp.ext)
nd5_sp_summary$herbivory_rank <- sapply(nd5_sp_summary$herb_dmg.sp.ext, function(x) { (length(herb_rank[herb_rank > x])/length(herb_rank)) })

# now calculate proportion of total propagules (median):
sum_propagules <- sum(10^nd5_sp_summary$log_med_0.50)
#nd5_sp_summary_all$propagule_fraction <- sapply(nd5_sp_summary_all$log_med_0.50, function(x) x^10/sum_propagules)

# generate cumulative sum of propagules as one increases herbivory_rank
cumsum_herb <- nd5_sp_summary[,c('herbivory_rank','log_med_0.50')] %>% arrange(herbivory_rank)
cumsum_herb$csum_props <- cumsum(10^cumsum_herb$log_med_0.50)/sum(10^cumsum_herb$log_med_0.50)

# compile useful quantiles:
focal_quantiles <- c(0.1,0.2,0.5,0.8)
csum_prop_qs <- sapply(focal_quantiles, function(x) last(cumsum_herb$csum_props[cumsum_herb$herbivory_rank<=x]))
d1 <- data.frame(focal_quantiles = focal_quantiles*100,
                 csum_prop_qs = round(csum_prop_qs*100))

# text for labels:
d2 <- data.frame(xlabs = paste0(d1[,1],'%'),
                 ylabs = paste0(d1[,2],'%'),
                 focal_quantiles = d1$focal_quantiles,
                 csum_prop_qs = d1$csum_prop_qs)

scale_plot1 <- ggplot() +
  geom_segment(data = d1, aes(x = focal_quantiles, y = 3.0, xend = focal_quantiles, yend = csum_prop_qs), lty = 2, lwd = 0.25, col = "darkorange2") +
  geom_segment(data = d1, aes(x = 3.5, y = csum_prop_qs, xend = focal_quantiles, yend = csum_prop_qs), lty = 2, lwd = 0.25, col = "darkorange2") +
  geom_text(data = d2, aes(x = focal_quantiles, y = 0, label = xlabs), size = 2.5, col = 'gray40') +
  geom_text(data = d2, aes(x = 0, y = csum_prop_qs, label = ylabs), size = 2.5, col = 'gray40') +
  geom_line(data = cumsum_herb, aes(x = herbivory_rank*100, y = csum_props*100), lwd = 1, col = "gray40") +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  xlab("% rank of patch-level herbivore load (high to low)") +
  ylab("% of total population bacterial load") +
  geom_segment(aes(x = 0, xend = 100, y = 0, yend = 100), col = 'gray80', lty = 3)

ggsave(scale_plot1, filename = here("figs/EL_scale_plot1_v2.pdf"), width = 2.5, height = 2)

