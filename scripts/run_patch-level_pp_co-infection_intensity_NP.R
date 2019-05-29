#! /usr/bin/Rscript
# Rscript to generate sub-plot level predicted bacterial abundance
# last updated 2018-SEP-18 by PTH

library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# load dataset to supply the newdata for predictive simulations
FOCAL.DAT.NAME <- here("data/112213.npd.txt")
focal.dat <- read.table(FOCAL.DAT.NAME,T,"\t")
focal.dat <- dplyr::filter(focal.dat, plot != '58b')
focal.dat$condition <- relevel(focal.dat$condition, ref = 'MOCK')

# load model to sample coefficient posteriors
FOCAL.MOD.NAME <- here("models/NP/Pseudomonadaceae_mod_skn4.rds")
focal.mod <- readRDS(FOCAL.MOD.NAME)

#### MAIN CALLS ####

# 1. generate newdat template
nd1 <- focal.dat[,c('plot','subplot','plant','n.leaves','n.mined.leaves','n.fruits')]
nd1$bASV <- factor('Pseudomonas_3')

# brute-force it since sapply won't handle data.frame col names very well...:
res <- data.frame()
for (k in 1:nrow(nd1)){
  res <- rbind(res,row_rep(nd1[k,]))
}
res$bASV <- factor(res$bASV)
nd2_l <- split(res, res$sample_id) # transform into list
#nd2 <- res

# extract posterior distributions of coefficients from focal.mod
coef_post <- data.frame(as.matrix(focal.mod))

# grab slope coefficients as vectors
a0 <- coef_post[,'b_Intercept'] # overall intercept of herb_dmg==0
b0 <- coef_post[,'b_herb_dmg1'] # overall slope effect of herb_dmg1==1
alpha <- coef_post[,'alpha'] # defines skew-normal shape parameter for estimating residual error
sigma0 <- coef_post[,'b_sigma_Intercept'] # defines posterior of standard deviation of residual error for herm_dmg==0
sigma1 <- coef_post[,'b_sigma_herb_dmg1'] # defines posterior of standard deviation of residual error for herm_dmg==1
a1 <- coef_post[,grep(paste0(nd2_l[[1]]$bASV,'\\.'),names(coef_post),value = T)[1]] # focal bASV-specific additive intercept
b1 <- coef_post[,grep(paste0(nd2_l[[1]]$bASV,'\\.'),names(coef_post),value = T)[2]] # focal bASV-specific additive slope

# calculate vector of posterior predicted log-ratio when herb_dmg==0
non_dmg_pplr <- (a0+a1) + (0) + rskew_normal(n = length(sigma0), mu = 0, sigma = sigma0, alpha = alpha)
yes_dmg_pplr <- (a0+a1) + (b0 + b1) + rskew_normal(n = length(sigma1), mu = 0, sigma = sigma1, alpha = alpha)

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
nd3 <- do.call(rbind, nd2_l)
nd3 <- reshape2::melt(nd3, id.vars = c('bASV','herb_dmg', 'sample_id','plot','subplot','plant','n.mined.leaves','n.fruits','leaf_id','plant_id'),
                      variable.name = "rep",
                      value.name = 'ln_bac_host')
write.csv(nd3, file = here("models/NP/patch_level_Pseudomonas_pplr.csv"), quote = F)
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
prev_ests <- read.csv(here("models/NP/ppcfu_data/prev-ests_Pseudomonadaceae.csv"))
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
write.csv(nd5, file = here("models/NP/patch_level_Pseudomonas_ppcfu.csv"), quote = F)

#### Summaries and plots ####

# read in data, if not already loaded:
#nd5 <- read.csv(file = here("models/NP/patch_level_Pseudomonas_ppcfu.csv"))

# first, sum up per subplot per rep:
nd5_sp_rep_plant <- dplyr::group_by(nd5, subplot, plant, rep) %>% summarise(sp_sum_pp_cfu = sum(pp_cfu),
                                                                            log_sp_sum_pp_cfu = log(sp_sum_pp_cfu, 10),
                                                                            n.leaves.tot = length(sample_id),
                                                                            herb_dmg.sum = sum(herb_dmg),
                                                                            n.fruits = unique(n.fruits))

# because the number of stems per subplot is 16 and is invariant across the dataset, we do not need to average anything.
# instead, we simply sum up the pp. cfu (exponentiated) as well as the n.mined.leaves to get our patch-level intensities.
nd5_sp_rep <- dplyr::group_by(nd5_sp_rep_plant, subplot, rep) %>% summarise(sp_sum_pp_cfu.sp = sum(sp_sum_pp_cfu),
                                                                            log_sp_sum_pp_cfu.sp = log(sp_sum_pp_cfu.sp, 10),
                                                                            n.leaves.sp = length(n.leaves.tot),
                                                                            herb_dmg.sp = sum(herb_dmg.sum),
                                                                            n.fruits.sp = sum(n.fruits),
                                                                            n.stems = length(n.leaves.tot)) %>% arrange(herb_dmg.sp)

# order subplots by total herbivore load:
nd5_sp_rep$subplot <- factor(nd5_sp_rep$subplot, levels = paste0(unique(nd5_sp_rep$subplot)))


# now capture quantiles of this quantity among the reps:
nd5_sp_summary <- dplyr::group_by(nd5_sp_rep, subplot) %>% summarise(log_med_0.025 = quantile(log_sp_sum_pp_cfu.sp, 0.025),
                                                                     log_med_0.975 = quantile(log_sp_sum_pp_cfu.sp, 0.975),
                                                                     log_med_0.50  = quantile(log_sp_sum_pp_cfu.sp, 0.5),
                                                                     log_med_0.75  = quantile(log_sp_sum_pp_cfu.sp, 0.75),
                                                                     log_med_0.25  = quantile(log_sp_sum_pp_cfu.sp, 0.25),
                                                                     n.leaves.sp = unique(n.leaves.sp),
                                                                     herb_dmg.sp = unique(herb_dmg.sp),
                                                                     n.fruits.sp = unique(n.fruits.sp),
                                                                     n.stems = unique(n.stems)) %>% arrange(herb_dmg.sp)
#
#
# # merge back with other info from subplots:
# sp_dat <- unique(nd5[,c('plot','subplot','n.mined.leaves','n.fruits')]) %>% group_by(subplot) %>%
#   summarise(sum.n.mined.leaves = sum(n.mined.leaves),
#             sum.n.fruits = sum(n.fruits))
#
# nd5_sp_summary_all <- merge(nd5_sp_summary, sp_dat, by = 'subplot', sort = F) %>% arrange(desc(sum.n.mined.leaves))
#
# #### plots ####
# nd5_sp_summary_all$subplot <- factor(nd5_sp_summary_all$subplot, levels = paste0(nd5_sp_summary_all$subplot))

co_infection_p1 <- ggplot(nd5_sp_summary, aes(x = herb_dmg.sp, y = log_med_0.50)) +
  geom_linerange(aes(ymin = log_med_0.025, ymax = log_med_0.975), size = 0.33, col = "gray80") +
  geom_linerange(aes(ymin = log_med_0.25, ymax = log_med_0.75), size = 0.75, col = "gray70") +
  geom_point(alpha = 0.5, col = "gray20") +
  xlab("mined leaves per patch") +
  ylab("pp. bacterial load per patch") +
  scale_y_continuous(limits = c(6.5,10.1), breaks = seq(7,10,1)) +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150,25))
#co_infection_p1

# distribution of herbivory per sub-plot
dp1 <- ggplot(nd5_sp_summary, aes(x = herb_dmg.sp)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  scale_x_continuous(limits = c(0,150), breaks = seq(0,150,25))

# distribution of fruit set per sub-plot
dp2 <- ggplot(nd5_sp_summary, aes(x = n.fruits.sp)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  #scale_x_continuous(limits = c(0,250), breaks = seq(0,250,50)) +
  # theme(axis.line = element_blank(),
  #       axis.text = element_blank()) +
  coord_flip()

# distribution of median ppcfu per patch:
dp3 <- ggplot(nd5_sp_summary, aes(x = log_med_0.50)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  scale_x_continuous(limits = c(6.5,10.1), breaks = seq(7,10,1)) +
  coord_flip()

# output top graphs:
# ggarrange(plotlist = list(co_infection_p1,dp3), widths = c(1,0.6), ncol = 2, align = 'hv') %>%
#   ggsave(filename = here("figs/Fig3a.pdf"), width = 3.5, height = 2)


ggarrange(plotlist = list(dp1, dp1, #placeholder
                          co_infection_p1, dp3),
          widths = c(1,0.5),
          heights = c(0.5, 1),
          align = 'hv',
          ncol = 2, nrow = 2) %>% ggsave(filename = here("figs/Fig3_NP_v2.pdf"), width = 3, height = 2.5)


# fruit-set versus n.mined.leaves
# co_infection_p2 <- ggplot(nd5_sp_summary, aes(x = herb_dmg.sp, y = n.fruits.sp)) +
#   #geom_linerange(aes(ymin = log_med_0.025, ymax = log_med_0.975), size = 0.33, col = "gray80") +
#   #geom_linerange(aes(ymin = log_med_0.25, ymax = log_med_0.75), size = 0.75, col = "gray70") +
#   geom_point(alpha = 0.5, col = "gray20") +
#   xlab("mined leaves per patch") +
#   ylab("fruits set per patch")
  #scale_x_continuous(limits = c(0,140), breaks = seq(0,125,25)) +
  #scale_y_continuous(limits = c(0,250), breaks = seq(0,250,50))
#co_infection_p2

# ggarrange(plotlist = list(co_infection_p1,co_infection_p2), nrow = 2, align = 'hv') %>%
#   ggsave(filename = here("figs/Fig3ab.pdf"), width = 2.5, height = 4)
#
# ## save plot of density plot of n.mined.leaves
# ggarrange(plotlist = list(dp1,co_infection_p2), heights = c(0.6,1), nrow = 2, align = 'hv') %>%
#   ggsave(filename = here("figs/Fig3b.pdf"), width = 2, height = 3)

## bookmark ##
# let's figure out the plots.. then we can sort out whether the specific results are correct

# ggarrange(plotlist = list(dp1, dp1, #placeholder
#                           co_infection_p1, dp3,
#                           co_infection_p2, dp2),
#           widths = c(1,0.5),
#           heights = c(0.5, 1, 1),
#           align = 'hv',
#           ncol = 2, nrow = 3) %>% ggsave(filename = here("figs/Fig3_all.pdf"), width = 3, height = 4)

# try with seed_plot1
seed_hist1 <- readRDS(file = here("figs/seed_hist1.rds"))
seed_plot1 <- readRDS(file = here("figs/seed_plot1.rds"))

# ggarrange(plotlist = list(dp1, dp1, #placeholder
#                           co_infection_p1, dp3,
#                           seed_plot1, seed_hist1),
#           widths = c(1,0.5),
#           heights = c(0.5, 1, 1),
#           align = 'hv',
#           ncol = 2, nrow = 3) %>% ggsave(filename = here("figs/Fig3_all_NP_v2.pdf"), width = 3, height = 4)

ggarrange(plotlist = list(seed_plot1, seed_hist1),
          widths = c(1,0.5),
          heights = c(0.5, 1, 1),
          align = 'hv',
          ncol = 2, nrow = 1) %>% ggsave(filename = here("figs/Fig3_NP_seeds_v2.pdf"), width = 2.75, height = 1.75)



# now make plot of quantile of herbivore damage versus percent contribution to bacterial propagules:
herb_rank <- sort(nd5_sp_summary$herb_dmg.sp)
nd5_sp_summary$herbivory_rank <- sapply(nd5_sp_summary$herb_dmg.sp, function(x) { (length(herb_rank[herb_rank > x])/length(herb_rank)) })


# now calculate proportion of total propagules (median):
sum_propagules <- sum(10^nd5_sp_summary$log_med_0.50)
#nd5_sp_summary_all$propagule_fraction <- sapply(nd5_sp_summary_all$log_med_0.50, function(x) x^10/sum_propagules)

# generate cumulative sum of propagules as one increases herbivory_rank
cumsum_herb <- nd5_sp_summary[,c('herbivory_rank','log_med_0.50')]
cumsum_herb <- dplyr::arrange(cumsum_herb, herbivory_rank)
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

ggsave(scale_plot1, filename = here("figs/NP_scale_plot1_v2.pdf"), width = 2.5, height = 2)

