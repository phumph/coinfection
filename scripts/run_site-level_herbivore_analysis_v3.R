#! /usr/bin/Rscript
# Rscript to generate site-level analysis of herbivore patterns
# last updated 2018-SEP-18 by PTH: initial draft of analysis
# last updated 2018-OCT-25 by PTH: modified to generate summary plots
# last updated 2018-OCT-30 by PTH: modified to run models for EL and NP
# last updated 2019-MAR-26 by PTH: added additional model comparisons using kfold and loo

library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# load dataset to supply the newdata for predictive simulations
#FOCAL.DAT.NAME <- here("data/EL.field.hormone.exp.master.txt")
FOCAL.DAT.NAME <- here("data/NP_field_herbivory_master_2018-OCT-25.txt")
SITE <- ifelse(length(grep('EL',FOCAL.DAT.NAME,value=T))>0,'EL','NP')
focal.dat <- read.table(FOCAL.DAT.NAME,T,"\t")
if (SITE == 'NP'){
  focal.dat$n.fruits <- as.numeric(as.vector(focal.dat$n.fruits))
}
focal.dat$subplot <- with(focal.dat, paste0(plot,'_',tx))
#focal.dat$plant <- with(focal.dat, paste0(sub.plot, stem.num))
focal.dat$condition <- relevel(focal.dat$condition, ref = 'CTR')

# plot of paired herbivore damage between sub-plots for each treatment:
focal.dat2 <- dplyr::group_by(focal.dat, condition, plot, tx) %>% summarise(tot.mined.leaves = sum(n.mined.leaves, na.rm = T),
                                                                            sd.mined.leaves = sd(n.mined.leaves, na.rm = T),
                                                                            tot.n.leaves = sum(n.leaves, na.rm = T),
                                                                            n.stems = length(stem.num),
                                                                            mean.n.leaves = tot.n.leaves / n.stems,
                                                                            se.mean.n.leaves = sd(n.leaves) / sqrt(n.stems),
                                                                            mean.mined.leaves = tot.mined.leaves/n.stems,
                                                                            se.mined.leaves = sd.mined.leaves/sqrt(n.stems),
                                                                            mean.height = mean(height, na.rm=T),
                                                                            se.mean.height = sd(height) / sqrt(n.stems),
                                                                            mean.prop.mined.leaves = mean(n.mined.leaves/n.leaves),
                                                                            region = unique(region))

# dcast for plotting:
focal.dat3a <- reshape2::dcast(focal.dat2, condition + plot + region ~ tx, value.var = 'mean.prop.mined.leaves')
#focal.dat3b <- reshape2::dcast(focal.dat2, condition + plot ~ tx, value.var = 'mean.mined.leaves')
focal.dat3b <- reshape2::dcast(focal.dat2, condition + plot + region ~ tx, value.var = 'tot.mined.leaves')
focal.dat3b2 <- reshape2::dcast(focal.dat2, condition + plot ~ tx, value.var = 'se.mined.leaves')
focal.dat3c <- reshape2::dcast(focal.dat2, condition + plot ~ tx, value.var = 'mean.n.leaves')
focal.dat3d <- reshape2::dcast(focal.dat2, condition + plot ~ tx, value.var = 'mean.height')
#focal.dat3e <- reshape2::dcast(focal.dat2, condition + plot ~ tx, value.var = 'sum.n.fruits')


# make plot with error bars:
names(focal.dat3b2)[c(3,4)] <- c('MOCK.se','TX.se')
#names(focal.dat3e)[c(3,4)] <- c('MOCK.seeds','TX.seeds')
#mined.dat1 <- cbind(focal.dat3b,focal.dat3b2[,c(3:4)],focal.dat3e[,c(3:4)])
mined.dat1 <- cbind(focal.dat3b,focal.dat3b2[,c(3:4)])

## plot paired data across conditions:

# mean proportion mined leaves
fd_pa <- ggplot(focal.dat3a, aes(x = MOCK, y = TX)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~ condition, ncol = 3)
  #facet_grid(region ~ condition)

# total mined leaves
fd_pb <- ggplot(focal.dat3b, aes(x = MOCK, y = TX)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~ condition, ncol = 3)

# total mined leaves broken down by region
fd_pb_reg <- ggplot(focal.dat3a, aes(x = MOCK, y = TX)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + facet_grid(region ~ condition)

# total mined leaves - histogram
fd_pb2 <- ggplot(focal.dat2, aes(x = tot.mined.leaves)) +
  geom_histogram(alpha = 0.4) +
  facet_wrap( ~ condition, nrow=3)# + coord_flip()

# mean n leaves - just for checking
fd_pc <- ggplot(focal.dat3c, aes(x = MOCK, y = TX)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~ condition, ncol = 3)

# mean stem height - just for checking
fd_pd <- ggplot(focal.dat3d, aes(x = MOCK, y = TX)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~ condition, ncol = 3)

#ggarrange(plotlist = list(fd_pa,fd_pb,fd_pc,fd_pd), nrow = 4)

# this doesn't work yet
# # do separate plot with error bars:
# fd_pb2 <- ggplot(mined.dat1, aes(x = MOCK, xmin = MOCK-MOCK.se, xmax = MOCK+MOCK.se, y = TX, ymin = TX-TX.se, ymax = TX+TX.se)) + geom_point() +
#   geom_abline(intercept = 0, slope = 1) + facet_wrap(~ condition, ncol = 3) +
#   geom_errorbar() +
#   geom_errorbarh() +
#   scale_x_continuous(limits = c(-0.5,9.5), breaks = seq(0,8,2)) +
#   scale_y_continuous(limits = c(-0.5,9.5), breaks = seq(0,8,2)) +
#   ggtitle(SITE)
#
# ggsave(fd_pb2, filename = paste0(here("figs/"),"field_herbivory_",SITE,"_by_TX.pdf"), width = 5, height = 2)


#### statistical models of herbivory distributions ####

# first, mean center and standardize predictors:
focal.dat2$mean.height.z <- (focal.dat2$mean.height - mean(focal.dat$height)) / sd(focal.dat$height)
focal.dat2$n.leaves.z <- (focal.dat2$mean.n.leaves - mean(focal.dat$n.leaves)) / sd(focal.dat$n.leaves)

# transform standard errors by dividing by standard deviation:
focal.dat2$se.mined.leaves.z <- focal.dat2$se.mined.leaves / sd(focal.dat$n.leaves)
focal.dat2$se.mean.height.z <- focal.dat2$se.mean.height / sd(focal.dat$height)
focal.dat2$se.n.leaves.z <- focal.dat2$se.mean.n.leaves / sd(focal.dat$n.leaves)


#### NB MODEL ####
# try NB model with flexible shape parameter as well as location parameters:

fd_pb2log <- ggplot(focal.dat2, aes(x = log(1+tot.mined.leaves))) +
  geom_histogram(alpha = 0.4) +
  facet_wrap( ~ condition, nrow=3)# + coord_flip()

# need to compare model predicted distribution to these histograms to inspect fit.
nb_3 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + condition + (condition|plot)),
            data = focal.dat2,
            family = negbinomial,
            cores = 4,
            iter = 10000,
            control = list(adapt_delta = 0.99))

# nb_2 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + condition + (1|plot), shape = ~ 1),
#             data = focal.dat2,
#             family = negbinomial,
#             iter = 5000,
#             control = list(adapt_delta = 0.95),
#             cores = 4)

# trying right-censored model since the response variable cannot exceed the number of total leaves per patch!
focal.dat2$censored <- 'right'
nb_1 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + condition + (1|plot)),
                 data = focal.dat2,
                 family = negbinomial,
                 iter = 1000,
                 control = list(adapt_delta = 0.95),
                 cores = 4)

# nb_1_cens <- brm(bf(tot.mined.leaves | cens(censored, tot.n.leaves) ~ offset(log(tot.n.leaves)) + condition + (1|plot)),
#                  data = focal.dat2,
#                  family = negbinomial,
#                  iter = 1000,
#                  control = list(adapt_delta = 0.95),
#                  cores = 4)

nb_0 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + (1|plot)),
            data = focal.dat2,
            family = negbinomial,
            iter = 2000,
            control = list(adapt_delta = 0.96),
            cores = 4)

p_0 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + (1|plot)),
            data = focal.dat2,
            family = poisson,
            iter = 2000,
            control = list(adapt_delta = 0.96),
            cores = 4)

p_1 <- brm(bf(tot.mined.leaves ~ offset(log(tot.n.leaves)) + condition + (1|plot)),
           data = focal.dat2,
           family = poisson,
           iter = 2000,
           control = list(adapt_delta = 0.96),
           cores = 4)


kfold(nb_0,nb_1,nb_3,p_0,p_1,K=10)

# save models
saveRDS(nb_0,file = here('models/herbivory/nb_0.rds'))
saveRDS(nb_1,file = here('models/herbivory/nb_1.rds'))
saveRDS(nb_2,file = here('models/herbivory/nb_2.rds'))
saveRDS(nb_3,file = here('models/herbivory/nb_3.rds'))



# load up models for inspection:
nb_0 <- readRDS(file = here('models/herbivory/nb_0.rds'))
nb_1 <- readRDS(file = here('models/herbivory/nb_1.rds'))
nb_2 <- readRDS(file = here('models/herbivory/nb_2.rds'))
nb_3 <- readRDS(file = here('models/herbivory/nb_3.rds'))


# compare models via LOO-IC


#### generate output of model parameters ####
# define focal model
focal_mod <- nb_1

# grab coefficients
fm1.c <- tidy(focal_mod, prob = 0.95)
fm2 <- cbind(term = fm1.c[,1], round(fm1.c[,-c(1)], 2)) # round and create output data.frame
fm2$coef <- with(fm2, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string

fm3 <- fm2[,c('term','coef')]
fm3 <- rbind(fm3[grep('^b_',fm3$term),],
             fm3[grep('^sd_',fm3$term),],
             fm3[grep('^shape',fm3$term),])

# output:
fm3$term <- c('Intercept (CTR)','beta JA','beta SA','sigma(plot)','NB shape')

# write lines to output:
con <- file(here("tables/FieldHerbivory_model_table_v3.tex"))
open(con, 'wr')
writeLines(kable(fm3, "latex", caption = "Model estimates for impact of hormone treatment on patch-level herbivore abundance in the field (site NP)", booktabs = T, escape = TRUE) %>%
             kable_styling(position = "center"),
           con = con, sep = ''
)
close(con)


#### posterior predictions ####

# spit out pp check broken down by condition:
#yrep <- posterior_predict(focal_mod, nsamples = 1000)
# make new factor levels:

# let's look at marginal effects:
me1 <- marginal_effects(focal_mod, probs = c(0.025,0.975))
me2 <- marginal_effects(focal_mod, probs = c(0.25,0.75))

# re-level:
me1[[1]]$condition <- factor(me1[[1]]$condition, levels = c('SA','JA','CTR'))
me2[[1]]$condition <- factor(me1[[1]]$condition, levels = c('SA','JA','CTR'))

# let's try to build simple marginal effects plot in ggplot with these data for export:
me_p1 <- ggplot() +
  geom_linerange(data = me1[[1]], aes(x=condition,ymin=lower__,ymax=upper__), lwd = 0.75, col = 'gray60') +
  geom_linerange(data = me2[[1]], aes(x=condition,ymin=lower__,ymax=upper__), lwd = 2, col = 'gray60') +
  geom_point(data = me1[[1]],aes(x=condition,y=estimate__), col = "black", size = 5,pch="|") +
  #geom_crossbar(data = me1[[1]], aes(x=condition,y=estimate__, ymin=estimate__,ymax=estimate__), col = "black", size = 0.25) +
  coord_flip() +
  scale_y_continuous(limits = c(0,150), breaks = seq(0,150,50))

# export plot:
ggsave(filename = here('figs/NP_herbivory_arginal_effects_p1.pdf'), width = 3, height = 1.5)


#### potentially obsolete below this line ####
the_plots <- as.factor(as.numeric(focal.dat2$plot)+100)
newdat1 <- data.frame(condition = focal.dat2$condition,
                      plot = the_plots,
                      tot.n.leaves = focal.dat2$tot.n.leaves)

yrep <- posterior_predict(focal_mod,
                          allow_new_levels = TRUE,
                          newdata = newdat1, nsamples = 1000) # generate new levels for plot ID:

yrep_dat <- as.data.frame(t(yrep))
newdat2 <- data.frame(newdat1,
                      yrep_dat)
newdat2 <- reshape2::melt(newdat2, value.var = 'count', variable.name = 'rep', id.vars = c('condition','plot','tot.n.leaves'))
newdat2$rep <- sapply(newdat2$rep, function(x) gsub('V','',paste0(x)))
newdat2$rep <- as.numeric(newdat2$rep)
newdat2$rep <- factor(newdat2$rep)

#newdat2 <- dplyr::filter(newdat2, rep < 20)


# try also to generate new random effects levels to see whether the effects will generalize.
# this model basically reflects the fact that fluctuations in plot-level intercepts may be responsible for the marginal effect of JA
# but that this won't generalize across new plots, based on the posterior distributions of model parameters.

# now generate plot:
# pp_plot <- ggplot(newdat2) +
#   geom_density(aes(x = value, fill = rep), alpha = 0.01, color = 'gray80', weight = 0.1) +
#   #geom_line(aes(x = value, group = 'rep'), stat="density", alpha=0.4) +
#   facet_wrap( ~ condition) + theme(legend.position='none')


# calculate posteior difference between JA and CTR plots per rep
# produce posterior predicted difference



# calculate posterior predicted means and over-plot observed:
newdat3 <- dplyr::group_by(newdat2, condition, rep) %>% summarise(means = mean(value))

# now summarise mean distribution from yrep data:
newdat4 <- dplyr::group_by(newdat3, condition) %>% summarise(q0.025 = quantile(means,0.025),
                                                             q0.25 = quantile(means,0.25),
                                                             q0.50 = quantile(means,0.50),
                                                             q0.75 = quantile(means,0.75),
                                                             q0.975 = quantile(means,0.975),
                                                             mu_mu = mean(means))

focal.dat.sums <- dplyr::group_by(focal.dat2, condition) %>% summarise(means = mean(tot.mined.leaves))
newdat4$condition <- factor(newdat4$condition, levels = rev(c('CTR','JA','SA')))
pp_means <- ggplot() +
  #geom_jitter(data = newdat3, aes(y = means, x = condition), width = 0.1, alpha = 0.2, col = 'gray40') +
  geom_linerange(data = newdat4, aes(x = condition, ymin = q0.025, ymax = q0.975),lwd = 0.75, col = 'gray60') +
  geom_linerange(data = newdat4, aes(x = condition, ymin = q0.25, ymax = q0.75), lwd = 2, col = 'gray60') +
  scale_y_continuous(limits = c(0, 205), breaks = seq(0,200,50)) +
  geom_crossbar(data = newdat4, aes(x = condition, y = q0.50, ymin = q0.50, ymax = q0.50), col = "black", lwd = 0.25) +
  xlab("") + ylab("") +
  #geom_line(data = newdat4, aes(x = condition, y = q0.50, group = condition), col = 'orange') +
  #geom_point(data = focal.dat.sums, aes(x = condition, y = means), col = "orange", size = 3) +
  #geom_point(data = focal.dat.sums, aes(x = condition, y = means), col = "white", size = 1) +
  coord_flip()

# make histogram figure to pair with it:
obs_hist <- ggplot(data = focal.dat2, aes(x = tot.mined.leaves)) +
  geom_histogram(bins = 30, fill = 'gray60', col = "white") +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0,150,50)) +
  facet_wrap(~ condition, ncol = 1)

## export plot:
ggarrange(plotlist = list(pp_means, obs_hist), nrow = 2, heights = c(0.33,1)) %>%
  ggsave(filename = here("figs/herb_hist_and_ppmeans.pdf"), width = 3, height = 5)


# model structure: mean.mined.leaves (with se estimate) ~ condition + condition:tx + (condition|plot)
pois_m1 <- brm(bf(tot.mined.leaves ~ n.leaves.z + mean.height.z + condition + condition:tx + (1|plot)),
               data = focal.dat2,
               family = poisson,
               cores = 4,
               iter = 10000)

pois_m1b <- brm(bf(tot.mined.leaves ~ me(n.leaves.z,se.n.leaves.z) + me(mean.height.z,se.mean.height.z) + condition + (1|plot)),
               data = focal.dat2,
               family = poisson,
               cores = 4,
               iter = 10000)

pois_m1c <- brm(bf(tot.mined.leaves ~ me(n.leaves.z,se.n.leaves.z) + me(mean.height.z,se.mean.height.z) + condition + (1|plot)),
                data = focal.dat2,
                family = zero_inflated_poisson,
                cores = 4,
                iter = 10000)

summary(pois_m1c)

# see if there is graphical support for bimodal (log) normal distribution of damage:
mean.mined.hist <- ggplot(focal.dat2, aes(x = mean.mined.leaves)) + geom_histogram(bins = 30) +
  facet_wrap(~ condition)# yep looks like it to me!

tot.mined.hist <- ggplot(focal.dat2, aes(x = tot.mined.leaves)) + geom_histogram(bins = 30) +
  facet_wrap(~ condition)# yep looks like it to me!



# try to plot histograms of stem-level sums;
ggplot(focal.dat, aes(x = n.mined.leaves)) + geom_histogram(bins = 20) +
  facet_wrap(~ condition)# yep looks like it to me!

# let's try poisson model with zero inflation:
zip1 <- brm(bf(n.mined.leaves ~ condition + (1|plot) + (1|subplot), zi =~ condition),
            data = focal.dat,
            family = zero_inflated_poisson,
            cores = 4,
            chains = 4)

zip2 <- brm(bf(n.mined.leaves ~ condition + (1|plot) + (1|subplot)),
            data = focal.dat,
            family = zero_inflated_poisson,
            cores = 4,
            chains = 4)

zip3 <- brm(bf(n.mined.leaves ~ condition + (condition|plot) + (1|subplot)),
            data = focal.dat,
            family = zero_inflated_poisson,
            cores = 4,
            chains = 4)

LOO(zip1,zip2,zip3)

# OK no real differences here.. but what about patch-level abundance variation? This isn't exactly helpful.

# OK there appears to be a good mixture; setup model in brms to estimate mixing proportions
the_mix0 <- mixture(gaussian,gaussian)
the_mix <- mixture(poisson,poisson)
the_mix_gamma <- mixture(gamma(),gamma())
the_mix2 <- mixture(negbinomial,negbinomial)
the_mix3 <- mixture(zero_inflated_poisson,poisson)

# try additional mixture models:
mix0 <- mixture(negbinomial, negbinomial)
mix0 <- mixture(poisson, poisson)
#
# prior <- c(
#   #prior(normal(0, 5), Intercept, nlpar = mu1),
#   #prior(normal(0, 5), Intercept, nlpar = mu2),
#   prior(dirichlet(2, 2), theta2)
# )

# need to constrain parameter space a bit more with stronger priors to keep estimation in bounds.
# can constrain thetas to be normal(0,4) on logit scale
# as well as the various means.. I know that the

# trying to identify model without plot-level random effect:
mixPOIS1 <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition),
               data = focal.dat2,
               family = mix0,
               #prior = prior,
               inits = 0,
               iter = 20000,
               chains = 4,
               cores = 4,
               control = list(adapt_delta = 0.98,
                              max_treedepth = 15))

hist(rpois(100,75))
#focal.dat2$log.tot.mined.leaves <- log(1+focal.dat2$tot.mined.leaves,10)

# prior <- c(
#   prior(normal(1, 2), Intercept, dpar = mu1),
#   prior(normal(1.5, 2), Intercept, dpar = mu2)
# )

mixGA <- brm(bf(log.tot.mined.leaves ~ 1, theta2 ~ condition),
               data = focal.dat2,
               family = the_mix0,
               #prior = prior,
               inits = 0,
               iter = 8000,
               chains = 4,
               cores = 4,
               control = list(adapt_delta = 0.999))

mixPP4 <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition),
              data = focal.dat2,
              family = the_mix3,
              #prior = prior,
              inits = 0,
              iter = 8000,
              chains = 4,
              cores = 4,
              control = list(adapt_delta = 0.999))

mixPP4b <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition + (1|plot)),
              data = focal.dat2,
              family = the_mix3,
              #prior = prior,
              inits = 0,
              iter = 8000,
              chains = 4,
              cores = 4,
              control = list(adapt_delta = 0.999))

# try to examine the posterior predictions more generally:
pp_check(mixPP4b,
         yrep = posterior_predict(mixPP4b, nsamples = 4000), type = "freqpoly_grouped", group = 'condition')

# need to try again the Gaussian mixture model, with different sigmas depending on



mixPP2 <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition + (1|plot)),
              data = focal.dat2,
              family = the_mix,
              #prior = prior,
              inits = 0,
              iter = 8000,
              chains = 4,
              cores = 4,
              control = list(adapt_delta = 0.999))

mixPP3 <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition + (condition|plot)),
              data = focal.dat2,
              family = the_mix,
              #prior = prior,
              inits = 0,
              iter = 8000,
              chains = 4,
              cores = 4,
              control = list(adapt_delta = 0.999))

## I'm tempted to include a zero inflation parameter in the poisson mixture model.


mixNBNB <- brm(bf(tot.mined.leaves ~ 1, theta2 ~ condition),
               data = focal.dat2,
               family = the_mix2,
               #prior = prior,
               inits = 0,
               iter = 8000,
               chains = 4,
               cores = 4,
               control = list(adapt_delta = 0.999))

summary(mix3)
marginal_effects(mix3)

pp_check(mix4)

summary(mix4)
marginal_effects(mix4)

# simply trying to run a negative binomial model on these counts:
nb1 <- brm(bf(tot.mined.leaves ~ tot.n.leaves + mean.height + condition + (1|plot)),
            data = focal.dat2,
            family = negbinomial,
            #prior = prior,
            inits = 0,
            iter = 8000,
            chains = 4,
            cores = 4,
            control = list(adapt_delta = 0.99))

# update and add covariates:
nb1b <- update(nb1, formula. = ~ . - tot.n.leaves, newdata = focal.dat2)
nb1c <- update(nb1b, formula. = ~ . - mean.height, newdata = focal.dat2)
nb1d <- update(nb1b, formula. = ~ . - mean.height - tot.n.leaves, newdata = focal.dat2)

# compare model fits:
LOO(nb1,nb1b,nb1c,nb1d)

# trying non-parametric test:
focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']),mean = 0, sd = 0.01)
wilcox.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']),mean = 0, sd = 0.01))

wilcox.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']),mean = 0, sd = 0.01))

wilcox.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']),mean = 0, sd = 0.01))

ks.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']),mean = 0, sd = 0.01))

ks.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='JA']),mean = 0, sd = 0.01))

ks.test(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='CTR']),mean = 0, sd = 0.01),
            focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']+rnorm(length(focal.dat2$tot.mined.leaves[focal.dat2$condition=='SA']),mean = 0, sd = 0.01))


nb2 <- brm(bf(tot.mined.leaves ~ condition + (condition|plot)),
            data = focal.dat2,
            family = negbinomial,
            #prior = prior,
            inits = 0,
            iter = 8000,
            chains = 4,
            cores = 4,
            control = list(adapt_delta = 0.99))


nb3 <- brm(bf(tot.mined.leaves ~ condition + (1|plot), shape ~ condition),
            data = focal.dat2,
            family = negbinomial,
            #prior = prior,
            inits = 0,
            iter = 8000,
            chains = 4,
            cores = 4,
            control = list(adapt_delta = 0.999))

pp_check(mix5, group = "condition")
summary(mix5)
LOO(mix4,mix5,mix6)
marginal_effects(mix5)

# somehow add zero inflation here.. gah!

#### previous analysis below this line ####


# do some exploratory plotting of aggregated herbivory levels across sub-plots:
subplot_data <- dplyr::group_by(focal.dat, subplot, condition, tx, plot) %>% summarise(n.mined.leaves.sum = sum(n.mined.leaves),
                                                                                       n.leaves.sum = sum(n.leaves),
                                                                                       height.avg = mean(height, na.rm = T))
subplot_data$group = paste0(subplot_data$subplot,'_',subplot_data$condition)
subplot_data$plot <- factor(subplot_data$plot)

ggplot(subplot_data, aes(x = tx, y = n.mined.leaves.sum, group = plot)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ condition)

# run model to estimate within-treatment effects:
el.ja.m1 <- brm(bf(n.mined.leaves ~ n.leaves + height + tx + (tx|plot) + (1|subplot)),
                  data = focal.dat[focal.dat$condition=='JA',],
                  family = negbinomial,
                  cores = 4)

# try model where counts are summed. Should get closer to normal approximation.
# this model seems over-parameterized.
el.ja.m2 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot) + (1|subplot)),
                data = subplot_data[subplot_data$condition=='JA',],
                family = gaussian,
                cores = 4)

el.ja.m3 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot)),
                data = subplot_data[subplot_data$condition=='JA',],
                family = gaussian,
                cores = 4,
                control = list(adapt_delta = 0.99))
plot(el.ja.m3)
el.ja.m4 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot)),
                data = subplot_data[subplot_data$condition=='CTR',],
                family = gaussian,
                cores = 4,
                control = list(adapt_delta = 0.99))

# ggplot(subplot_data, aes(x = n.leaves.sum, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = n.mined.leaves.sum, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = n.leaves.sum, y = n.mined.leaves.sum, col = condition)) + geom_point()
#
#
# ggplot(subplot_data, aes(x = condition, y = n.leaves.sum, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = condition, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = condition, y = n.mined.leaves.sum, col = condition)) + geom_point()


