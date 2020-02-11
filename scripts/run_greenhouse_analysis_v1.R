#! /usr/bin/Rscript
# run_greenhouse_analysis_v1.R
# script to run greenhouse analyses
# last updated: 2019-AUG-07 by PTH

# header
library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

library(parallel)
options(mc.cores = parallel::detectCores())

# step 1: find data and get into analyzable format
gh1 <- read.csv(file = here("data/greenhouse_doublings_data.csv"))

# run model:
all_strains_bm1 <- brm(bf(n.doublings ~ 0 + strain.id + plant.tx:strain.id + (1|plant.id)),
                       data = gh1,
                       family = gaussian(),
                       iter = 8000,
                       cores = getOption("mc.cores", 1L))

# add lOO to these models; save:
all_strains_bm1 <- add_ic(all_strains_bm1, ic = 'loo', reloo=T)
saveRDS(all_strains_bm1, file = here("models/gh_bm1_v3.rds"))

# load for plotting:
if(!exists("all_strains_bm1")) {
  all_strains_bm1 <- readRDS(file = here("models/gh_bm1_v3.rds"))
}

# now collect and export coefficients:
library(broom)
bm1_coefs <- tidy(all_strains_bm1)
focal_coefs <- bm1_coefs[grep(':plant.tx', bm1_coefs$term),]
post_bm1 <- data.frame(all_strains_bm1)

focal_coefs <- post_bm1[ , names(post_bm1) %in% grep('.plant.tx', names(post_bm1), value=T)]
coef_res <- data.frame(t(apply(focal_coefs, 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))
coef_res$strain.id <- sapply(row.names(coef_res), function(x) gsub('b_strain.id', '', x))
coef_res$strain.id <- sapply(coef_res$strain.id, function(x) gsub('.plant.tx', '', x))

# re-order strain.id to match phylogeny:
strain_order <- c("20A","22B","26B","46B","02A","20B","36A","29A","39A","33E","03A","46A")
coef_res$strain.id <- factor(coef_res$strain.id, levels = rev(strain_order))
names(coef_res) <- c("q0.025","q0.25","q0.50","q0.75","q0.975","strain.id")

gh2 <- dplyr::group_by(gh1, strain.id, plant.tx) %>% summarise(log.mu.cfu.g = log(mean(10^log.avg.g),10),
                                                               log.med.cfu.g = log(median(10^log.avg.g),10))
gh2$strain.id <- factor(gh2$strain.id, levels = rev(strain_order))

#### PLOTS ####

gh_p1 <- ggplot(gh2, aes(y = strain.id, x = factor(plant.tx))) +
  geom_tile(aes(fill = log.med.cfu.g)) +
  #scale_fill_gradient2(low = "steelblue", high = "darkorange2", midpoint = median(fam_sums$log_med_0.50)) +
  #scale_fill_gradientn(colors = cividis(20)) +
  #scale_fill_gradientn(colors = magma(20)) +
  #scale_fill_gradientn(colors = plasma(20)) +
  scale_fill_gradientn(colors = parula(64)) +
  #scale_fill_gradient(low = "white", high = "black") +
  theme(legend.position = 'top') + xlab("")

# now for the posterior of n.doublings different in JA versus MOCK:
doubling_diffs <- ggplot(coef_res) +
  geom_linerange(aes(x = strain.id,
                     #y = log2_med_med,
                     ymax = q0.975,
                     ymin = q0.025),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 0.5) +
  geom_linerange(aes(x = strain.id,
                     #y = log2_med_med,
                     ymax = q0.25,
                     ymin = q0.75),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 1.25) +
  geom_point(aes(x = strain.id, y = q0.50), pch = '|', size = 3, col = "white") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, col = "black", lty = 'dotted') +
  coord_flip() +
  theme(axis.text.y = element_blank()) + xlab("")

ggarrange(plotlist = list(gh_p1, doubling_diffs),
          ncol = 2,
          widths = c(0.4,1),
          common.legend = T, align = 'h') %>%
  ggsave(filename = here("figs/Fig4_parula_v2.pdf"), width = 3, height = 3)

#### posterior simulations of communities ####

# steps:
  # 1. sample from joint posterior
  # 2. concatenate abundances of predicted doublings per leaf class
  # 3. use baseline density as 10^5
  # 4. create community matrixes of medians across leaf types.. or just use marginal effects to simualate this.

# re-load up model:
if(!exists("all_strains_bm1")) {
  all_strains_bm1 <- readRDS(file = here("models/gh_bm1_v3.rds"))
}
strain_order <- c("20A","22B","26B","46B","02A","20B","36A","29A","39A","33E","03A","46A")

# version 1.0: ignoring group-level term and just drawing from the joint posterior of the intercept and slope for each strain (as well as from the residual distribution)
# version 2.0: this would add additional variance in the form of plant-level intercept effects. Could sample new arbitrary levels for this

# construct newdat
plant.tx <- c(1, 2)
newdat <- expand.grid(strain.id = strain_order, plant.tx = plant.tx)
newdat <- cbind(newdat,
                plant.id = c(rep('A',length(strain_order)),rep('B',length(strain_order)))
)

# add spp column:
newdat <- merge(newdat, data.frame(unique(gh1[,names(gh1) %in% c('strain.id','spp')])), sort = F)

# I want two leaves, each of one treatment, and to sample new levels of the random effects each simulation:
ppgh1 <- brms::posterior_predict(all_strains_bm1,
                                 newdata = newdat,
                                 allow_new_levels = TRUE,
                                 sample_new_levels = 'gaussian',
                                 nsamples = 200)

# transform into abundances:
ppgh2 <- round(2^ppgh1)

# add back to original data.frame and plot distributions
ppgh3 <- cbind(newdat, t(ppgh2))
ppgh3b <- reshape2::melt(ppgh3, id.vars = c('strain.id','plant.id','plant.tx','spp'), variable.name = 'rep', value.name = 'ppcfu')

gh_breakdown0 <- ggplot(ppgh3b, aes(x = strain.id, y = log(ppcfu,10), col = factor(plant.tx))) +
  geom_jitter(position = position_dodge(1), alpha = 0.2) +
  geom_boxplot(position = position_dodge(1), alpha = 0.2) +
  facet_wrap(~ spp, scales = 'free_x') +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") +
  scale_y_continuous(limits = c(2.25, 12.5), breaks = seq(2.5,12.5,2.5))

ppgh4a <- dplyr::group_by(ppgh3b, plant.tx, rep) %>% summarise(sum_ppcfu = sum(ppcfu))
ppgh4a$A <- 1
gh_breakdown1 <- ggplot(ppgh4a, aes(x = factor(plant.tx), y = log(sum_ppcfu,10), col = factor(plant.tx))) +
  geom_jitter(width = 0.1, alpha = 0.2) + geom_boxplot(alpha = 0.2) +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") + facet_wrap(~ A) +
  scale_y_continuous(limits = c(2.25, 12.5), breaks = seq(2.5,12.5,2.5))

ppgh4b <- dplyr::group_by(ppgh3b, plant.tx, spp, rep) %>% summarise(sum_ppcfu = sum(ppcfu))

# posterior predicted total bacteria in JA treated leaves:
gh_breakdown2 <- ggplot(ppgh4b, aes(x = factor(plant.tx), y = log(sum_ppcfu,10), col = factor(plant.tx))) +
  geom_jitter(width = 0.1, alpha = 0.2) + geom_boxplot(alpha = 0.2) +
  facet_wrap(~ spp) +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") +
  scale_y_continuous(limits = c(2.25, 12.5), breaks = seq(2.5,12.5,2.5))

ggarrange(plotlist = list(gh_breakdown1, gh_breakdown2, gh_breakdown0),common.legend = T, widths = c(0.265,0.455,1.125), ncol = 3) %>%
  ggsave(filename = here("figs/gh_structure_v1.pdf"), width = 6, height = 2.5)

# try to also distinguish among the community profiles when simulated...
# under assumption of no interactions among strains, JA treatment should increase total Pseudomonas count while
# also shifting the relative abundances of clades and strains within clades, favoring a minority of P. syringae strains.

# do this tomorrow! See if I can make some sense of these simulated outcomes at the community level.

# sum together abundances of P. syringae and P. fluorescens:


#### PP simulations of Genus, Clade, and Strain-level abundance patterns --> community structure ####

# use model coefficients to do a proper estimation.

# construct newdat
plant.tx <- c(1,2)
plant.id <- c(1:50)

# replicate this across population of n=50 plants, to constitute a reasonable sample over which to take medians.
# I don't want these predictions to feature sampling noise, only the uncertainty in the parameter estiamtes themselves.
newdat <- expand.grid(strain.id = strain_order, plant.tx = plant.tx, plant.id = plant.id)
# add spp column:
newdat <- merge(newdat, data.frame(unique(gh1[,names(gh1) %in% c('strain.id','spp')])), sort = F)
newdat$unique.id <- with(newdat,paste0(plant.id,'_',plant.tx))
newdat2 <- split(newdat, newdat$unique.id)
R <- 200 # define number of posterior simulations to generate per simulated plant sample

# generate posterior draws for each plant.id
gh_pp1 <- lapply(newdat2, function(x) brms::posterior_predict(all_strains_bm1,
                                                              newdata = x,
                                                              nsamples = 200,
                                                              allow_new_levels = TRUE,
                                                              sample_new_levels = 'gaussian'))

# add back to original newdat:
for (i in seq_along(newdat2)){
  newdat2[[i]] <- cbind(newdat2[[i]], t(2^gh_pp1[[i]]))
}

# bind all rows back together
gh_pp2 <- do.call(rbind, newdat2)

# melt for summary statistics and diversity calculations:
gh_pp3 <- reshape2::melt(gh_pp2, id.vars = c('strain.id','plant.tx','plant.id','spp','unique.id'),
                         variable.name = 'rep',
                         value.name = 'ppcfu')

gh_pp4 <- dplyr::group_by(gh_pp3, strain.id, plant.tx, spp, rep) %>% summarise(log_med_ppcfu = log(median(ppcfu),10))

# now sum up across strains within each clade for each rep:
gh_pp5 <- dplyr::group_by(gh_pp4, plant.tx, spp, rep) %>% summarise(log_sum_med_ppcfu = log(sum(10^log_med_ppcfu),10))

# now, naturally, we can sum up both clades and reveal the difference at this phylogenetic level:
gh_pp6 <- dplyr::group_by(gh_pp4, plant.tx, rep) %>% summarise(log_sum_med_ppcfu = log(sum(10^log_med_ppcfu),10))
# add mock factor for plot dimension consistency:
gh_pp6$Genus <- 'Pseudomonas'

# PLOTS
# now I have medians across the R replicates for each strain for each treatment. Let's plot:
gh_breakdown_strain <- ggplot(gh_pp4, aes(x = strain.id, y = log_med_ppcfu, col = factor(plant.tx))) +
  geom_jitter(position = position_dodge(1), alpha = 0.2, size = 0.5) +
  geom_boxplot(position = position_dodge(1), alpha = 0.5, outlier.size = 0) +
  facet_wrap(~ spp, scales = 'free_x') +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") +
  xlab("plant tx") +
  scale_y_continuous(limits = c(4, 10), breaks = seq(4,10,1))

# now I have the distribution of sums of the median abundances per strain, broken down by clade and plant.tx.
# plot:
gh_breakdown_clade <- ggplot(gh_pp5, aes(x = factor(plant.tx), y = log_sum_med_ppcfu, col = factor(plant.tx))) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.size = 0) +
  facet_wrap(~ spp, scales = 'free_x') +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") +
  xlab("plant tx") +
  scale_y_continuous(limits = c(6, 10), breaks = seq(6,10,1))

gh_breakdown_fam <- ggplot(gh_pp6, aes(x = factor(plant.tx), y = log_sum_med_ppcfu, col = factor(plant.tx))) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.size = 0) +
  facet_wrap(~ Genus, scales = 'free_x') +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  ylab("pp log10 CFU g-1 leaf") +
  xlab("plant tx") +
  scale_y_continuous(limits = c(6, 10), breaks = seq(6,10,1))

# let's put it together and see what the plot looks like:
ggarrange(plotlist = list(gh_breakdown_fam, gh_breakdown_clade, gh_breakdown_strain),common.legend = T, widths = c(0.265,0.455,1.125), ncol = 3) %>%
  ggsave(filename = here("figs/gh_structure_v2.pdf"), width = 6, height = 2.5)


# let's calculate the number of doublings that separate the plant treatments when summed across all strains/clades:
gh_pp7 <- reshape2::dcast(gh_pp6, rep + Genus ~ plant.tx, value.var = 'log_sum_med_ppcfu')
gh_pp7$log2diff <- log(10^(gh_pp7[,4]),2) - log(10^(gh_pp7[,3]),2)

# now generate plot which examines relative abundance across plant treatment types.
# this is to show the "community reshaping" point we make throughout the paper
# Do we use relative abundances? Do we plot change in rank abundances?

# generate community matrix:
gh_pp4$med_ppcfu <- round(10^gh_pp4$log_med_ppcfu)
gh_pp4d <- reshape2::dcast(gh_pp4, plant.tx + rep ~ strain.id, value.var = 'med_ppcfu')
gh_pp4d[,-c(1:2)] <- gh_pp4d[,-c(1:2)]/rowSums(gh_pp4d[,-c(1:2)])

# first, generate plot of relative abundance differences between treatments for each strain
gh_pp4d2 <- reshape2::melt(gh_pp4d, id.vars = c('plant.tx','rep'), variable.name = 'strain.id', value.name = 'freq')
gh_pp4d2$log_f <- log(gh_pp4d2$freq,10) # compute log frequency

# summarise with quantiles for each strain by treatment
gh_pp4d3 <- dplyr::group_by(gh_pp4d2, strain.id, plant.tx) %>% summarise(mu_log_f = mean(log_f),
                                                                         q0.025   = quantile(log_f, probs = c(0.025)),
                                                                         q0.25    = quantile(log_f, probs = c(0.25)),
                                                                         q0.50    = quantile(log_f, probs = c(0.5)),
                                                                         q0.75    = quantile(log_f, probs = c(0.75)),
                                                                         q0.975   = quantile(log_f, probs = c(0.975)))

# re-shape df and plot:
# merge separate subsets of the data!
gh_pp4d4 <- merge(gh_pp4d3[gh_pp4d3$plant.tx==1,], gh_pp4d3[gh_pp4d3$plant.tx==2,], by = 'strain.id', sort = F)

# merge with spp to plot by color:
gh_pp4d4 <- merge(gh_pp4d4, unique(gh1[,names(gh1) %in% c('strain.id','spp')]))

# now plot:
rfp1 <- ggplot(gh_pp4d4) +
  geom_errorbar(aes(x = q0.50.x, y = q0.50.y, ymin = q0.025.y, ymax = q0.975.y), col = "gray40", width = 0) +
  geom_errorbarh(aes(x = q0.50.x, y = q0.50.y, xmin = q0.025.x, xmax = q0.975.x), col = "gray40", height = 0) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = q0.50.x, y = q0.50.y, fill = spp), color = "black", pch = 21) +
  scale_fill_manual(values = c('white','black')) +
  ylab("med frequency, JA+") +
  xlab("med frequency, mock") +
  scale_y_continuous(limits = c(-5.1,0), breaks = seq(-5,0,1)) +
  scale_x_continuous(limits = c(-5.1,0), breaks = seq(-5,0,1)) +
  guides(fill=guide_legend(
    keywidth=0.2,
    keyheight=0.25,
    default.unit="cm")) + theme(legend.text=element_text(size=8))

ggsave(rfp1, filename = here("figs/gh_rel_freq_plot.pdf"), width = 3, height = 2.5)

# this plot shows re-shaping of compositional patterns that get reflected in diversity and evenness calculations to follow:
# now time to calculate the statistics:
# gh_pp4d4
# need to split this by rep and plant.tx
gh_pp5s <- split(gh_pp4d, gh_pp4d$rep)

# call up accessory function
quantize <- function(x){
  res <- c(mean = mean(x$value),
           quantile(x$value, probs=c(0.025,0.25,0.5,0.75,0.975)))
  res <- data.frame(stat = x$variable[1], t(res))
}

# need to lapply diversity calculator here:
gh_div <- lapply(gh_pp5s, calc_div)
gh_div2 <- do.call(rbind, gh_div)
gh_div3 <- reshape2::melt(gh_div2, id.var = 'rep')

# lapply quantize to each element
gh_div3b <- split(gh_div3, gh_div3$variable)
gh_div4 <- do.call(rbind, lapply(gh_div3b,quantize))

# now plot:
H <- ggplot(gh_div4[gh_div4$stat %in% c('H0','H1'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("median Shannon (H')") +
  xlab("plant tx") +
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2,0.5))

Hd <- ggplot(gh_div4[gh_div4$stat %in% c('Hd'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("delta H'") +
  xlab("") +
  geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  scale_y_continuous(limits = c(-2,0.1), breaks = seq(-2,0,0.5))

E <- ggplot(gh_div4[gh_div4$stat %in% c('E0','E1'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("median evenness (J')") +
  xlab("plant tx") +
  scale_y_continuous(limits = c(0,0.4), breaks = seq(0,0.4,0.1))

Ed <- ggplot(gh_div4[gh_div4$stat %in% c('Ed'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("delta J'") +
  xlab("") +
  geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  scale_y_continuous(limits = c(-0.3,0.05), breaks = seq(-0.3,0,0.1))

B <- ggplot(gh_div4[gh_div4$stat %in% c('SJ'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("community divergence (S-J)") +
  xlab("dataset") +
  #geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  scale_y_continuous(limits = c(0,0.30), breaks = seq(0,0.3,0.1))

# put it all together:
ggarrange(plotlist = list(H, Hd, E, Ed, B), ncol = 5, align = 'hv') %>% ggsave(filename = here("figs/gh_div_plots_v1.pdf"), width = 4, height= 2)

# calculate posterior summary statistics
Hd_bpv <- length(gh_div2$Hd[gh_div2$Hd < 0]) / length(gh_div2$Hd)
Ed_bpv <- length(gh_div2$Ed[gh_div2$Ed < 0]) / length(gh_div2$Ed)

#### RE-DO OF PLOTS WITH ALL DATA ####

## compute doublings and grab medians for clade-level data
gh_pp5b <- reshape2::dcast(gh_pp5, rep + spp ~ plant.tx, value.var = 'log_sum_med_ppcfu')
gh_pp5b$log2diff <- log(10^(gh_pp5b[,4]),2) - log(10^(gh_pp5b[,3]),2)

# here is where we calculate summary statistics at clade-level to report in main text:
quantile(gh_pp5b$log2diff[gh_pp5b$spp=='Psyr'], probs = c(0.025,0.25,0.5,0.75,0.975))
quantile(gh_pp5b$log2diff[gh_pp5b$spp=='Pfluo'], probs = c(0.025,0.25,0.5,0.75,0.975))

clade_doublings <- data.frame(rbind(t(quantile(gh_pp5b$log2diff[gh_pp5b$spp=='Psyr'], probs = c(0.025,0.25,0.5,0.75,0.975))),
                                    t(quantile(gh_pp5b$log2diff[gh_pp5b$spp=='Pfluo'], probs = c(0.025,0.25,0.5,0.75,0.975)))),
                              strain.id = c('Psyr','Pfluo'))
names(clade_doublings) <- names(coef_res)

## compute doublings and grab medians for genus-level data
gh_pp6b <- reshape2::dcast(gh_pp6, rep + Genus ~ plant.tx, value.var = 'log_sum_med_ppcfu')
gh_pp6b$log2diff <- log(10^(gh_pp6b[,4]),2) - log(10^(gh_pp6b[,3]),2)

# add to doublings data.frame and re-level factor to plot correctly:
genus_doublings <- data.frame(t(quantile(gh_pp6b$log2diff, probs = c(0.025,0.25,0.5,0.75,0.975))), strain.id='genus')
names(genus_doublings) <- names(coef_res)

# add back to coef_res
coef_res <- rbind(coef_res, clade_doublings, genus_doublings) # levels are already in correct order magically!


# now add abundance data to gh2 in order to produce heatmap
# need to take median across reps, call it log.med.cfu.g:
gh_pp6$Genus <- paste0('clade')
abund_genus <- dplyr::group_by(gh_pp6, Genus, plant.tx) %>% summarise(log.med.cfu.g = median(log_sum_med_ppcfu))
#abund_genus$plant.tx <- factor(abund_genus$plant.tx)

# do same per clade:
abund_clade <- dplyr::group_by(gh_pp5, spp, plant.tx) %>% summarise(log.med.cfu.g = median(log_sum_med_ppcfu))
#abund_clade$plant.tx <- factor(abund_clade$plant.tx)

# change factor names
names(abund_genus)[1] <- 'strain.id'
names(abund_clade)[1] <- 'strain.id'

# put it all back together:
gh3 <- dplyr::bind_rows(gh2[,-3], abund_clade, abund_genus)

# re-order strains:
gh3$strain.id <- factor(gh3$strain.id, levels = rev(c('clade','Psyr','Pfluo',strain_order)))

gh_p2 <- ggplot(gh3, aes(y = strain.id, x = factor(plant.tx))) +
  geom_tile(aes(fill = log.med.cfu.g)) +
  #scale_fill_gradient2(low = "steelblue", high = "darkorange2", midpoint = median(fam_sums$log_med_0.50)) +
  scale_fill_gradientn(colors = parula(64)) +
  #scale_fill_gradientn(colors = magma(20)) +
  #scale_fill_gradientn(colors = plasma(20)) +
  #scale_fill_gradientn(colors = viridis(20)) +
  #scale_fill_gradient(low = "white", high = "black") +
  theme(legend.position = 'top') + xlab("")

# now for the posterior of n.doublings different in JA versus MOCK:
doubling_diffs2 <- ggplot(coef_res) +
  geom_linerange(aes(x = strain.id,
                     #y = log2_med_med,
                     ymax = q0.975,
                     ymin = q0.025),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 0.5) +
  geom_linerange(aes(x = strain.id,
                     #y = log2_med_med,
                     ymax = q0.25,
                     ymin = q0.75),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 1.25) +
  geom_point(aes(x = strain.id, y = q0.50), pch = '|', size = 3, col = "white") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, col = "black", lty = 'dotted') +
  coord_flip() +
  theme(axis.text.y = element_blank()) + xlab("")

ggarrange(plotlist = list(gh_p2, doubling_diffs2),
          ncol = 2,
          widths = c(0.4,1),
          common.legend = T, align = 'h') %>%
  #ggsave(filename = here("figs/Fig4_bw_v3.pdf"), width = 3.125, height = 3.125)
ggsave(filename = here("figs/Fig4_parula_w_clades.pdf"), width = 3.125, height = 3.125)
