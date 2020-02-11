# script to prepare data for plotting for Family-level and bASV-level data for Figs. 2 and relevant appendixes
# Last updated: 2018-11-09 PTH

# Updates
# 2018-08-15:
  # 1. Produced median predicted CFU for damage type for all included Families. Do also for all-bacteria models.
  # 2. Produced log2-fold difference (i.e., number of cell divisions difference)
  # 3. Did this for both EL and NP model sets
  # 4. Produced Fig. 2a (all-bacteria and Family-level plot)

# 2018-11-09
  # 1. Produce supplemental figure showing EL versus NP comparisons of (a) median abundance (+/- damage) and (b) doublings diff.

##########################

# load libraries and headers
source(file.path("./phy_header.R"))
source(file.path("./phy_functions.R"))

# load all relevant data for Family and bASV-level data at rep-level.
FAM_FILES_DIR_EL <- "../models/EL/ppcfu_data/"
FAM_FILES_DIR_NP <- "../models/NP/ppcfu_data/"
FAM_PREFIX <- 'fam-sum-all*'
FAM_REP_PREFIX <- 'fam-sum-rep*'

ALL_FILES_DIR <- "../models/all-bacteria/"
ALL_PREFIX <- 'sum-all_'
ALL_REP_PREFIX <- 'sum-all-rep_'

# family-level taxonomy:
#bTAX <- read.csv(here("data/bTAX_table_final.csv"))
bTAX <- read.csv(file.path("../data/bTAX_table_26-JUN-2018.csv"))
FAMS <- paste0(read.csv(file.path("../models/all_fams.csv"), header = F)[,1])
FAMS <- FAMS[!FAMS %in% c('Rhizobiales_Incertae_Sedis','Enterococcaceae')]
fam_tax <- unique(bTAX[,c('Phylum','Class','Order','Family')]) %>% dplyr::filter(Family %in% FAMS) %>% arrange(desc(Phylum), desc(Class), desc(Order), desc(Family))

# for summary plots, put EL and then NP as columns and then the log-2-fold difference (i.e., number of doublings) that separate the two categories.
# Better to plot the heatmap of the log-10 CFU/g and then the posterior of the number of doublings (posterior) as the linechart. That could look quite nice..

# first, load fam-sum-all data to grab mean or median per leaf category
fam_sum_files_EL <- Sys.glob(paste0(FAM_FILES_DIR_EL,FAM_PREFIX))
fam_sum_files_NP <- Sys.glob(paste0(FAM_FILES_DIR_NP,FAM_PREFIX))

allB_sum_EL <- Sys.glob(paste0(ALL_FILES_DIR,ALL_PREFIX,'EL*'))
allB_sum_NP <- Sys.glob(paste0(ALL_FILES_DIR,ALL_PREFIX,'NP*'))

fam_rep_files_EL <- Sys.glob(paste0(FAM_FILES_DIR_EL,FAM_REP_PREFIX))
fam_rep_files_NP <- Sys.glob(paste0(FAM_FILES_DIR_NP,FAM_REP_PREFIX))

allB_rep_EL <- Sys.glob(paste0(ALL_FILES_DIR,ALL_REP_PREFIX,'EL*'))
allB_rep_NP <- Sys.glob(paste0(ALL_FILES_DIR,ALL_REP_PREFIX,'NP*'))

fam_sums <- do.call(rbind, lapply(as.list(fam_sum_files_EL), read.csv))
fam_sums$dataset <- 'EL'

fam_sums2 <- do.call(rbind, lapply(as.list(fam_sum_files_NP), read.csv))
fam_sums2$dataset <- 'NP'

# add in allBact sets:
all_sum <- do.call(rbind, lapply(as.list(c(allB_sum_EL,allB_sum_NP)), read.csv))
all_sum$Family <- 'all_bacteria'
all_sum$dataset <- c(rep('EL',2),rep('NP',2))

# stitch allBact, EL, and NP sums data together:
fam_sums <- rbind(all_sum, fam_sums, fam_sums2)

# re-define factor order of Family
#fam_sums$Family <- factor(fam_sums$Family, levels = c('all_bacteria',paste0(fam_tax$Family)))

# ditch Enterococcaceae
fam_sums <- dplyr::filter(fam_sums, Family %in% c('all_bacteria',FAMS)) %>% arrange(log_mu_mu)
#fam_sums$Family <- factor(fam_sums$Family, levels = unique(paste0(fam_sums$Family)))
fam_sums$Family <- factor(fam_sums$Family, levels = c(rev(unique(paste0(fam_tax$Family))),'all_bacteria'))

# re-level based on NP Families sorted:
# need to order by phylogeny, prioritzed by abundance

# now construct log2-fold differences between dmg classes for side-by-side plotting:
fam_sums <- dplyr::arrange(fam_sums, dataset, Family, herb_dmg)

## below this is probably obsolete:
# fs1 <- reshape2::dcast(fam_sums, Family + dataset ~ herb_dmg, value.var = c('log_med_0.50'))
# fs2 <- reshape2::dcast(fam_sums, Family + dataset ~ herb_dmg, value.var = c('log_med_0.025'))
# fs3 <- reshape2::dcast(fam_sums, Family + dataset ~ herb_dmg, value.var = c('log_med_0.975'))
# fs4 <- reshape2::dcast(fam_sums, Family + dataset ~ herb_dmg, value.var = c('log_med_0.75'))
# fs5 <- reshape2::dcast(fam_sums, Family + dataset ~ herb_dmg, value.var = c('log_med_0.25'))
# fam_sums_cast <- cbind(fs1,fs2[,c(3:4)],fs3[,c(3:4)],fs4[,c(3:4)],fs5[,c(3:4)])
# names(fam_sums_cast) <- c('Family','dataset','log_med_0.50_0','log_med_0.50_1','log_med_0.025_0','log_med_0.025_1','log_med_0.975_0','log_med_0.975_1','log_med_0.75_0','log_med_0.75_1','log_med_0.25_0','log_med_0.25_1')
#
# # now calculate:
# fam_sums_cast <- dplyr::mutate(fam_sums_cast,
#                                log2_med_med = log(10^(log_med_0.50_1 - log_med_0.50_0),2),
#                                log2_med_0.975 = log(10^(log_med_0.975_1 - log_med_0.975_0),2),
#                                log2_med_0.025 = log(10^(log_med_0.025_1 - log_med_0.025_0),2),
#                                log2_med_0.75 = log(10^(log_med_0.75_1 - log_med_0.75_0),2),
#                                log2_med_0.25 = log(10^(log_med_0.25_1 - log_med_0.25_0),2))
#
# fam_sums_cast$Family <- factor(fam_sums_cast$Family, levels = c(rev(unique(paste0(fam_tax$Family))),'all_bacteria'))

#### DIFFERENCE IN DOUBLINGS ####
# add rep-level information for posterior difference in doublings between damage classes:
fam_reps <- do.call(rbind, lapply(as.list(fam_rep_files_EL), read.csv))
fam_reps$dataset <- 'EL'

fam_reps2 <- do.call(rbind, lapply(as.list(fam_rep_files_NP), read.csv))
fam_reps2$dataset <- 'NP'

# add in allBact sets:
all_rep <- do.call(rbind, lapply(as.list(c(allB_rep_EL, allB_rep_NP)), read.csv))
all_rep$Family <- 'all_bacteria'
all_rep$dataset <- c(rep('EL',400),rep('NP',400))
all_rep$prev <- NA
all_rep$Family <- 'all_bacteria'
all_rep$dataset <- c(rep('EL',400),rep('NP',400))

# stitch allBact, EL, and NP sums data together:
fam_reps <- rbind(all_rep, fam_reps, fam_reps2)

fam_reps <- transform(fam_reps,
                      herb_dmg = factor(herb_dmg),
                      Family = factor(Family),
                      dataset = factor(dataset))

## let's see if this is the same as if we were to subtract the two distributions with rep-level data:
fam_reps <- dplyr::filter(fam_reps, Family %in% c('all_bacteria',FAMS))

# need to calculate difference between distributions:
# define function to calculate quantile of log2-fold difference:
log2_diff_calc <- function(x, the_val = 'log_med_tot_leaf_pp_cfu', tax_col = 'Family', R = 200){
  require(dplyr)
  taxa <- paste0(unique(x[,tax_col]))
  datasets <- paste0(unique(x$dataset))
  res <- data.frame()
  for (f in seq_along(taxa)){
    for (d in seq_along(datasets)){
      x2 <- x[x[,tax_col] == taxa[f] & x[,'dataset'] == datasets[d],]
      test.mat <- as.matrix(data.frame(x1 = x2[x2$herb_dmg==1,the_val],
                                       x0 = x2[x2$herb_dmg==0,the_val]))
      x.diffs <- log(10^(as.vector(outer(test.mat[,1],test.mat[,2],'-'))),2)
      qs <- quantile(x.diffs,
                     probs = c(0.025,0.25,0.5,0.75,0.975))
      res <- rbind(res,
                   data.frame(taxa = taxa[f],
                              dataset = datasets[d],
                              t(as.vector(qs)))
      )
    }
  }
  names(res) <- c(paste0(tax_col),'dataset','q0.025','q0.25','q0.50','q0.75','q0.975')
  return(res)
}


fam_rep_res <- log2_diff_calc(fam_reps)
fam_rep_res$Family <- factor(fam_rep_res$Family, levels = c(rev(unique(paste0(fam_tax$Family))),'all_bacteria'))

# make heatmap of log10 abundance broken down by dataset and leaf type:
fp1 <- ggplot(fam_sums, aes(y = Family, x = factor(herb_dmg))) + geom_tile(aes(fill = log_med_0.50)) +
  #scale_fill_gradient2(low = "steelblue", high = "darkorange2", midpoint = median(fam_sums$log_med_0.50)) +
  #scale_fill_gradientn(colors = viridis(20)) +
  #scale_fill_gradientn(colors = cividis(20)) +
  #scale_fill_gradientn(colors = plasma(20)) +
  #scale_fill_gradientn(colors = kovesi.diverging_linear_bjr_30_55_c53(64)) +
  scale_fill_gradientn(colors = parula(64)) +
  #scale_fill_gradientn(colors = ocean.haline(64)) +
  facet_wrap(~ dataset) +
  theme(legend.position = 'top') + xlab("")

# make linerange plot of doublings:
doublings_fam_level <- ggplot(fam_rep_res) +
  geom_linerange(aes(x = Family,
                     #y = log2_med_med,
                     ymax = q0.975,
                     ymin = q0.025),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 0.5) +
  geom_linerange(aes(x = Family,
                     #y = log2_med_med,
                     ymax = q0.25,
                     ymin = q0.75),
                 #col = 'steelblue',
                 col = 'gray40',
                 size = 1.25) +
  geom_point(aes(x = Family, y = q0.50), pch = '|', size = 3, col = "white") +
  theme_phy1() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, col = "black", lty = 'dotted') +
  facet_wrap(~ dataset) + coord_flip() + ylab("median doublings gained from herbivory") + xlab("") +
  scale_y_continuous(limits = c(-2,7), breaks = seq(-2,6,2)) +
  theme(axis.text.y = element_blank())

ggarrange(plotlist = list(fp1, doublings_fam_level), ncol = 2, widths = c(0.85,1.5), common.legend = T, align = 'h') %>%
  ggsave(filename = file.path("../figs/Fig2a_with_all_bact_parula.pdf"), width = 6, height = 2.75)


#### SI Figure-making ####

## 1. produce xyplot of EL v NP median abundances +/- credible posterior intervals for both damage classes
# use fam_sums; need to cast for plotting:
# filter and merge:
tmpEL <- dplyr::filter(fam_sums, dataset == 'EL')
tmpNP <- dplyr::filter(fam_sums, dataset == 'NP')
fam_sums_cast <- merge(tmpEL, tmpNP, by = c('Family','herb_dmg'))

# define number of colors
coln <- length(unique(fam_sums_cast$Family))

# make plot:
med_comp_sites <- ggplot(fam_sums_cast, aes(x = log_med_0.50.x, y = log_med_0.50.y, fill = Family)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = log_med_0.025.y, ymax = log_med_0.975.y), size = 0.5) +
  geom_errorbarh(aes(xmin = log_med_0.025.x, xmax = log_med_0.975.x), size = 0.5) +
  geom_errorbar(aes(ymin = log_med_0.25.y, ymax = log_med_0.75.y), size = 1) +
  geom_errorbarh(aes(xmin = log_med_0.25.x, xmax = log_med_0.75.x), size = 1) +
  geom_point(size = 2, pch = 21) +
  facet_wrap(~ herb_dmg) +
  scale_fill_manual(values = plasma(coln)) +
  theme(legend.position = 'none') +
  # guides(fill=guide_legend(
  #   keywidth=0.2,
  #   keyheight=0.25,
  #   default.unit="cm")) + theme(legend.text=element_text(size=8)) +
  scale_x_continuous(limits = c(3.8,9.4), breaks = seq(4,9,1)) +
  scale_y_continuous(limits = c(3.8,9.4), breaks = seq(4,9,1)) +
  xlab("med pp. log10 CFU g-1, site EL") +
  ylab("med pp. log10 CFU g-1, site NP")

## 2. now compare relative changes (i.e., number of doublings difference between DMG classes):
# use fam_rep_res
tmpEL2 <- dplyr::filter(fam_rep_res, dataset == 'EL')
tmpNP2 <- dplyr::filter(fam_rep_res, dataset == 'NP')
fam_rep_res_cast <- merge(tmpEL2, tmpNP2, by = c('Family'))

# make plot:
doub_comp_sites <- ggplot(fam_rep_res_cast, aes(x = q0.50.x, y = q0.50.y, fill = Family)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_errorbar(aes(ymin = q0.025.y, ymax = q0.975.y), size = 0.5) +
  geom_errorbarh(aes(xmin = q0.025.x, xmax = q0.975.x), size = 0.5) +
  geom_errorbar(aes(ymin = q0.25.y, ymax = q0.75.y), size = 1) +
  geom_errorbarh(aes(xmin = q0.25.x, xmax = q0.75.x), size = 1) +
  geom_point(size = 2, pch = 21) +
  scale_fill_manual(values = plasma(coln)) +
  guides(fill=guide_legend(
    keywidth=0.2,
    keyheight=0.25,
    default.unit="cm")) + theme(legend.text=element_text(size=8)) +
  scale_x_continuous(limits = c(-2,6.5), breaks = seq(-2,6,1)) +
  scale_y_continuous(limits = c(-2,6.5), breaks = seq(-2,6,1)) +
  xlab("log2 difference, site EL") +
  ylab("log2 difference, site NP")

# put together, make relatively big for ease of viewing in SI:
# top panels
ggsave(med_comp_sites, filename = file.path("../figs/SI_med_comp.pdf"), width = 5, height = 2.5)

# bottom panels, with legend:
ggsave(doub_comp_sites, filename = file.path("../figs/SI_doub_comp.pdf"), width = 4, height = 2.5)

# ggarrange(plotlist = list(med_comp_sites,doub_comp_sites),common.legend = T,ncol = 2, widths = c(1.8,1), legend = 'right') %>%
#   ggsave(filename = here("figs/SI_med_doub_comp.pdf"), width = 6.5, height = 2)

#### OBSOLETE ####
# now moving on to Pseudomonas:
# for heatmap for median of median log abundance between leaf classes:
# Ps_sum_EL <- read.csv(here("models/EL/ppcfu_data/bASV-all_Pseudomonadaceae.csv"))
# Ps_sum_NP <- read.csv(here("models/NP/ppcfu_data/bASV-all_Pseudomonadaceae.csv"))
# Ps_sum_EL$dataset <- 'EL'
# Ps_sum_NP$dataset <- 'NP'
# Ps_sum <- rbind(Ps_sum_EL,Ps_sum_NP)
#
# # need to order the bASVs in terms of putative taxonomy:
# Ps_tree_order <- rev(paste0(read.csv(here("blast-res/bASVs/out/bASV_trees/Ps_tree_order_v1.txt"),F)[,1]))
# Ps_sum <- dplyr::filter(Ps_sum, bASV %in% Ps_tree_order)
# Ps_sum$bASV <- factor(Ps_sum$bASV, levels = Ps_tree_order)
#
# # try plot:
# ps1 <- ggplot(Ps_sum, aes(y = bASV, x = factor(herb_dmg))) + geom_tile(aes(fill = log_med_0.50)) +
#   #scale_fill_gradient2(low = "steelblue", high = "darkorange2", midpoint = median(fam_sums$log_med_0.50)) +
#   #scale_fill_gradientn(colors = viridis(20)) +
#   scale_fill_gradientn(colors = plasma(20)) +
#   facet_wrap(~ dataset) +
#   theme(legend.position = 'top') + xlab("")
#
#
# # for log2-diff plot:
# # load up rep-level bASV sim data
# Ps_dat_EL <- read.csv(here("models/EL/ppcfu_data/bASV-rep_Pseudomonadaceae.csv"))
# Ps_dat_NP <- read.csv(here("models/NP/ppcfu_data/bASV-rep_Pseudomonadaceae.csv"))
# Ps_dat_EL$dataset <- 'EL'
# Ps_dat_NP$dataset <- 'NP'
# Ps_dat <- rbind(Ps_dat_EL,Ps_dat_NP)
# Ps_dat <- dplyr::filter(Ps_dat, bASV %in% Ps_tree_order)
#
# # next step is to run through each bASV and calculate log2-diff between damaged and un-damaged:
# Ps_dat_res <- log2_diff_calc(Ps_dat, the_val = 'log_med_pp_cfu', tax_col = 'bASV', R = 200)
# Ps_dat_res$bASV <- factor(Ps_dat_res$bASV, levels = Ps_tree_order)
#
# # make linerange plot of doublings:
# doublings_Ps <- ggplot(Ps_dat_res) +
#   geom_linerange(aes(x = bASV,
#                      #y = log2_med_med,
#                      ymax = q0.975,
#                      ymin = q0.025),
#                  #col = 'steelblue',
#                  col = 'gray40',
#                  size = 0.5) +
#   geom_linerange(aes(x = bASV,
#                      #y = log2_med_med,
#                      ymax = q0.25,
#                      ymin = q0.75),
#                  #col = 'steelblue',
#                  col = 'gray40',
#                  size = 1.25) +
#   geom_point(aes(x = bASV, y = q0.50), pch = '|', size = 3, col = "white") +
#   theme_phy1() +
#   theme(legend.position = "none") +
#   geom_hline(yintercept = 0, col = "black", lty = 'dotted') +
#   facet_wrap(~ dataset) + coord_flip() + ylab("median doublings gained from herbivory") + xlab("") +
#   scale_y_continuous(limits = c(-2,7), breaks = seq(-2,6,2)) +
#   theme(axis.text.y = element_blank())
#
# ggarrange(plotlist = list(ps1, doublings_Ps), ncol = 2, widths = c(0.85,1.5), common.legend = T, align = 'h') %>% ggsave(filename = here("figs/Fig2c.pdf"), width = 7.5, height = 5)
#
#
# # Hmm.. Ok this probably belongs as a supplemental table. But maybe I can summarise the information somehow?
# # Let's try to do this.. just map boolean factor onto the Ps_sum data..
# Ps_tax_info <- read.csv(here("blast-res/bASVs/out/bASV_trees/Ps_guide_for_Fig2.txt"))
# names(Ps_tax_info)[1] <- 'bASV'
# Ps_tax_info <- dplyr::filter(Ps_tax_info, bASV %in% Ps_tree_order)
#
# Ps_sum$group <- 0
# Ps_sum$group[Ps_sum$bASV %in% Ps_tax_info$bASV[Ps_tax_info$Ps==1]] <- 'syringae'
# Ps_sum$group[Ps_sum$bASV %in% Ps_tax_info$bASV[Ps_tax_info$Pf==1]] <- 'fluorescens'
#
# ggplot(Ps_sum, aes(x = factor(herb_dmg), y = log_mu_mu, col = group)) + geom_boxplot(position = position_dodge(1)) + geom_jitter(position = position_dodge(1)) + facet_wrap(~ dataset)
#
# # not really much to see here.. it's basically a single P. syringae bASV that dominates
#
# #### OBSOLETE AND/OR TESTING ####
#
# # load up tax-guide for Pseudomonas
# # flag P. syringae, P. fluorescens, or other strain names by group
# # Ps_guide <- read.csv(here("blast-res/bASVs/out/Pseudomonadaceae_tax-guide.csv"))
# # Ps_guide2 <- dplyr::filter(Ps_guide, V5 > 99)
# #
# # Ps_guide$unique_id[grep('fluorescens', Ps_guide$species.group)]
# # a1 <- unique(Ps_guide$unique_id[grep('fluorescens', Ps_guide$species.group)])
# # a2 <- unique(Ps_guide$unique_id[grep('syringae', Ps_guide$species.group)])
# # table(a1 %in% a2)
#
# # calculate log-2-fold difference posterior between DMG+ and DMG- leaves
# bASV_log2 <- reshape2::dcast(bASV_rep, bASV + rep ~ herb_dmg, value.var = 'med_pp_cfu')
# names(bASV_log2) <- c('bASV','rep','herb_dmg_0', 'herb_dmg_1')
#
# # calculate log-2-fold difference
# bASV_log2$l2d <- log(bASV_log2$herb_dmg_1 / bASV_log2$herb_dmg_0,2)
# ggplot(bASV_log2, aes(x = l2d, col = bASV)) + geom_density()
# ggplot(bASV_rep, aes(x = factor(herb_dmg), y = log_med_pp_cfu, col = bASV)) + geom_boxplot()
