# script to examine patterns of bacterial abundance in hormone treatments at NP
# last updated: 2018-NOV-01 by PTH

# take chunks from run_brmsfit_bASV.R
library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

DAT <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv", header = T)
FAMS <- paste0(read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/models/all_fams.csv",header = F)$V1)
FAMS <- FAMS[!FAMS %in% c('Enterococcaceae','Rhizobiales_Incertae_Sedis')]

DAT <- DAT[DAT$bASV_count>0,]
DAT$log_ratio <- log(DAT$bASV_count / DAT$host)

options(mc.cores = parallel::detectCores())

# first, examine all bacteria model and generate plot:

DAT2 <- droplevels(dplyr::filter(DAT, Family %in% FAMS))
DAT3 <- unique(DAT2[,names(DAT2) %in% c('fungi','host','bac','herb_dmg','Row.names','plot_num','sub_plot_tx','sp_cond','plot_cond','tx_full','stem_tx')])

## let's plot fungus plot for kicks:
# DAT3$herb_dmg <- factor(DAT3$herb_dmg)
# ggplot(DAT3, aes(x = log(bac/host), col = herb_dmg, y = log(fungi/host))) + geom_point()

# generate plot of hormone treatment info:
# subset to only include non-herbivore damaged leaves
DAT4 <- dplyr::filter(DAT3, herb_dmg == 0)
DAT4$plot_cond <- factor(DAT4$plot_cond, levels = c('M','JA','SA'))
DAT4$bact <- 'all_bacteria'

# plot all bacteria
bac_horm_p1 <- ggplot(DAT4, aes(x = plot_cond, y = log(bac/host), col = sub_plot_tx)) +
  facet_wrap(~ bact) +
  #geom_jitter(width = 0.1,alpha=0.4) +
  geom_point(position = position_dodge(0.85), alpha = 0.4) +
  geom_boxplot(position = position_dodge(0.85), alpha = 0.4) + theme(legend.position = 'none') +
  scale_color_manual(values = c('gray40','dodgerblue')) +
  scale_y_continuous(limits = c(-11,5), breaks = seq(-10,5,5))# shrink facet text

# now let's look at each Family:
# sum up by family:
FAMDAT <- dplyr::group_by(DAT2, Row.names, herb_dmg, sub_plot_tx, plot_cond, plot_num, Family) %>%
  filter(herb_dmg == 0) %>%
  summarise(tot_fam_count = sum(bASV_count), host = unique(host)) %>%
  arrange(desc(tot_fam_count))

# sort FAMS by abundance, averaged across all samples:
FAM_order <- dplyr::group_by(FAMDAT, Family) %>% summarise(means = mean(log(tot_fam_count/host))) %>% arrange(desc(means))
FAMDAT$Family <- factor(FAMDAT$Family, levels = paste0(FAM_order$Family))
FAMDAT$plot_cond <- factor(FAMDAT$plot_cond, levels = c('M','JA','SA'))

# plot:
fam_horm_p1 <- ggplot(FAMDAT, aes(x = plot_cond, y = log(tot_fam_count/host), col = sub_plot_tx)) +
  #facet_grid(plot_cond ~ Family) +
  facet_wrap(~ Family, nrow = 2, ncol = 7) +
  #geom_jitter(width = 0.05, alpha=0.4, position = position_dodge(0.75)) +
  #stat_summary(geom = 'point', fun.y = "median", aes(col = sub_plot_tx), pch =  '-', position = position_dodge(0.75)) +
  geom_boxplot(alpha = 0.4, position = position_dodge(0.85)) +
  geom_point(alpha=0.4, position = position_dodge(0.85)) +
  scale_color_manual(values = c('gray40','dodgerblue')) + theme(legend.position = 'top') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_y_continuous(limits = c(-11,5), breaks = seq(-10,5,5))# shrink facet text
#fam_horm_p1

# trying separate plot incorportating plot-wise variation:
# FAMDAT$plot_num <- factor(FAMDAT$plot_num)
# ggplot(FAMDAT, aes(x = plot_num, y = log(tot_fam_count/host), col = sub_plot_tx)) +
#   geom_boxplot(alpha = 0.4, position = position_dodge(0.85)) +
#   geom_point(alpha=0.4, position = position_dodge(0.85)) +
#   facet_wrap(~ Family, nrow = 2, ncol = 7)


# output plots for supplemental info
ggsave(fam_horm_p1, filename = here("figs/hormone_effects_on_fams.pdf"),width = 8.5, height = 4.5)
ggsave(bac_horm_p1, filename = here("figs/hormone_effects_on_all_bact.pdf"),width = 1.4, height = 2)


#### repeat for EL ####

