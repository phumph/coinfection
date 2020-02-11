# Script to compare gamma at treated versus un-treated stems of both damage classes.
# Compares values for all bacteria and then each of the 14 major taxa across the dataset.
# This plot will allow us to see the relative distribution of differences between treated and untreated stems
# and visually compare this to difference between damaged and undamaged classes at the same time.
# point will be to emphasize the small, if non-existent effect of the hormones versus the herbivore damage.

# header
source(file.path("./phy_header.R"))
source(file.path("./phy_functions.R"))

#### run for allBact ####

# load data
ASVs <- read.csv(file = file.path("../data/ASV_table_26-JUN-2018.csv"), row.names = 1)
bTAX <- read.csv(file = file.path("../data/bTAX_table_26-JUN-2018.csv"), row.names = 1)
mNP  <- read.csv(file = file.path("../data/NP_sample_data_final.csv"), row.names = 1)
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

# now bring in individual Family taxa for joint plotting:
DAT <- read.csv(file.path("../data/NPD_long_v1.csv"), header = T)
FAMS <- paste0(read.csv(file.path("../models/all_fams.csv"), header = F)$V1)
FAMS <- FAMS[!FAMS %in% fams_to_remove]

DAT <- DAT[DAT$bASV_count>0,]
DAT$log_ratio <- log(DAT$bASV_count / DAT$host)

# combined data for plotting:
focal_cols <- c('bASV_count','bASV','host','plot_cond','sub_plot_tx','stem_tx','plot_num','herb_dmg','Family')
DAT2 <- DAT[,focal_cols]
DAT2 <- droplevels(DAT2[DAT2$Family %in% FAMS,])
# sum across bASVs at the Family level:
DAT3 <- dplyr::group_by(DAT2, Family, plot_cond, sub_plot_tx, stem_tx, plot_num, herb_dmg) %>%
  summarise(Family_count = sum(bASV_count), host = mean(host))

DAT3$log_ratio <- log(DAT3$Family_count / DAT3$host)
NPD2 <- NPD[,c('plot_cond','sub_plot_tx','stem_tx','plot_num','herb_dmg','log_ratio')]
NPD2$Family <- 'allBacteria'
DAT4 <- rbind(data.frame(DAT3[,c('plot_cond','sub_plot_tx','stem_tx','plot_num','herb_dmg','log_ratio','Family')]),NPD2)

# decide which Families to include:
# Fams_to_exclude <- c('Cytophagaceae','Enterobacteriaceae','Hyphomicrobiaceae','Methylobacteriaceae','Sphingobacteriaceae')
# DAT5 <- dplyr::filter(DAT4, !Family %in% Fams_to_exclude)
DAT5 <- DAT4 # try without excluding any

# # generate plot for Mock:
# bvh_m <- ggplot(DAT5[DAT5$plot_cond=='M',], aes(x = sub_plot_tx, y = log_ratio, group = plot_num)) +
#   geom_hline(yintercept=0,col='black') +
#   facet_grid(Family ~ herb_dmg) +
#   geom_line(col='gray40',alpha=0.6) +
#   geom_point(col='gray40',alpha=0.6)
#
# bvh_ja <- ggplot(DAT5[DAT5$plot_cond=='JA',], aes(x = sub_plot_tx, y = log_ratio, group = plot_num)) +
#   geom_hline(yintercept=0,col='black') +
#   facet_grid(Family ~ herb_dmg) +
#   geom_line(col='gray40',alpha=0.6) +
#   geom_point(col='gray40',alpha=0.6)
#
# bvh_sa <- ggplot(DAT5[DAT5$plot_cond=='SA',], aes(x = sub_plot_tx, y = log_ratio, group = plot_num)) +
#   geom_hline(yintercept=0,col='black') +
#   facet_grid(Family ~ herb_dmg) +
#   geom_line(col='gray40',alpha=0.6) +
#   geom_point(col='gray40',alpha=0.6)

#### Combine counts per stem_tx and plot all at once ####
# re-level stem_tx:
DAT5$stem_tx <- factor(DAT5$stem_tx, levels = c('MC','JATX','SATX'))

# order the Families
fam_order <- c('allBacteria','Pseudomonadaceae','Enterobacteriaceae','Oxalobacteraceae',
               'Comamonadaceae','Sphingomonadaceae','Rhizobiaceae','Methylobacteriaceae',
               'Hyphomicrobiaceae','Streptococcaceae','Paenibacillaceae','Sphingobacteriaceae',
               'Flavobacteriaceae','Cytophagaceae','Microbacteriaceae')

DAT5$Family <- factor(DAT5$Family, levels = fam_order)

# creat composite factor to put all on same x-axis:
DAT5$tx_dmg <- factor(paste0(DAT5$stem_tx,'_',DAT5$herb_dmg))
DAT5$tx_dmg <- factor(DAT5$tx_dmg, levels = c('MC_0','JATX_0','SATX_0','MC_1','JATX_1','SATX_1'))

# calculate medians
DAT6 <- dplyr::group_by(DAT5, Family, tx_dmg, herb_dmg, stem_tx) %>% summarise(med_lr = median(log_ratio))

bvh_all <- ggplot(data = DAT5, aes(x = tx_dmg, y = log_ratio, col = factor(herb_dmg))) +
  geom_jitter(alpha=0.6, width=0.1) +
  geom_hline(yintercept=0, col='black') +
  facet_wrap( ~ Family, ncol = 5) +
  scale_color_manual(values = c('gray40','darkorange2')) +
  geom_crossbar(data = DAT6, aes(y = med_lr, x = tx_dmg, ymin = med_lr, ymax = med_lr), col='black',alpha=0.6) +
  xlab("sample-level treatment") + ylab("gamma")

ggsave(bvh_all, file = file.path("../figs/bvh_all_v1.pdf"),width = 6, height = 5)
