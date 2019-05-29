# run_diversity_checks.R
# this script uses posterior predictive distribtutions of abundance for each bASV
# to generate estimates of community metrics such as diversity (H'), evenness (J'),
# and pairwise dissimilarity (S--J divergence) assuming independtly assembling communities

# this assumption is likely violated in practice but this gives us an initial glimpse into the magnitude of
# emergent community-level statistics arising from different 'responses'

# first step: load ppcfu files for specific model set
# run one by one; capture relevant columns into growing simulated-communities
# calculate focal diversity statistics per row
# parse into sets and examine sampling distributions of both the metrics and their summaries across the R=200 reps

# second step: Calculate Beta-diversity by averaging across leaves in the sample sets per rep
# calculate Shannon--Jensen divergence for each R=200 pairs of average probability vectors

# Capturing data routine:
library(here)
source(here("scripts/phy_header.R"))

#### FUNCTION DEFS ####

# function to grab bASV sample-level ppcfu data:
grab_ppcfu <- function(SET = 'EL', OUTFILE = here("models/div/"), base_mod_dir = here("models/"), prefix = 'ppcfu_*'){
  require(dplyr)
  require(reshape2)
  # grab all relevant files with prefix from focal dir:
  the_files <- Sys.glob(paste0(base_mod_dir,SET,'/ppcfu_data/',prefix))
  df <- data.frame()
  for (i in seq_along(the_files)){
    f <- read.csv(file = the_files[i])
    f <- dplyr::arrange(f, bASV, rep, herb_dmg, samp_id)
    f$pp_cfu <- round(f$pp_cfu * f$hits) # zero-out all posterior predicted prevalence == 0 hits
    f2 <- reshape2::dcast(f, rep + herb_dmg + samp_id ~ bASV, value.var = 'pp_cfu')
    rm(f)
    if(i==1){
      df <- f2
    } else{
      df <- cbind(df,f2[,-c(1:3)])
    }
  }
  write.csv(df, file = paste0(OUTFILE,'div_count_matr_',SET,'.csv'), row.names = F, quote = F)
}

# second function for doing similar operations on pre-computed medians across sample types for each of R reps:
grab_ppcfu_meds <- function(SET, OUTFILE = here("models/div/"), base_mod_dir = here("models/"), prefix = 'bASV-rep_*', output = FALSE){
  require(dplyr)
  require(reshape2)
  # grab all relevant files with prefix from focal dir:
  the_files <- Sys.glob(paste0(base_mod_dir,SET,'/ppcfu_data/',prefix))
  df <- data.frame()
  for (i in seq_along(the_files)){
    f <- read.csv(file = the_files[i])
    f <- dplyr::arrange(f, bASV, rep, herb_dmg)
    f2 <- reshape2::dcast(f, rep + herb_dmg ~ bASV, value.var = 'med_pp_cfu')
    rm(f)
    if(i==1){
      df <- f2
    } else{
      df <- cbind(df,f2[,-c(1:2)])
    }
  }
  if (output == TRUE){
    write.csv(df, file = paste0(OUTFILE,'div_med_matr_',SET,'.csv'), row.names = F, quote = F)
  }
  return(df)
}

# family-level results:
grab_ppcfu_fams <- function(SET, OUTFILE = here("models/div/"), base_mod_dir = here("models/"), prefix = 'fam-sum-rep_*', output = FALSE){
  require(dplyr)
  require(reshape2)
  # grab all relevant files with prefix from focal dir:
  the_files <- Sys.glob(paste0(base_mod_dir,SET,'/ppcfu_data/',prefix))
  df <- data.frame()
  for (i in seq_along(the_files)){
    f <- read.csv(file = the_files[i])
    f <- dplyr::arrange(f, rep, herb_dmg)
    f2 <- reshape2::dcast(f, rep + herb_dmg ~ Family, value.var = 'med_tot_leaf_pp_cfu')
    f2 <- data.frame(f2)
    rm(f)
    if(i==1){
      df <- data.frame(f2)
    } else{
      df <- cbind(df,f2[,-c(1:2)])
      names(df)[length(names(df))] <- last(names(f2))
    }
  }
  if (output == TRUE){
    write.csv(df, file = paste0(OUTFILE,'div_med_matr_Family_',SET,'.csv'), row.names = F, quote = F)
  }
  return(df)
}


# function to run inside lapply for calculating diversity statistics and beta diversity:
run_calc_div <- function(div, SET, BASE = 'bASV', do.plot = TRUE){
  div2 <- split(div,div$rep)
  divres <- lapply(div2,calc_div)
  divres2 <- do.call(rbind, divres)

  # generate quantiles of each column of divres2
  quantize <- function(x){
    res <- c(mean = mean(x$value),
             quantile(x$value,probs=c(0.025,0.25,0.5,0.75,0.975)))
    res <- data.frame(stat = x$variable[1],t(res))
  }

  divres3 <- reshape2::melt(divres2, id.var = 'rep')

  # lapply quantize to each element
  divres3b <- split(divres3, divres3$variable)
  divres_stats <- do.call(rbind,lapply(divres3b,quantize))
  divres_stats$dataset <- SET

  # now plot dist
  sh.p1 <- ggplot(dplyr::filter(divres3, variable %in% c('H0','H1')), aes(y = value, x = variable)) +
    geom_jitter(alpha = 0.2, col='gray40', width = 0.1) +
    geom_boxplot(alpha=0.4, col='gray40') + xlab("") + ylab("Shannon (H')") #+
  #scale_y_continuous(limits = c(6, 8), breaks = seq(6,8,0.5))

  sh.dp1 <- ggplot(dplyr::filter(divres3, variable %in% c('Hd')), aes(y = value, x = variable)) +
    geom_jitter(alpha = 0.2, col='gray40', width = 0.1) +
    geom_boxplot(alpha=0.4, col='gray40') + xlab("") + ylab("difference in H'") +
    geom_hline(yintercept = 0, lty = 3, col = 'gray60')
  #scale_y_continuous(limits = c(6, 8), breaks = seq(6,8,0.5))

  ev.p1 <- ggplot(dplyr::filter(divres3, variable %in% c('E0','E1')), aes(y = value, x = variable), col = 'gray40') +
    geom_jitter(alpha = 0.2, col='gray40', width = 0.1) +
    geom_boxplot(alpha=0.4, col='gray40') + xlab("") + ylab("evenness (J')") #+
    #scale_y_continuous(limits = c(0.75, 1.0), breaks = seq(0.75,1.0,0.05))

  ev.dp1 <- ggplot(dplyr::filter(divres3, variable %in% c('Ed')), aes(y = value, x = variable), col = 'gray40') +
    geom_jitter(alpha = 0.2, col='gray40', width = 0.1) + ylab("difference in J'") +
    geom_boxplot(alpha=0.4, col='gray40') + xlab("") +
    geom_hline(yintercept = 0, lty = 3, col = 'gray60')
  #scale_y_continuous(limits = c(0.75, 1.0), breaks = seq(0.75,1.0,0.05))

  # now plot beta div and put all in one plot:
  beta.p1 <- ggplot(dplyr::filter(divres3, variable %in% c('SJ')), aes(y = value, x = variable), col = 'gray40') +
    geom_jitter(alpha = 0.2, col='gray40', width = 0.1) +
    geom_boxplot(alpha=0.4, col='gray40') + xlab("") + ylab("beta diversity (S-J)") #+
    #scale_y_continuous(limits = c(0.1, 0.3), breaks = seq(0.1,0.3,0.05))
  if(do.plot==TRUE){
    ggarrange(plotlist = list(sh.p1,ev.p1,sh.dp1,ev.dp1,beta.p1), ncol = 5, widths = c(1,1,0.75,0.75,0.75),align='hv') %>%
      ggsave(filename = paste0(here("models/div/"),'div_plots_',BASE,'_',SET,'.pdf'),width = 6.5, height = 2.5)
  }
  return(divres_stats)
}

#### MAINS ####
ELdiv <- grab_ppcfu_meds(SET='EL', output = TRUE)
NPdiv <- grab_ppcfu_meds(SET='NP', output = TRUE)

run_calc_div(ELdiv, SET='EL', BASE = 'bASV')
run_calc_div(NPdiv, SET='NP', BASE = 'bASV')

# now do this again but at family level:
ELfam <- grab_ppcfu_fams(SET='EL', output = TRUE)
NPfam <- grab_ppcfu_fams(SET='NP', output = TRUE)

ELdivres <- run_calc_div(ELfam, SET='EL', BASE = 'Family', do.plot = FALSE)
NPdivres <- run_calc_div(NPfam, SET='NP', BASE = 'Family', do.plot = FALSE)

divres_all <- rbind(ELdivres,NPdivres)

# now plot:
H <- ggplot(divres_all[divres_all$stat %in% c('H0','H1'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("median Shannon (H')") +
  xlab("herbivory") +
  scale_y_continuous(limits = c(1.90,3.5), breaks = seq(2.0,3.5,0.5)) +
  facet_wrap(~ dataset)

E <- ggplot(divres_all[divres_all$stat %in% c('E0','E1'),]) +
  geom_linerange(aes(x = stat, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = stat, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = stat, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("median evenness (J')") +
  xlab("herbivory") +
  scale_y_continuous(limits = c(0.5,0.9), breaks = seq(0.5,0.9,0.1)) +
  facet_wrap(~ dataset)

Hd <- ggplot(divres_all[divres_all$stat %in% c('Hd'),]) +
  geom_linerange(aes(x = dataset, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = dataset, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = dataset, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("delta H'") +
  xlab("dataset") +
  facet_wrap(~stat) +
  geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  scale_y_continuous(limits = c(-1.2,0.0), breaks = seq(-1.2,0.0,0.4))

Ed <- ggplot(divres_all[divres_all$stat %in% c('Ed'),]) +
  geom_linerange(aes(x = dataset, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = dataset, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = dataset, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("delta J'") +
  xlab("dataset") +
  geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  facet_wrap(~stat) +
  scale_y_continuous(limits = c(-0.3,0.0), breaks = seq(-0.3,0.0,0.1))

B <- ggplot(divres_all[divres_all$stat %in% c('SJ'),]) +
  geom_linerange(aes(x = dataset, ymin = X2.5., ymax = X97.5.), size = 0.75, col = "gray27") +
  geom_linerange(aes(x = dataset, ymin = X25., ymax = X75.), size = 1.5, col = "gray27") +
  geom_point(aes(x = dataset, y = X50.), size = 5, shape = 95, col = "white") +
  theme_phy1() +
  ylab("community divergence (S-J)") +
  xlab("dataset") +
  #geom_hline(yintercept = 0, lty= 3, col = 'gray40') +
  facet_wrap(~stat) +
  scale_y_continuous(limits = c(0,0.30), breaks = seq(0,0.3,0.1))

# put it all together:
ggarrange(plotlist = list(H,Hd,B,E,Ed), ncol = 3, nrow = 2, widths = c(1,0.66,0.66), heights = c(1,1)) %>%
  ggsave(filename = here("figs/Fig2_diversity_v1.pdf"), width = 3.5, height = 4)


# go back and plot the change in frequencies across all bASVs and Families between herbivory types using computed pp medians:
# this plot will be for supplemental and shows the data basis on which the diversity calculations are made/justified

x <- ELfam
plot_rel_freqs <- function(x){
  x2 <- cbind(x[,c(1:2)], x[,-c(1:2)] / rowSums(x[,-c(1:2)])) # generate (log) frequency table
  x3 <- reshape2::melt(x2, id.vars = c('herb_dmg','rep'), variable.name = 'Family', value.name = 'freq')
  x3$log_f <- log(x3$freq,10) # compute log frequency
  # calculate quantiles
  x4 <- dplyr::group_by(x3, Family, herb_dmg) %>% summarise(mu_log_f = mean(log_f),
                                                                           q0.025 = quantile(log_f, probs = c(0.025)),
                                                                           q0.25 = quantile(log_f, probs = c(0.25)),
                                                                           q0.50 = quantile(log_f, probs = c(0.5)),
                                                                           q0.75 = quantile(log_f, probs = c(0.75)),
                                                                           q0.975 = quantile(log_f, probs = c(0.975)))
  # do some more wrangling to get into plotting shape:
  x4b <- merge(x4[x4$herb_dmg==0,], x4[x4$herb_dmg==1,], by = 'Family', sort = F)
  # now plot:
  coln <- length(unique(x4b$Family))
  the_plot <- ggplot(x4b) +
    geom_errorbar(aes(x = q0.50.x, y = q0.50.y, ymin = q0.025.y, ymax = q0.975.y), col = "gray40") +
    geom_errorbarh(aes(x = q0.50.x, y = q0.50.y, xmin = q0.025.x, xmax = q0.975.x), col = "gray40") +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x = q0.50.x, y = q0.50.y, fill = Family), color = "black", pch = 21) +
    scale_fill_manual(values = plasma(coln)) +
    ylab("med frequency, DMG +") +
    xlab("med frequency, DMG -") +
    scale_x_continuous(limits = c(-3.5,0), breaks = seq(-3,0,1)) +
    scale_y_continuous(limits = c(-3.5,0), breaks = seq(-3,0,1)) +
    guides(fill=guide_legend(
      keywidth=0.2,
      keyheight=0.25,
      default.unit="cm")) + theme(legend.text=element_text(size=8))

  return(the_plot)
}

ELfam_plot <- plot_rel_freqs(ELfam) + ggtitle("EL")
NPfam_plot <- plot_rel_freqs(NPfam) + ggtitle("NP")

ggarrange(plotlist = list(ELfam_plot, NPfam_plot), common.legend = T, align = 'hv', nrow = 2, legend = 'right') %>%
  ggsave(filename = here("figs/rel_freq_plot_Fams.pdf"), width = 3.5, height = 4.5)

# done!
