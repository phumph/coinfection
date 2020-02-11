# script to examine efficacy of surface sterilization procedure
# init: 2019-AUG-28
# PTH

# ------ #
# header #
# ------ #

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)

# ---------------- #
# data preparation #
# ---------------- #

# first, load up output from DADA2 which contains samples from Diversity plots
data.path <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/16S_seq/" # change to local data path for repeating this analysis.
outdir <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/figs/"
TAB <- read.csv(paste0(data.path,"seqs_through_pipe_v1.csv"),T)               # output from 00_DADA2_pipeline
ASV <- read.csv(paste0(data.path,"dada2_ASV_table_v3.csv"),T, row.names = 1)  # output from 00_DADA2_pipeline
TAX <- read.csv(paste0(data.path,"dada2_taxa_Silva_wSpp_uniques_v1.csv"),T)   # output from 00_DADA2_pipeline

# grab subset of ASV table corresponding to DIV plots
# this is the set which got the non-sterilization control treatment
# for testing efficacy of this procedure
DIV <- ASV[grep('D',row.names(ASV)),]
DIV <- DIV[!grepl('NPB',row.names(DIV)), ]

DIV$CONTROL <- 0
DIV$CONTROL[grep('C',row.names(DIV))] <- 1
DIV$sample_id <- sapply(row.names(DIV), function(x) gsub('C','',x))
DIV$row.sum <- rowSums(DIV[, !names(DIV) %in% c('CONTROL','sample_id')])
DIVC <- DIV[DIV$sample_id %in% (DIV$sample_id[DIV$CONTROL==1]), ]

# select only those with paired data:
paired_samples <- rowSums(table(DIVC$sample_id, DIVC$CONTROL))
to_exclude <- names(paired_samples)[paired_samples < 2]

if (to_exclude > 0) {
  DIVC <- DIVC[!DIVC$sample_id %in% to_exclude, ]
}

divc_hits <- DIVC %>%
  dplyr::select(-CONTROL, -row.sum) %>%
  dplyr::summarise_if(is.numeric, sum)

# get those with no hits
no_hits <- names(divc_hits)[divc_hits==0]
TAB2 <- dplyr::filter(TAB, X %in% row.names(DIVC))
TAX2 <- dplyr::filter(TAX, !X %in% no_hits)
rownames(TAX2) <- TAX2[,'X']

# find cp and mt:
cp <- grep("Chloroplast", TAX2[,'unique.id'])
mt <- grep("Mitochondria", TAX2[,'unique.id'])

###
### APPROACH 1: assume cp and mt seqs are correct; sum and calculate gamma
###

# sum host-derived reads in DIVC
DIVC$host <- rowSums(DIVC[, names(DIVC) %in% paste0(TAX2[cp,'X']) | names(DIVC) %in% TAX2[mt,'X']])
DIVC$all  <- rowSums(DIVC[, !names(DIVC) %in% tail(names(DIVC), 4)])
DIVC$bact <- DIVC$all - DIVC$host
DIVC$gamma <- log(DIVC$bact / DIVC$host)

# do Bayesian model to report posterior of coefficient estimate;
# re-sample logCFU model to estimate median effect size of surface sterilization
DIVC$CONTROL <- factor(DIVC$CONTROL)
library(brms)
gamma_m2 <- brm(bf(gamma ~ CONTROL + (1 | sample_id)),
                data = DIVC,
                iter = 8000,
                warmup = 4000,
                cores = 4)

summary(gamma_m2)

# preparing data for plotting
DIVCb <-
  DIVC %>%
  dplyr::select(sample_id, CONTROL, gamma) %>%
  tidyr::spread(key = 'CONTROL', value = 'gamma', fill = NA) %>%
  tidyr::gather(key = 'CONTROL', value = 'gamma', -sample_id) %>%
  tidyr::spread(key = 'CONTROL', value = 'gamma') %>%
  dplyr::mutate(diff = `1` - `0`)

pp_diff <- data.frame(coef = as.matrix(gamma_m2)[,2])

# generate plots for supplement:
# panel (a): paired gamma
ss_p1 <- ggplot(DIVC, aes(x = CONTROL, y = gamma, group = sample_id)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlab('surface sterilized') +
  ylab(expression(gamma)) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        plot.margin = unit(c(1,0,1,1),'lines')) +
  scale_y_continuous(limits = c(-10,5), breaks = seq(-10,5,2.5)) +
  scale_x_discrete(labels = c('yes','no'))

# panel (b): distribution of differences over-plotted with marginal effect
ss_p2 <- ggplot(DIVCb, aes(x = 1, y = diff)) +
  geom_jitter(width = 0.05, alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(1,0,0,0),'lines')) +
  xlab('') +
  ylab(expression(Delta ~ gamma)) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  scale_y_continuous(limits = c(-2, 5), breaks = seq(-2,5,1)) +
  scale_x_continuous(limits = c(0.8,1.2))

ss_p3 <- ggplot(pp_diff, aes(x = coef)) +
  geom_density(fill = 'gray40', col = NA, alpha = 0.5) +
  theme_minimal() +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0,1,0,-1),'lines')) +
  xlab('') +
  ylab('') +
  scale_x_continuous(limits = c(-2, 5), breaks = seq(-2,5,1)) +
  geom_vline(xintercept = quantile(pp_diff$coef, probs = c(0.025, 0.975)), lty = 3, col = 'gray40') +
  geom_vline(xintercept = quantile(pp_diff$coef, probs = c(0.5)), lty = 1, col = 'gray40')

ss_plot_all <- ggpubr::ggarrange(plotlist = list(ss_p1, ss_p2, ss_p3), ncol = 3,
                  align = 'h',
                  widths = c(1,0.4,0.25), labels = c('a','b',''))

ggsave(ss_plot_all, file = file.path(outdir,'ss_plot_all.pdf'),
       width = 4,
       height = 2.5,
       device = 'pdf',
       units = 'in')

###
### APPROACH 2: manually inspect cp and mt sequences
###

# export .fasta of cp and mt sequences
cp_seq_name <- paste0('>seq',rownames(TAX2)[cp])
cp_seqs <- rownames(TAX2)[cp]

mt_seq_name <- paste0('>seq',rownames(TAX2)[mt])
mt_seqs <- rownames(TAX2)[mt]

# create interleaved character vectors of sequences and names:
MT <- character(length(mt_seqs) * 2)
MT[c(TRUE, FALSE)] <- mt_seq_name
MT[c(FALSE, TRUE)] <- mt_seqs
writeLines(MT, con = paste0(data.path,"mtASV_seqs_DIVC.fasta"), sep = "\n")

CP <- character(length(cp_seqs) * 2)
CP[c(TRUE, FALSE)] <- cp_seq_name
CP[c(FALSE, TRUE)] <- cp_seqs
writeLines(CP, con = paste0(data.path,"cpASV_seqs_DIVC.fasta"), sep = "\n")

library(taxize)

# import blastn hit tables of results
mtbr <- read.csv(paste0(data.path,"mt_hit_table_DIVC.csv"), F)
names(mtbr) <- c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")

cpbr <- read.csv(paste0(data.path,"cp_hit_table_DIVC.csv"), F)
names(cpbr) <- c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")

write.hit.num <- function(df){
  df[,'hit.num'] <- NA
  uniques <- unique(df[,'seqid'])
  for(s in 1:length(uniques)){
    len <- length(df[,'seqid'][df[,'seqid'] == paste0(uniques[s])])
    df[,'hit.num'][df[,'seqid'] == paste0(uniques[s])] <- c(1:len)
  }
  return(df)
}

# write hit number in order for later ranking taxonomic matches
mtbr <- write.hit.num(mtbr)
cpbr <- write.hit.num(cpbr)

# define function to assign taxonomy to the GI numbers of the hits
get_taxonomy <- function(x) {
  #ENTREZ_KEY="90cc1824a010404743ee8240935b44464207"
  paste0(genbank2uid(x, key = ENTREZ_KEY)[[1]][1]) # taxonomy ID
}

# now get taxID from NCBI using taxize function genbank2uid
# grab unique gi numbers:
mtbr_gi <- data.frame(gi = unique(mtbr[,'subject']), tax_ID = NA) # n = 593
cpbr_gi <- data.frame(gi = unique(cpbr[,'subject']), tax_ID = NA) # n = 243

mtbr_gi[,'tax_ID'] <- sapply(mtbr_gi[,'gi'], get_taxonomy) # this is slow, since function makes n queries to ncbi. ~3 min.
cpbr_gi[,'tax_ID'] <- sapply(cpbr_gi[,'gi'], get_taxonomy) # this is slow, since function makes n queries to ncbi. ~2 min.

mt_class.all <- classification(unique(mtbr_gi[,'tax_ID']), callopts = list(), return_id = TRUE, db = 'ncbi')
#cpbr_gi2 <- cpbr_gi[cpbr_gi[,'tax_ID']!="NA",]
cp_class.all <- classification(unique(cpbr_gi[,'tax_ID']), callopts = list(), return_id = TRUE, db = 'ncbi')

# turn list into data.frame
mt_df.all <- do.call(rbind, mt_class.all)
cp_df.all <- do.call(rbind, cp_class.all)

# merge to match GI number with tax_ID, and then grab taxonomic name that matches tax_ID
# mtbr.all  <- merge(mtbr,mtbr_gi, by = 'gi', sort = F)
mtbr[,'tax_ID'] <- sapply(mtbr[,'subject'], function(x) mtbr_gi[,'tax_ID'][match(x, mtbr_gi[,'gi'])])
mtbr[,'name']   <- sapply(mtbr[,'tax_ID'], function(x) mt_df.all[,'name'][match(x, mt_df.all[,'id'])])

cpbr[,'tax_ID'] <- sapply(cpbr[,'subject'], function(x) cpbr_gi[,'tax_ID'][match(x, cpbr_gi[,'gi'])])
cpbr[,'name']   <- sapply(cpbr[,'tax_ID'], function(x) cp_df.all[,'name'][match(x, cp_df.all[,'id'])])


# ------------------------------------------ #
# analyzing effects of surface sterilization #
# ------------------------------------------ #

# plotting library size differences
lib_size_diffs <- ggplot(DIVC, aes(x = factor(CONTROL), y = log(row.sum,10), group = sample_id)) +
  geom_point() +
  geom_line()


a1 <- lmerTest::lmer(log(row.sum,10) ~ CONTROL + (1 | sample_id), data = DIVC)
summary(a1)


# OK so first finding: sterilization decreases number of reads. That's fairly clear from these data.
# Next I need to determine gamma by identifying cp + mt sequences in this set.
# I'll basically copy down the routine from the post-processing script
# and re-apply it here to these samples.

#
