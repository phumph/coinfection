### collapse_OTUs.R
### PTH 10 Oct 2017
# purpose: collapse MED NODES that are probably pseudo-counts of 16S copies from the same genome
# inputs:
  # 1. OTU table
  # 2. fasta file with 16S sequences from all MED NODES
  # 3. taxonomy file that groups MED NODES into family-level taxonomic classification

# WHAT THIS SCRIPT DOES:
# calculates pairwise distance between resulting sequences (raw nucleotide differences, as a proportion).
# takes user-imputted OTU table and calculates point-estimates of Pearson (and Spearman) rank correlation between NODES in the OTU table
# collapses NODES in OTU table according to nucleotide distance % cutoff (≤1% for this case)
# generates lists of OTUs that were collapsed into each new OTU

# outputs:
  # 1. collapsed OTU table
  # 2. list of all NODES that were collapsed
  # 2. family-level fasta alignments are retained in the pwd.

#### step 0. ####
# load libraries
require(reshape2)
require(dplyr)
require(picante)
require(ape)
require(ggplot2)
require(gridExtra)

#### step 1: ####
# load fasta alignment of representative sequences for each MED NODE
ntt <- read.csv("../data/phy/node_to_tax.txt",header=F,colClasses='character') # MED NODE-to-taxonomy mapping file
tax <- read.table("../data/phy/TAXONOMY_HIERARCHY.txt",T,"\t") # taxonomic hierarchy of MED NODES mapped to taxonomy
seq <- read.FASTA("../data/phy/node_reps_aln.filter.fasta")

# use taxonomy mapper to assign taxon names to MED NODE IDs
names(ntt) <- c("NODE","TAX")
seq_name_map <- merge(data.frame(NODE = names(seq)),ntt, by = "NODE", sort = FALSE)
# cbind(seq_name_map, names(seq)) # OK they match; re-assign names
names(seq) <- seq_name_map[,'TAX']

#### step 2: ####
# first, for EL
# load OTU table data and compute Spearman and Pearson correlation among samples for all OTUs
otu_tab <- read.table("../data/EL_v4.tsv",T,"\t")

# extract matching columns in names(seq)
otu_tab <- otu_tab[,names(otu_tab) %in% names(seq)] # we don't care if we lose sample IDs because this is about OTU-wise correlations

# remove all OTUs with total counts of zero
removes <- names(otu_tab[,colSums(otu_tab)==0])
otu_tab <- otu_tab[,!(names(otu_tab) %in% removes)]


# second, from NPB
otu_tab_NPB <- read.table("../data/NPB_v4.tsv",T,"\t")
# extract matching columns in names(seq)
otu_tab_NPB <- otu_tab_NPB[,names(otu_tab_NPB) %in% names(seq)] # we don't care if we lose sample IDs because this is about OTU-wise correlations

# remove all OTUs with total counts of zero
removes2 <- names(otu_tab_NPB[,colSums(otu_tab_NPB)==0])
otu_tab_NPB <- otu_tab_NPB[,!(names(otu_tab_NPB) %in% removes)]

# remove sequences for which we have no OTU information
# to_prune_from_seq <- names(seq)[names(seq) %in% names(otu_tab)==FALSE]
# seq <- seq[!(names(seq) %in% to_prune_from_seq)]
# 
# 
# # generate specific file for NPB
# to_prune_from_seq2 <- names(seq)[names(seq) %in% names(otu_tab_NPB)==FALSE]
# seq2 <- seq[!(names(seq) %in% to_prune_from_seq)]

#### step 3: ####
# calculate distances
dist1 <- dist.dna(seq, model = "raw") # calculates raw distance between sequences. No evolutionary model implied.
# define depth at which to cut the tree. 0.01 corresponds to 1% (raw) sequence divergence.
# note: we don't care that more distant sequences will have inflated divergence owing to poor alignment of these sequences.
# this is because the whole goal here is to distinguish closely-releated sequences, which will be aligned properly with one another,
# even if poorly aligned with other properly-aligned subsets.
cuts_0.01 <- cutree(hclust(dist1, method = "complete"), h = 0.009)

#### step 4: OUTPUT DISTANCE STATISTICS ####
# NOTE: done for EL; not done for NPB!!!
# calculate pairwise correlations among OTUs BEFORE binning into 1% OTUs
library(picante)
cp1 <- cor.table(log(1+otu_tab), cor.method = "pearson")$r
cs1 <- cor.table(log(1+otu_tab), cor.method = "spearman")$r

# grab only lower triangle
# functions from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# # Get upper triangle of the correlation matrix
# get_upper_tri <- function(cormat){
#   cormat[lower.tri(cormat)]<- NA
#   return(cormat)
# }

cp2 <- get_lower_tri(cp1)
cs2 <- get_lower_tri(cs1)

# melt correlation tables and merge with dist.dna object 'dist1' above
cp3 <- reshape2::melt(cp2)
cs3 <- reshape2::melt(cs2)
dist2 <- get_lower_tri(as.matrix(dist1))
dist3 <- reshape2::melt(dist2)

# merge distances with correlations for outputting
dist3[,'contr'] <- paste0(dist3[,'Var1'],'_',dist3[,'Var2'])
cp3[,'contr'] <- paste0(cp3[,'Var1'],'_',cp3[,'Var2'])
cs3[,'contr'] <- paste0(cs3[,'Var1'],'_',cs3[,'Var2'])
cordat1 <- merge(cp3, cs3, by = 'contr', sort = FALSE)
cordat2 <- merge(cordat1, dist3, by = 'contr', sort = FALSE)
cordat3 <- cordat2[,c(1,2,3,4,7,10)]
names(cordat3) <- c('contr','contr.x','contr.y','pearson','spearman','dna_dist')

# only examine correlations across small distances (≤3%):
close_rels <- dplyr::filter(cordat3, dna_dist <= 0.03, dna_dist > 0.0000) # grab all ≤ 3% sequences

# make plots as output
cplot1 <- ggplot(data = close_rels, aes(x = dna_dist, y = spearman)) + geom_point() + theme_bw() + ggtitle("spearman")
cplot2 <- ggplot(data = close_rels, aes(x = dna_dist, y = pearson)) + geom_point() + theme_bw() + ggtitle("pearson")
cplot3 <- ggplot(data = close_rels, aes(x = pearson, y = spearman, color = dna_dist)) + geom_point() + theme_bw() + theme(legend.position = "bottom") +
  scale_color_continuous(low = "black", high = "dodgerblue")

pdf(file = "../figs_working/pearson_spearman_dna_dist.pdf",width = 7, height = 5)
  grid.arrange(cplot1,cplot2,cplot3) # save this plot..
dev.off()

# now plot distributions of Pearson coefficients for each level of dna_dist..
# first increase phylogenetic depth of comparison by linking up higher-level taxonomy information
close_rels2 <- merge(close_rels, tax[,c(6,8)], by.x = 'contr.y', by.y = 'TAX_NAME')

# define binary distance factor

close_rels2[,'nt_dist'] <- c("1bp")
close_rels2[,'nt_dist'][close_rels2[,'dna_dist']>0.005] <- c("2bp")

# write.csv(file = "../data/close_rels_Pseudo_Sphingo_1.csv", dplyr::filter(close_rels2, dna_dist < 0.01, K %in% c('Pseudomonadaceae', 'Sphingomonadaceae')), quote = F)

distr_plot1 <- ggplot(data = dplyr::filter(close_rels2, dna_dist < 0.01, K %in% c('Pseudomonadaceae', 'Sphingomonadaceae')), aes(x = pearson)) + theme_bw() +
  geom_histogram(bins = 40, fill = "gray50") + facet_grid(factor(nt_dist) ~ K)

pdf(file = "../figs_working/pearson_dist_pseudo_sphingo.pdf", width = 5, height = 4)
  print(distr_plot1)
dev.off()

# examining whether total count of the OTU has any bearing on this search for respectable criterion..
# add colSum value to data.frame
the_col_sums <- data.frame(count = sort(colSums(otu_tab)))
the_col_sums[,'IDs'] <- row.names(the_col_sums)

close_rels2 <- merge(close_rels, the_col_sums, by.x = 'contr.x', by.y = 'IDs', sort = FALSE)
close_rels3 <- merge(close_rels2, the_col_sums, by.x = 'contr.y', by.y = 'IDs', sort = FALSE)

# sort by dna_dist
write.table(file = "../data/phy/OTU_correlations.txt", dplyr::arrange(close_rels3, dna_dist), sep = "\t")

# look at Pseudomonas only for example
pseudo <- close_rels3[grep("Pseudomonas", close_rels3[,'contr.y']),]
pplot1 <- ggplot(data = pseudo, aes(x = pearson, y = spearman, color = dna_dist, size = count.y)) + geom_point() + theme_bw() + theme(legend.position = "bottom") +
  scale_color_gradientn(colors = c("black", "red", "dodgerblue", midpoint = 0.01))

sphingo <- close_rels3[grep("Sphingomonas", close_rels3[,'contr.y']),]
sphplot1 <- ggplot(data = sphingo, aes(x = pearson, y = spearman, color = dna_dist, size = count.y)) + geom_point() + theme_bw() + theme(legend.position = "bottom") +
  scale_color_gradientn(colors = c("black", "red", "dodgerblue", midpoint = 0.01))

pdf(file = "pearson_cor_pseudo_sphingo.pdf", width = 8, height = 5)
  grid.arrange(pplot1, sphplot1, ncol = 2)
dev.off()

# now plot with only ≥99% relatives.

#### step 5: GENERATE 1% OTUS ####
# debugging
cuts <- cuts_0.01
orig_otu <- otu_tab_NPB

# write function to pool counts
pool_by_cluster <- function(cuts, orig_otu){
  k <- 1
  otu_map <- NULL
  
  # determine whether there are hits in 'cuts' that aren't represented in orig_otu
  cuts <- cuts[labels(cuts) %in% names(orig_otu)]
  cuts <- sort(cuts)
  ks <- unique(cuts)
  
  
  # define results data.frame
  res_tab <- data.frame(matrix(ncol = length(table(cuts)),
                               nrow = length(orig_otu[,1])))
  #for(k in 1:69){
  for(k in 1:length(unique(cuts))){
    # grab all names of first cut
    group_k_hits <- names(cuts[cuts %in% ks[k]])
    
    if (length(group_k_hits)==1){
      res_tab[,k] <- orig_otu[,names(orig_otu) %in% group_k_hits]
      names(res_tab)[k] <- group_k_hits
      otu_map <- c(otu_map, group_k_hits)
      next
    } else {
    
    # grab sub-table from orig_otu
    sub_tab <- orig_otu[,names(orig_otu) %in% group_k_hits]
    
    # grab abundances of each of these k hits from original otu_table
    k_sums <- sort(colSums(sub_tab), decreasing = T)
    sOTU_name <- names(k_sums[1])
    
    # create summed vector
    res_tab[,k] <- rowSums(sub_tab)
    names(res_tab)[k] <- sOTU_name
    otu_map <- c(otu_map, rep(sOTU_name,length(group_k_hits)))
    }
  }
  return(list('OTUs_new' = res_tab,
              'OTU_map' =  otu_map))
}

# generate pooled 1% table
otu_tab_EL_new  <- pool_by_cluster(cuts = cuts_0.01, orig_otu = otu_tab)
otu_tab_NPB_new <- pool_by_cluster(cuts = cuts_0.01, orig_otu = otu_tab_NPB)

# combine with previous OTU table meta-data
EL_old_tab  <- read.table("../data/EL_v4.tsv",T,"\t")
NPB_old_tab <- read.table("../data/NPB_v4.tsv",T,"\t")

meta_col_names <- c('samples','is.in_sample_prep_sheet', 'ROW_TOTALS','contam', 'total', 'prop_c', 'PlotNum', 'PlotCond', 'SubPlot', 'stem.tx','LeafPos', 'sample_mass', 'herb.dmg', 'LeafCFU', 'logLeafCFU')
EL_new_otu <- data.frame(EL_old_tab[,names(EL_old_tab) %in% meta_col_names], otu_tab_EL_new[[1]])

# compare rowSums between old and new OTU table (post-pooling)
EL_new_otu[,'ROW_TOTALS_1perc'] <- rowSums(otu_tab_EL_new[[1]])

# new totals are slightly lower owing to the fact that several OTUs in the prior table were not in the sequence file. That is expected--we don't know how to pool them.
# but this is an annoying discrepancy, and shouldn't happen if we have sequences for all OTUs...  


## now do the same for NPB
NPB_meta_col_names <- c('samples', 'ROW_TOTALS','contam', 'total', 'prop_c', 'PlotNum', 'PlotCond', 'SubPlot', 'stem.tx','sample_mass', 'herb.dmg')
NPB_new_otu <- data.frame(NPB_old_tab[,names(NPB_old_tab) %in% NPB_meta_col_names], otu_tab_NPB_new[[1]])
NPB_new_otu[,'ROW_TOTALS_1perc'] <- rowSums(otu_tab_NPB_new[[1]])

# plot(NPB_new_otu[,'ROW_TOTALS'],NPB_new_otu[,'ROW_TOTALS_1perc'])
# why are some samples losing so many sequenes? High abundance OTUs removed for some reason? I need to look into this more carefully..

# write tables: EL first
write.table(file = "../data/EL_otu_table_1percent.txt", otu_tab_EL_new[[1]], sep = '\t')
write.table(file = "../data/EL_otu_table_1percent_otu_map.txt", otu_tab_EL_new[[2]], sep = '\t')
write.table(file = "../data/EL_otu_v5.txt", EL_new_otu, sep = '\t')

# write tables: now NPB
write.table(file = "../data/NPB_otu_table_1percent.txt", otu_tab_NPB_new[[1]], sep = '\t')
write.table(file = "../data/NPB_otu_table_1percent_otu_map.txt", otu_tab_NPB_new[[2]], sep = '\t')
write.table(file = "../data/NPB_otu_v5.txt", NPB_new_otu, sep = '\t')
