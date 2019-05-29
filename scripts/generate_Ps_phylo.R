# generate_Ps_phylo.R
# last updated: 2018-OCT-25 by PTH

# Description:
# uses rentrez library to access, manipulate and generate phylogeny with Pseudomonadaceae full-length 16S sequences

library(here)
#install.packages("rentrez")
#devtools::install_github("renkun-ken/rlist")
library(rlist)
library(rentrez)
source(here("scripts/phy_header.R"))

# try to search gene database for 16S with search term matching Pseudomonadaceae
# filter results by those that are full-length. Let's see if this is possible:
Ps_16S <- entrez_search(db="nuccore", term="16S*[Gene] AND Pseudomonas[ORGN]", use_history=TRUE, retmax = 5700)

# now let's fetch some sequences, write them to file, then decide which to keep:
for(seq_start in seq(1,5500,50)){
  recs <- entrez_fetch(db="nuccore", web_history=Ps_16S$web_history,
                       rettype="fasta", retmax=50, retstart=seq_start)
  cat(recs, file="Ps_16S.fasta", append=TRUE)
  cat(seq_start+49, "sequences downloaded\r")
}

# load sequences back in; discard those that aren't terribly long
library(ape)
Ps_hits <- ape::read.dna(here('Ps_16S.fasta'), format = 'fasta')

# grab sequences hitting focal groups:
focal_clades <- c('fluorescens','syringae','viridiflava','putida','stutzeri')
Ps_hits2 <- Ps_hits[c(sample(grep(focal_clades[1],names(Ps_hits)),50,replace=FALSE),
                      grep(focal_clades[2],names(Ps_hits)),
                      grep(focal_clades[3],names(Ps_hits)),
                      sample(grep(focal_clades[4],names(Ps_hits)),20,replace=F))]

# now exclude all sequences that hit 23S, since this is NOT right locus:
Ps_hits2 <- Ps_hits2[c(1:length(Ps_hits2))[!c(1:length(Ps_hits2)) %in% c(grep('23S',names(Ps_hits2)))]]

# got 120 sequences; let's look at length distribution:
# filter all reads out with fewer than 1000 nucleotides:
Ps_len <- do.call(rbind,lapply(Ps_hits2, length))

hist(Ps_len[,1],breaks=100)
cut1 <- 500 # assign cutoff to near-full-length sequences
cut2 <- 1550

# grab indexes of all lengths passing cutoffs
cuts <- (Ps_len[,1]>cut1) & (Ps_len[,1]<cut2)
Ps_hits3 <- Ps_hits2[cuts] # 530 sequences retained, from focal groups. Nice!

# write these to file, add back in our focal bASVs, then re-load and align:
write.dna(Ps_hits3, format = 'fasta', file = here("Ps_hits3.fasta"))

# now, need to bring in actual bASV sequences
Ps_bASVs <- read.dna(here("blast-res/bASVs/Pseudomonadaceae_bASVs copy.fasta"), format = 'fasta')
write.dna(Ps_bASVs, format = 'fasta', file = here('Ps_hits3.fasta'), append = TRUE)

# run alignment using mafft (default parameters)
#system("mafft --adjustdirectionaccurately Ps_hits3.fasta > Ps.aln.mafft.fasta")
system("linsi --adjustdirectionaccurately Ps_hits3.fasta > Ps.aln.mafft.fasta") # slower, more accurate

# now, produce phylogeny of this alignment using phangorn:
library(phangorn)
Ps.aln <- read.phyDat(here("Ps.aln.mafft.fasta"), format = 'fasta')

# remove redundantly named sequences
# do by grabbing all non-uniques, removing the first one to match:
non_uniques <- names(table(names(Ps.aln))[table(names(Ps.aln))>1])
Ps.aln2 <- subset(Ps.aln, subset = !(names(Ps.aln) %in% names(Ps.aln)[match(non_uniques,names(Ps.aln))]))

# quick and dirty dist matrix and NJ tree as starter for ML inference
Ps.d1 <- dist.ml(Ps.aln2)
Ps.NJ <- NJ(Ps.d1)

# store output:
pdf(file = here("Ps_tree_nj1.pdf"), width = 10, height = 6)
  plot(Ps.NJ, cex = 0.05)
dev.off()

# now: generate tree:
fitStart = pml(Ps.NJ,as.phyDat(Ps.aln2),k=4)
fit = optim.pml(fitStart, model="GTR", optGamma=TRUE, rearrangement="stochastic")

# or, choose best model:
mt <- modelTest(as.phyDat(Ps.aln2), tree=Ps.NJ, multicore=TRUE)
bestmodel <- mt$Model[which.min(mt$AICc)]
env = attr(mt, "env")
fitStart = eval(get(bestmodel, env), env)
fit = optim.pml(fitStart, rearrangement = "stochastic",
                optGamma=TRUE, optInv=TRUE, model="GTR")
bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE) # finally, fit tree via bootstrapping

# save it after all that work!
saveRDS(fit, file = here("Ps.bs.tree.fit"))
saveRDS(bs, file = here("Ps.bs.tree"))

# tree_plot2 <- plotBS(midpoint(fit$tree), bs, p = 50, type="p", cex = 0.2)
# ggsave(tree_plot2, filename = here("figs/Ps_tree_bs.pdf"),width = 10, height = 8)

# need to trim off several of these taxa; something is wrong with these sequences.


