# R script to export all sequences from bASVs from each Family in set of families passed from file
# inputs: (1) sequences + taxonomy file for all bASVs
#         (2) FAMS files
# Last updated: PTH 21-JUL-18

library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

bTAX <- read.csv(here("data/bTAX_table_26-JUN-2018.csv"))
FAMS <- paste0(read.csv(here("models/all_fams.csv"),header = F)$V1)
OUTDIR <- here("blast-res/bASVs/")

NPD1 <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/NPD_long_v1.csv")
ELD1 <- read.csv("/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/coinfection/data/ELD_long_v1.csv")

# exclude non-zero rows so we discard any bASVs that happen to have no hits
NPD1 <- dplyr::filter(NPD1, bASV_count>0)
ELD1 <- dplyr::filter(ELD1, bASV_count>0)

# export list of bASVs to use across all datasets
bASVs <- paste0(unique(unique(ELD1$bASV[ELD1$Family %in% FAMS]),
                       unique(NPD1$bASV[NPD1$Family %in% FAMS])))

# run through each Family and export .fasta of bASVs in that family
# ensuring that the bASV is actually in the list of bASVs that we will consider (bASVs)
bTAX2 <- dplyr::filter(bTAX,unique.id %in% bASVs)

# now run through Families and export .fasta of all bASVs each
for (f in seq_along(FAMS)){
  # sub-set seq-tax data for only focal Family
  fam_tmp <- dplyr::filter(bTAX2, Family == FAMS[f])

  # now export sequences as .fasta, for all bASVs:
  cat(paste0(">",fam_tmp$unique.id,"\n",fam_tmp$X), file = paste0(OUTDIR,FAMS[f],"_bASVs.fasta"), sep = "\n")

  # export single sequence per Family for Family-level phylogeny:
  fs <- sample(seq(1,length(fam_tmp[,1])),1)
  if(f == 1){
    cat(paste0(">",fam_tmp$Family[fs],"\n",fam_tmp$X[fs]), file = paste0(OUTDIR,"Family-seqs.fasta"), sep = "\n")
  } else{
    cat(paste0(">",fam_tmp$Family[fs],"\n",fam_tmp$X[fs]), file = paste0(OUTDIR,"Family-seqs.fasta"), sep = "\n", append = TRUE)
  }
}

