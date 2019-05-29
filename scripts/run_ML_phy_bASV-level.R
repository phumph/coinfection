#!/usr/bin/env Rscript
# R script to automate phylogeny estimation (ML) from BLASTed sequences per Family
# Last updated: PTH 21-JUL-18

# load required libraries
library(knitr)
library(markdown)
library(rmarkdown)
library(here)
library(taxize)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# define blast-res source dir
BLAST_RES_DIR <- here("blast-res/bASVs/")
BLAST_SOURCE_DIR <- here("blast-res/bASVs/")
blast_res_filenames <- Sys.glob(paste0(BLAST_RES_DIR,"*.csv"))
OUTDIR <- here("blast-res/bASVs/out/")

# read files:
BFs <- list()
for (f in 1:length(blast_res_filenames)){
  BFs[[f]] <- read.csv(blast_res_filenames[f], header = F)
}

# grab column header with subject sequence:
seq_col <- "V4"
tax_col <- "V6"

# function to run classification on unique TAX_ID, then merge back to original list
grab_tax <- function(x,...){
  # re-name TAX_ID column by match to tax_col
  names(x)[grep(paste0(tax_col),names(x))] <- 'id'
  # grab taxonomy from NCBI for each unique TAX_ID in list element
  tax_uniques <- classification(unique(x$id), callopts = list(), return_id = TRUE, db = 'ncbi')
  # grab lowest rank and turn results into data.frame
  tax_uniques2 <- do.call(rbind,tax_uniques)
  tax_uniques2$query <- gsub('\\..+','',row.names(tax_uniques2)) # grab input lowest-level tax_ID (strip off redundant info from row.names)

  # toss out unimportant taxonomic levels
  tax_uniques2 <- dplyr::filter(tax_uniques2, !(rank %in% c('superkingdom','phylum','class','order','no rank')))

  # turn into single-row-per-query data.frame:
  tax_final <- reshape2::dcast(tax_uniques2, query ~ rank, value.var = 'name')

  # now merge back to original list
  return(merge(x, tax_final, by.x = 'id', by.y = 'query', sort = F, all.x = T))
}

# function to export fastas given column IDs from data.frame
## TO-DO: purge redundant sequences!!!
export_fasta <- function(x, seq_col, id_col, fam){
  # generate file-name for set of bASVs
  filename <- paste0(fam,"_hits.fasta")
  # grab only non-redundant sequences:
  x2 <- x[unique(x[,seq_col]),]
  x2[,id_col] <- gsub(' ','_',x2[,id_col]) # replace spaces with underscores for seq_names
  x2[,id_col] <- gsub('\\(|\\)|\\[|\\]','\\.',x2[,id_col]) # replace additional illegal characters (for RAxML) with underscores
  # now make all names non-redundant
  the_names <- NULL
  for(n in seq_along(x2[,id_col])){
    if (n == 1){
      the_names[n] <- x2[n,id_col]
      next
    } else {
      the_names[n] <- ifelse(!(x2[n,id_col] %in% the_names[1:(n-1)]),x2[n,id_col],paste0(x2[n,id_col],'_',length(grep(x2[n,id_col],the_names))+1))
    }
  }
  x2[,id_col] <- the_names
  cat(paste0(">",x2[,id_col],"\n",paste0(x2[,seq_col])), file = paste0(OUTDIR,filename), sep = "\n") # dump lines of sequence to file
}

# grab taxonomy and append to all list elements
BFs2 <- lapply(BFs, grab_tax, tax_col) # takes about 6 minutes on MacBook Pro 2012

# now export sequences for making phylogenies
fam_names <- gsub(BLAST_RES_DIR,'',blast_res_filenames)
fam_names <- gsub('_bASVs.fasta.csv','',fam_names)

for(i in seq_along(BFs2)){
  export_fasta(BFs2[[i]], seq_col, id_col = 'species', fam = fam_names[i])
}

# also spit out taxonomy guide..
# grab taxonomy columns:
export_tax_guide <- function(y,x){
  tax_col_names <- names(y)[!names(y) %in% names(x)]
  tax_guide <- unique(y[,tax_col_names]) %>% dplyr::arrange(V1)
  names(tax_guide)[1] <- 'unique_id'
  guide_filename <- paste0(unlist(strsplit(paste0(unique(y$family)[unique(y$family) %in% fam_names]), "_"))[1],"_tax-guide.csv")
  write.csv(tax_guide, file = paste0(OUTDIR,guide_filename), row.names = F)
}

# need to fix outfile.. not properly referencing family.

# STILL NEED TO DO THIS PART!
lapply(BFs2, export_tax_guide, BFs)

## now merge original queries with parsed results for alignment and tree-making.
# list of input files:
#
# blast_in_filenames <- Sys.glob(paste0(BLAST_RES_DIR,"*.fasta"))
#
