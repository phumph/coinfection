# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# if (length(args)<5) {
#   stop("Positional arguments are: [1] input dir where models are stored; [2] output dir where reports will be put; [3] Rmarkdown template file that generates the report; [4] prevalence/abundance stats file' [5] OTU table", call.=FALSE)
# }

## R script to automate report generation from brmsfit object, stored as .rds object.

#IN_DIR <- args[1]
#IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/'
IN_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/NPB_v1/'
 
#OUT_DIR <- args[2]
#OUT_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/out_v1/'
OUT_DIR <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/models/NPB_out_v1/'
 
#REP_TEMPLATE <- args[3]
REP_TEMPLATE <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/scripts/OTU_rep_template.Rmd'
 
#STATS_FILE <- args[4]
STATS_FILE <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/abund_stats_all_sites.txt'
 
#OTU_TAB <- args[5]
#OTU_TAB <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/EL_otu_v5.txt'
OTU_TAB <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/data/NPB_otu_v5.txt'

# load required libraries
library(knitr)
library(markdown)
library(rmarkdown)

# call up directory and grob all model files 
mod_fil <- rev(Sys.glob(paste0(IN_DIR, '*.rds')))

# grab summary prevalence statistics for relevant OTUs
otu_stats <- read.table(paste0(STATS_FILE), T, '\t')
otu_tab <- read.table(paste0(OTU_TAB), T, '\t')

## start main report loop:
for (mod in 1:length(mod_fil)){
#for (mod in 1:2){  
  print(paste("Working on file", mod_fil[mod]))
  OTU <- readRDS(paste0(mod_fil[mod]))
  
  # grab short filename for output
  otu_name <- sub('.rds','',sub(pattern = paste0(IN_DIR), '', mod_fil[mod]))
  
  # for pdf reports
  rmarkdown::render(input = paste0(REP_TEMPLATE),
    output_format = "pdf_document",
    output_file = paste0("modrep_",otu_name,"_",Sys.Date(), ".pdf"),
    output_dir = paste0(OUT_DIR))
}