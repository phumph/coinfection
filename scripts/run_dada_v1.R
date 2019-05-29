
# Rscript for performing DADA2 pipeline on EL and NPB sample sets
# PTH
# Nov 15 2017
# inputs:
  # --dir, fastq file input directory [==DIR]
  # --filt, output directory for filtered files (defaults to DIR) [==FILT]
  # --QC, output directory for QC information (defaults to DIR) [==QC_OUT]
  # --out, output directory for summary statistics and ASV table output (defaults to DIR) [==OUT]
  # --params, parameter input file for dada [==PARAMS]:
    # maxEE = maximum number of expected errors
  
suppressPackageStartupMessages(library(optparse))
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
suppressPackageStartupMessages(library("devtools"))
devtools::install_github("benjjneb/dada2")
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(bayesplot))


# need to setup argument parsing
#### PARSING COMMAND LINE ARGUMENTS ####
option_list <- list(
  make_option(c("-i", "--dir"), action="store", type="character", default=NA, 
              help="Input file of processed counts (output of counts_parser.R"), 
  
  make_option(c("-o", "--output"), action="store", type="character", default=NA, 
              help="Output directory (full path; must exist) for figures and data"), 
  
  # make_option(c("-e", "--exclude"), action="store", type="character", default=NA,
  # help="Time point(s) to exclude, comma- or semi-colon-separated"),
  
  make_option(c("-e", "--exclude"), action="store", type="integer", default=NA,
              help="Time point to exclude from fitness calculation (only single values supported for now..)"),
  
  make_option(c("-c", "--cutoff"), action="store", type="integer", default=20,
              help="Define cutoff for initial time-point read count for inclusion [default %default]"),
  
  make_option(c("-r", "--iterations"), action="store", type="integer", default=3,
              help="Define the number of iterations to run of the model [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))




## temporary variable assignments for debugging:
QC_PATH <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/16S_seq/demult/dada2_v1/QC/"
DIR     <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/16S_seq/demult/out"
FILT    <- "/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/16S_seq/filt/"

runQC <- TRUE
setwd(paste0(DIR))

# collect .fastq files
fnFs <- sort(list.files(DIR, pattern="-R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(DIR, pattern="-R2.fastq", full.names = TRUE))

# Get sample names
sample.names <- sapply(strsplit(fnFs, "/"), tail, n=1) # remove directory information from sample names
sample.names <- sapply(sample.names, function(x) paste(unlist(strsplit(x, "-"))[1]))

if (runQC == TRUE){
  # export quality profiles of N=n samples:
  n <- 50
  ns <- sample(length(fnFs),n)
  the_F_files <- fnFs[ns]
  the_R_files <- fnRs[ns]
  the_names <- sample.names[ns]
  
  for(fil in 1:length(the_F_files)){
    QP1 <- plotQualityProfile(paste0(path,the_F_files[fil]))
    QP2 <- plotQualityProfile(paste0(path,the_R_files[fil]))
    png(filename = paste0(outpath,"QProf.",the_names[fil], ".png"), width = 640, height = 960)
      grid.arrange(QP1,QP2,nrow=2)
    dev.off()
  }
} else if (doFilt == TRUE) {
  
  filtpath <- FILT
  filtFs <- paste0(filtpath, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq")
  filtRs <- paste0(filtpath, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq")
  
  # need to supply full path:
  out <- filterAndTrim(fnFs[1:4], filtFs[1:4], fnRs[1:4], filtRs[1:4], truncLen=c(150,150),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=FALSE, multithread=FALSE)

  
} else{
  # just do DADA on filt files, with error objects re-learned.
  
}
# took about 3 min to run on 50 samples.

  head(out)

outpath <- '/Users/phumph/Dropbox/Phyllosphere_project/analysis_phy/16S_seq/demult/dada2_v1/QC/'

# calibrate the base-wise error model for fwd and rev reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# plot error profiles and parametric fit to data
EP1 <- plotErrors(errF, nominalQ=TRUE)
EP2 <- plotErrors(errR, nominalQ=TRUE)

# save file to QC outpath
png(filename = paste0(outpath,"EProf.png"), width = 640, height = 960)
grid.arrange(EP1,EP2,nrow=2)
dev.off()
```

Dereplicate the filtered fastq files:
```{r message=FALSE}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

Run dada on the forward and reverse samples:
```{r}
dadaFs <- dada(derepFs, err=errF, pool=TRUE, selfConsist=FALSE, multithread=TRUE) # 60m
dadaRs <- dada(derepRs, err=errR, pool=TRUE, selfConsist=FALSE, multithread=TRUE) # 40m
```

Merge:
```{r}
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
```

Make sequence table:
```{r}
seqtab.all <- makeSequenceTable(mergers)
dim(seqtab.all) # 360 samples, 1936 sequences (with chimeras)
```

Identify and remove chimeras:
```{r}
bim <- isBimeraDenovo(seqtab.all, verbose=TRUE)
seqtab <- seqtab.all[,!bim]
```

Assign taxonomy:
```{r}
taxa <- assignTaxonomy(seqtab.nochim, "Training/silva_nr_v128_train_set.fa.gz", multithread=TRUE) # if lots of NAs, may need to revcomp the refdb. Use tryRC=TRUE
taxa <- addSpecies(taxa, "Training/silva_species_assignment_v128.fa.gz")

# export for the inspection
write.csv(taxa, file="dada2_taxa_Silva_wSpp_v1.csv", quote=F)
```

