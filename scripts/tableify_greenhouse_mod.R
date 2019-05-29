## script to generate table-form output of greenhouse analysis/
## Last updated: 2018-Sep-24 by PTH

library(here)
source(here("scripts/phy_header.r"))
source(here("scripts/phy_functions.r"))

# script to produce table for CFU model
# last updated 2018-SEP-20 by PTH

# load focal greenhouse model
#gh.m1 <- readRDS(file = here("models/gh_bm1.rds"))
gh.m1 <- readRDS(file = here("models/gh_bm1_v3.rds"))

gh1.c <- tidy(gh.m1, prob=0.95)
gh2 <- cbind(term = gh1.c[,1], round(gh1.c[,-c(1)], 2)) # round and create output data.frame
gh2$coef <- with(gh2, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string

gh3 <- gh2[,c('term','coef')]
gh3 <- gh3[grep('^b_',gh3$term),]
gh3$stat <- 'intercept'
gh3$stat[grep('plant.tx2$',gh3$term)] <- 'slope'
gh3$strain.id <- sapply(gh3$term, function(x) gsub(':plant.tx2$','',gsub('^b_strain.id','',x)))

# merge back with original gh data to find clade:
gh1 <- read.csv(file = here("data/greenhouse_doublings_data.csv"))
strain_info <- unique(gh1[,names(gh1) %in% c('strain.id','spp')])
gh4 <- merge(gh3[,names(gh3) %in% c('strain.id','stat','spp','coef')], strain_info, by = 'strain.id')
gh4 <- reshape2::dcast(gh4[,c('strain.id','stat','spp','coef')], spp + strain.id ~ stat, value.var = 'coef')

# re-level based on order from plot:
strain_order <- c("20A","22B","26B","46B","02A","20B","36A","29A","39A","33E","03A","46A")
gh4$strain.id <- factor(gh4$strain.id, levels = rev(strain_order))
gh4 <- dplyr::arrange(gh4, desc(strain.id))

con <- file(here("tables/GH_model_table_v4.tex"))
open(con, 'wr')
writeLines(kable(gh4[,-1], "latex", caption = "Strain-level model terms for experimental infection", booktabs = T, escape = TRUE) %>%
             kable_styling(position = "center") %>%
             group_rows("P. syringae", 1, 6) %>%
             group_rows("P. fluorescens", 7, 12),
           con = con, sep = ''
)
close(con)


