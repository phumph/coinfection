# script to produce table for CFU model
# last updated 2018-SEP-20 by PTH

library(here)
source(here("scripts/phy_header.r"))
source(here("scripts/phy_functions.r"))

# load CFU mod
cfu_mod <- readRDS(here("models/cfu_m1.rds"))
cc1 <- tidy(cfu_mod, prob = 0.95)[1:3,]
cc2 <- cbind(term = cc1[,1], round(cc1[,-c(1)], 2)) # round and create output data.frame
cc2$coef <- with(cc2, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string
cc3 <- cc2[,c('term','coef')]
cc3$term <- c('Intercept','Slope','$\\sigma_{\\varepsilon}$')
names(cc3) <- c('Term','Estimate')
con <- file(here("tables/CFU_model_table_v2.tex"))
open(con, 'wr')
writeLines(kable(cc3, "latex", caption = "Model coefficients linking 16S counts to CFU", booktabs = T, escape = FALSE) %>%
             kable_styling(position = "center"), con = con, sep = ''
)
close(con)

