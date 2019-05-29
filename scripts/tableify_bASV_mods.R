## script to generate table-form output of focal models from Family-level bASV analysis.
## Also output, for each site, are model results for the allBact focal models.
## Last updated: 2018-Sep-24 by PTH

library(here)
source(here("scripts/phy_header.r"))
source(here("scripts/phy_functions.r"))

# load model list file
mods <- read.csv(here("models/best_mod_list_v1.txt"))

# re-order rows according to match with Fig 2:
FAMS <- paste0(read.csv(here("models/all_fams.csv"), header = F)[,1])
FAMS <- FAMS[!FAMS %in% c('Rhizobiales_Incertae_Sedis','Enterococcaceae')]
fam_tax <- unique(bTAX[,c('Phylum','Class','Order','Family')]) %>% dplyr::filter(Family %in% FAMS) %>% arrange(desc(Phylum), desc(Class), desc(Order), desc(Family))
fam_order <- c('all-bact',paste0(fam_tax$Family))
mods$Family <- factor(mods$Family, levels = fam_order)
mods <- dplyr::arrange(mods, dataset, Family)

# function to run through models and grab relevant coefficients
mod_coef_extractor <- function(x){
  if (x$dataset == 'EL'){
    if (x$Family == 'all-bact'){
      site_dir <- here("models/all-bacteria/")
    } else {
      site_dir <- here("models/EL/")
    }
  } else {
    if (x$Family == 'all-bact'){
      site_dir <- here("models/all-bacteria/")
    } else{
      site_dir <- here("models/NP/")
    }
  }
  if (x$Family != 'all-bact'){
    the_mod <- with(x,readRDS(grep(paste0(Family), grep(paste0(best_mod),Sys.glob(paste0(site_dir,'*')),value=T),value=T)))
  } else{
    the_mod <- with(x,readRDS(grep(paste0(dataset),grep(paste0(best_mod),Sys.glob(paste0(site_dir,'*')),value=T),value=T)))
  }

  t1 <- tidy(the_mod, prob = 0.95) # grab all terms
  t2 <- t1[grep('^b_|sigma',t1$term),] # grab population-level effects, including sigma terms

  if (x$best_mod %in% c('ga2','ga4','skn2','skn4')){
    # generate posterior samples
    ps1 <- posterior_samples(the_mod)
    # grab focal columns, exponentiate
    ps2 <- ps1[,grep('^b_sigma',names(ps1))] %>% exp()
    # summarise col by col
    ps3 <- t(apply(ps2,2,function(x) quantile(x,probs = c(0.025,0.975))))
    # take posterior means:
    ps_mu <- apply(ps2,2,mean)
    # now add std.error
    ps_se <- apply(ps2,2,function(x) sd(x))
    # put all together:
    psf <- data.frame(names(ps_mu),ps_mu,ps_se,ps3)
    names(psf) <- c('term','estimate','std.error','lower','upper')
    # merge with t2, removing previous versions:
    t2b <- t2[!t2$term %in% psf$term,]
    t2 <- rbind(t2b,psf)
  }

  # add additional effects like alpha and variance terms for intercept and slope of bASV:
  if (length(grep('skn',paste0(x$best_mod)))>0){
    t2 <- rbind(t2,t1[grep('alpha',t1$term),])
  }
  if (last(unlist(strsplit(paste0(x$best_mod),''))) %in% c('2','3','4')) {
    t2 <-  rbind(t2,t1[grep('^sd_',t1$term),])
  }

  t3 <- cbind(term = t2[,1], round(t2[,-c(1,3)], 2)) # round and create output data.frame
  t3$coef <- with(t3, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string
  t4 <- data.frame(site = x$dataset,
                   Family = x$Family,
                   model = paste0(x$best_mod),
                   term = t3$term,
                   coef = t3$coef)
  return(t4)
  #sprintf("Output: %s %s %s", x$dataset, x$Family, x$best_mod)
}

# # run this via lapply on the data.frame given by mods
# # for debugging:
xx <- split(mods, list(mods$dataset, mods$Family))
x <- xx[[11]]

mod_res <- lapply(split(mods, list(mods$dataset, mods$Family)), mod_coef_extractor)
mod_res2 <- do.call(rbind, mod_res)

# debug with loop to figure out where it breaks:
# xx <- split(mods, list(mods$dataset, mods$Family))
# xxx <- list()
# for (i in seq_along(xx)){
#   xxx[[i]] <- mod_coef_extractor(xx[[i]])
# }

# need to make sure model factors are all the same:
mod_res2$term <- sapply(mod_res2$term, function(x) if (length(grep('herb_dmg',paste0(x))>0)) { gsub('herb_dmg1','herb_dmg',x)} else { paste0(x) })

# change 'sigma' to 'b_sigma_Intercept' so it collapses with single sigma value in final table:
mod_res2$term[mod_res2$term=='sigma'] <- sapply(mod_res2$term[mod_res2$term=='sigma'], function(x)  gsub('sigma','b_sigma_Intercept',x))

# generate final factor
mod_res2$term <- factor(mod_res2$term, levels = paste0(unique(mod_res2$term)))

mod_res3 <- reshape2::dcast(mod_res2, site + Family + model ~ term, value.var = 'coef')
mod_res3[is.na(mod_res3)] <- paste0('--')

# modify the res for printing:
mod_res4 <- mod_res3[,-1]

# re-order:
mod_res4 <- mod_res4[,c(1,2,3,4,9,10,7,5,6,8)]
#names(mod_res4) <- c("Family","model","$\\alpha_{0}$","$\\beta_{0}$","$s_{\\alpha}$","$s_{\\beta}$","$\\sigma_{0}$","$\\sigma_{\\beta}$","$\\sigma_{\\gamma}$")
names(mod_res4) <- c("Family","model","$\\alpha_{0}$","$\\beta_{0}$","$\\tau_{\\alpha}$","$\\tau_{\\beta}$","$\\tau_{p}$","$\\sigma_{0}$","$\\sigma_{1}$","$\\alpha$")

# now produce tables as .tex documents to import into working draft:
#install.packages("kableExtra")
# build table call:
con <- file(here("tables/Family-level_and_all-bact_model_table_v3.tex"))
open(con, 'wr')
writeLines(kable(mod_res4, "latex", caption = "Model results for Family-level bacterial abundance", booktabs = T, escape = FALSE) %>%
             kable_styling(latex_options = c("scale_down"), position = "center") %>%
             group_rows("site EL", 1, 15) %>%
             group_rows("site NP", 16, 30), con = con, sep = ''
)
close(con)

