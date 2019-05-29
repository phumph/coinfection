# tableify_bASV-v-hormones_mods.R
# this script takes input directory of model results and outputs coefficient estimates
# for specified models.
# last updated: 2019-APR-04 by PTH

# header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))

# define names of models to grab:
the_mods <- c('h1d','h1d2','h1u','h1u2')
INDIR <- here("models/bact_v_hormones/NP/")

# load up models into list, where each list element has each model:
mod_files <- lapply(as.list(the_mods), function(x) Sys.glob(paste0(INDIR,'*',x,'.rds')))

# define functions to load/extract coefficients from focal models:
parse_models <- function(x,add_loo=FALSE){
  # define dependency function
  coef_extractor <- function(y,...){
    # load model:
    MOD <- readRDS(y)
    # extract coefficients and compile results into useful format:
    # optional: add LOO-IC information and re-write:
    if (add_loo==TRUE){
      MOD <- add_ic(MOD, ic = "loo", reloo = TRUE)
      saveRDS(MOD,file=y)
    }

    # extract outer interval from
    coef1 <- tidy(MOD,prob=0.95,par_type = "non-varying")
    coef1[,-1] <- round(coef1[,-1],3)

    # add model name, strip from filename:
    strsplit(y,'_') %>% unlist() %>% last() -> modstr
    modstr <- gsub('.rds','',modstr)
    coef1$model <- modstr

    # gather inner 50 percentile of coefficient posterior
    coef2 <- tidy(MOD,prob=0.5,par_type = "non-varying")
    coef2[,-1] <- round(coef2[,-1],3)

    coefs_comb <- cbind(coef1[grep('stem_',coef1$term),c('term','lower')],
                        coef2[grep('stem_',coef2$term),c('lower','estimate','upper')],
                        coef1[grep('stem_',coef1$term),c('upper','model')]
    )
    taxon <- strsplit(y,'/') %>% unlist()
    coefs_comb$taxon <- strsplit(grep('.rds',taxon,value=T),'_') %>% unlist() %>% first()
    names(coefs_comb) <- c('term','q0.025','q0.25','estimate','q0.75','q0.975','model','taxon')
    rm(MOD)
    return(coefs_comb)
  }

  # take input vector of model files as input
  # make into list
  # lapply function to load, extract coefficients, remove model:
  all_coefs <- lapply(as.list(x),coef_extractor)
  all_coefs2 <- do.call(rbind,all_coefs)

return(all_coefs2)
}

# perform routine to gather all coefs:
all_mod_coefs <- lapply(mod_files, function(x) parse_models(x, add_loo=TRUE))

# assemble all into one df:
all_mod_coefs2 <- do.call(rbind,all_mod_coefs)

# filter out bad taxa:
all_mod_coefs2 <- all_mod_coefs2[!all_mod_coefs2$taxon %in% c('Paenibacillaceae','Rhizobiaceae','Streptococcaceae') ,]

# generate plot:
ggplot(all_mod_coefs2[all_mod_coefs2$model %in% c('h1d2','h1u2'),], aes(x = taxon, y = estimate)) +
  geom_linerange(aes(ymin = q0.025, ymax = q0.975, x = taxon), col = 'gray40', size = 0.5) +
  geom_linerange(aes(ymin = q0.25, ymax = q0.75, x = taxon), col = 'gray40', size = 1) +
  geom_point(col="black") +
  coord_flip() +
  geom_hline(yintercept = 0) +
  facet_grid(model ~ term)

# try to add LOO and re-save objects. Then I can conduct the loo comparison.
all_mod_coefs3 <- lapply(mod_files, function(x) parse_models(x, add_loo=FALSE))

# compare loo for each model object again:
# use mod_files to load all models for a given taxon and compare the models.
the_fams <- unique(all_mod_coefs[[1]]$taxon)
# remove certain families with too few data points:
the_fams <- the_fams[!the_fams %in% c('Paenibacillaceae','Rhizobiaceae','Streptococcaceae')]

# load models
loo_parser <- function(mod_names=c('h0d','h1d','h0d2','h1d2','h0u','h1u','h0u2','h1u2'),...){
  mod_files2 <- lapply(as.list(mod_names), function(x) Sys.glob(paste0(INDIR,'*',x,'.rds')))
  # load models
  for (f in seq_along(the_fams[-c(1:2)])){
    the_fam_mod_files <- lapply(mod_files2,function(x) grep(the_fams[f],x,value=T))
    the_fam_mods <- lapply(the_fam_mod_files,readRDS)
    # extract loo and store:
    the_loos <- lapply(the_fam_mods,loo,reloo=T)
    names(the_loos) <- mod_names
    the_loos2 <- do.call(rbind, sapply(the_loos, function(x) with(x, data.frame(model_name, looic, se_looic)), simplify = F)) %>% dplyr::arrange(looic)
    the_loos2$model_name <- mod_names
    the_loos2$taxon <- the_fams[f]

    # export loo values
    write.csv(the_loos2, file = paste0(here("models/bact_v_hormones/NP/"),"LOO-IC_",the_fams[f],"_bASVs-v-hormones_v1.csv"),quote=F,row.names=F)

    # compare LOO-IC values
    # loo_comps <- brms::compare_ic(the_fam_mods[[1]],
    #                               the_fam_mods[[2]],
    #                               the_fam_mods[[3]],
    #                               the_fam_mods[[4]]
    # )
  }
}

loo_parser()
