### Insect herbivory alters the structure of leaf microbiomes
### Humphrey & Whiteman
### last updated: PTH 15-JUN-2018
## Function definition file
## phy_functions.R


logit <- function(x){
  return(log(x /(1-x)))
}


#### POSTERIOR PREDICTIVE CHECK PLOTTING FUNCTIONS ####
pp.max <- function(sim1, dat1, fact1 = 'herb_dmg', response1 = 'log_ratio'){

  sim1 <- data.frame(sim1)

  # join data and sim data:
  dat2 <- cbind(dat1, sim1)

  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]

  # calculate maximum value of all simulations:
  max.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, max)
  max.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, max)

  # calculate max value of observed zeros by grouping factor:
  lvl1.max <- max(lvl1[,paste0(response1)])
  lvl2.max <- max(lvl2[,paste0(response1)])

  # construct data.frame
  res1 <- data.frame(max.val = c(as.vector(max.sim.1),as.vector(max.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(max.val = c(lvl1.max,lvl2.max),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))

  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.max < max.sim.1))/length(max.sim.1),3)
  q.2 <- round(length(which(lvl2.max < max.sim.2))/length(max.sim.2),3)

  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)

  # plot:
  plot1 <- ggplot() +
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = max.val), bins = 25, col = "gray40", fill = "gray40") +
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = max.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") +
    xlab(expression(max(Y))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_phy1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  return(list('plot' = plot1, 'data' = res2))
}



pp.stat <- function(sim1, dat1, fact1, response1, FUN = max, xlabs = "max"){

  sim1 <- data.frame(sim1)

  # join data and sim data:
  dat2 <- cbind(dat1, sim1)

  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]

  # calculate maximum value of all simulations:
  stat.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, FUN)
  stat.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, FUN)

  # calculate max value of observed zeros by grouping factor:
  lvl1.stat <- FUN(lvl1[,paste0(response1)])
  lvl2.stat <- FUN(lvl2[,paste0(response1)])

  # construct data.frame
  res1 <- data.frame(stat.val = c(as.vector(stat.sim.1),as.vector(stat.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(stat.val = c(lvl1.stat,lvl2.stat),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))

  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(stat.sim.1 >= lvl1.stat))/length(stat.sim.1),3)
  q.2 <- round(length(which(stat.sim.2 >= lvl2.stat))/length(stat.sim.2),3)

  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)

  # plot:
  plot1 <- ggplot() +
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = stat.val), bins = 25, col = "gray40", fill = "gray40") +
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = stat.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") +
    xlab(xlabs) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_phy1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  return(list('plot' = plot1, 'data' = res2))
}

#### functions for posterior predictions for log ratio and CFU models ####

grab_rds <- function(...){
  fam_files <- Sys.glob(paste0(IN_DIR, FAMS[m],'*'))
  fam_mods <- list()
  for (f in 1:length(fam_files)){
    fam_mods[[f]] <- readRDS(file = fam_files[f])
  }
  return(fam_mods)
}

grab_mod_names <- function(x){
  #gsub('.rds','',unlist(strsplit(gsub(IN_DIR,'',x),'_'))[3])
  gsub('.rds','',grep('.rds',unlist(strsplit(gsub(IN_DIR,'',x),'_')), value = T))
}

bernoulli_trials<-function(p)
{
  U <- runif(length(p), 0, 1) # draw p numbers from uniform distribution
  outcomes <- U < p # determine if draw is less than p, taken as a vector (i.e., element-wise comparison between elements of U and p). If the uniform number is less than p, score as success.
  return(as.numeric(outcomes)) # return binary vector of same length as p.
}

sim_bernoulli <- function(a0,b0,a1,b1,...){
  # generate vector of probabilities to use in bernoulli trials
  prob_vec <- replicate(R, c(rbeta(n = h0, shape1 = a0, shape2 = b0), rbeta(n = h1, shape1 = a1, shape2 = b1)))
  # now conduct bernoulli trials, one for each element of prob_vec
  berns <- bernoulli_trials(prob_vec)
  return(berns)
}

# draw log_ratio posterior samples from best_mod
leaf_level_pplr <- function(x, the_mod, ...){
  y_rep <- brms::posterior_predict(the_mod, newdata = x, nsamples = R)
  return(y_rep)
}

# function def to generate predicted CFUs from draws from joint posterior of cfu_mod; returns predicted CFU on linear (i.e. count) scale.
# PROBABLY OBSOLETE!!!
cfu_post_pred.old <- function(x, IDs, log=FALSE, ...){
  l <- length(unique(IDs))
  jp <- data.frame(cfu_mod)
  cfu_jpd <- jp[sample(1:length(jp[,1]),l,replace=T),] # draw from joint posterior for each row of the
  cfu_jpd2 <- cfu_jpd[IDs,] # replicate each draw from joint posterior across all sub-items for each leaf sample
  cfu <- cfu_jpd[,1] + cfu_jpd[,2] * x + rnorm(n = 1, mean = 0, sd = cfu_jpd[,3])

  ifelse(log==TRUE,return(10^cfu),return(cfu))
}

# adding pp_cfu using tapply:
# PROBABLY OBSOLETE!!!
cfu_post_pred_lowest <- function(x, log=FALSE, ...){
  #l <- length(unique(IDs))
  jp <- data.frame(cfu_mod)
  cfu_jpd <- jp[sample(1:length(jp[,1]),1,replace=T),] # draw from joint posterior
  deviates <- sapply(cfu_jpd[,3], function(s) rnorm(n = length(x), mean = 0, sd = s))
  cfu <- cfu_jpd[,1] + cfu_jpd[,2] * x + deviates

  ifelse(log==FALSE,return(10^cfu),return(cfu))
}

# function to generate joint posterior draws per samp_id for each rep, with observation-level normal deviates per sigma draw per sample
# the output matrix will be of the same dimension as input matrix
pp_cfu_per_leaf <- function(x, samp_col = 'samp_id', ...){
  jp <- data.frame(cfu_mod)
  dat1 <- expand.grid(samp_id = c(1:max(as.numeric(x[,samp_col]))),
                      rep = c(1:max(as.numeric(x$rep))))
  n <- dim(dat1)[1]
  cfu_jpd <- cbind(dat1, jp[sample(1:length(jp[,1]),n,replace=T),]) # draw from joint posterior
  the_rows <- rep(c(1:dim(cfu_jpd)[1]), each = length(unique(x$bASV))) # add to original data.frame
  x2 <- cbind(x, cfu_jpd[the_rows,-c(1,2)])
  x2$deviate <- sapply(x2$sigma, function(x) rnorm(n = 1, mean = 0, sd = x)) # draw normal deviates for each row:
  x2$log_pp_cfu <- x2$b_Intercept + x2$b_ln_bac_host * x2$ln_bac_host + x2$deviate # now calculate pp_cfu from model coefficients
  x2$pp_cfu <- 10^x2$log_pp_cfu # exponentiate y_rep values
  return(x2)
}


pp_cfu_per_leaf_all <- function(x, samp_col = 'samp_id', the_fact = 'sp_id', ...){
  jp <- data.frame(cfu_mod)

  if (is.numeric(x[,samp_col])==FALSE){
    dat1 <- expand.grid(samp_id = unique(x[,samp_col]),
                        rep = c(1:max(as.numeric(x$rep))))
  } else {
    dat1 <- expand.grid(samp_id = c(1:max(as.numeric(x[,samp_col]))),
                        rep = c(1:max(as.numeric(x$rep))))
  }
  n <- dim(dat1)[1]
  cfu_jpd <- cbind(dat1, jp[sample(1:length(jp[,1]),n,replace=T),]) # draw from joint posterior
  the_rows <- rep(c(1:dim(cfu_jpd)[1]), each = length(unique(x[,the_fact]))) # add to original data.frame
  x2 <- cbind(x, cfu_jpd[the_rows,-c(1,2)])
  x2$deviate <- sapply(x2$sigma, function(x) rnorm(n = 1, mean = 0, sd = x)) # draw normal deviates for each row:
  x2$log_pp_cfu <- x2$b_Intercept + x2$b_ln_bac_host * x2$ln_bac_host + x2$deviate # now calculate pp_cfu from model coefficients
  x2$pp_cfu <- 10^x2$log_pp_cfu # exponentiate y_rep values
  return(x2)
}


#### ACCESSORY FUNCTIONS ####
row_rep <- function(x, fact = 'n.leaves', fact2 = 'n.mined.leaves'){
  colz <- names(x)[!names(x) %in% fact]
  x2 <- x[rep(row.names(x), x[,fact]), colz]
  x2$herb_dmg <- 0
  if (as.vector(x[,fact2])>0){
    x2$herb_dmg[row.names(x2) %in% sample(row.names(x2),x[,fact2])] <- 1
  }
  x2$leaf_id <- seq(1,x[,fact],1)
  x2$plant_id <- paste0(x2$subplot,'_',
                        x2$plant)
  x2$sample_id <- paste0(x2$plant_id,'_',
                         x2$leaf_id)
  return(x2)
}

#### DIVERSITY FUNCTIONS ####

shannon <- function(x, is.prop=FALSE){
  if (is.prop == TRUE){
    tmp1 <- x[x>0]
    H <- -sum(tmp1*log(tmp1,2))
  } else {
    tmp1 <- x/sum(x)
    tmp2 <- tmp1[x>0]
    H <- -sum(tmp2*log(tmp2,2))
  }
  return(H)
}

evenness <- function(x, is.prop=FALSE){
  if (is.prop == TRUE){
    tmp1 <- x[x>0]
    H <- -sum(tmp1*log(tmp1,2))
  } else {
    tmp1 <- x/sum(x)
    tmp2 <- tmp1[x>0]
    H <- -sum(tmp2*log(tmp2,2))
    Hmax <- log(length(tmp2),2)
  }
  return(H/Hmax)
}

SJ.divergence <- function(x,y,base=2,sq=F){
  # remove double 0 entries
  xx <- x[!((x==0) & (y==0))]
  yy <- y[!((x==0) & (y==0))]

  if(sum(xx)>1){
    # convert to prop
    xx <- xx/sum(xx)
    yy <- yy/sum(yy)
  }

  # calculate average probability vector
  m <- 0.5 * (xx + yy)

  # calculate vectors of xx/m and yy/m and set zero values to 1
  xm <- xx/m
  xm[xm==0] <- 1

  ym <- yy/m
  ym[ym==0] <- 1

  SJ <- 0.5 * (sum(xx * log(xm,base)) + sum(yy * log(ym,base)))
  if(sq==F){return(SJ)} else{return(sqrt(SJ))}
}

# function to calculate all diversity indexes on two-row data.frame input x (supplied via lapply)
calc_div <- function(x){
  s0 <- shannon(x[1,-c(1:2)])
  s1 <- shannon(x[2,-c(1:2)])
  e0 <- evenness(x[1,-c(1:2)])
  e1 <- evenness(x[2,-c(1:2)])
  sj <- SJ.divergence(x[1,-c(1:2)],x[2,-c(1:2)])
  data.frame(rep = x[1,1], H0 = s0, H1 = s1, Hd = s1 - s0, E0 = e0, E1 = e1, Ed = e1 - e0, SJ = sj)
}

