# run_fruit_set_analysis.R

#header
library(here)
source(here("scripts/phy_functions.R"))
source(here("scripts/phy_header.R"))

# load NP data
npb1 <- read.table(here("data/112213.npd.txt"),T,"\t")

# remove redundant plot ID
npb1 <- dplyr::filter(npb1, plot != '58b')
npb1 <- droplevels(npb1)

# set reference level of condition to MOCK
npb1$condition <- relevel(npb1$condition, ref = 'MOCK')

# aggregate by subplot (i.e., patch)
npb2 <- dplyr::group_by(npb1, plot, subplot, condition, tx) %>% summarise(n.mines.sum = sum(n.mines,na.rm=T),
                                                                          n.leaves.sum = sum(n.leaves,na.rm=T),
                                                                          height.avg = mean(height, na.rm =T),
                                                                          n.mined.leaves.sum = sum(n.mined.leaves),
                                                                          n.yellow.leaves.sum = sum(n.yellow.leaves),
                                                                          n.fruits.sum = sum(n.fruits,na.rm=T),
                                                                          mean.stem.height = mean(height, na.rm = T))

# run Bayesian model here:
fruits.bm1 <- brm(bf(n.fruits.sum ~ mean.stem.height + n.mined.leaves.sum + (1|plot)),
                  family = negbinomial(),
                  data = npb2)

fruits.bm2 <- brm(bf(n.fruits.sum ~ mean.stem.height + n.leaves.sum + n.mined.leaves.sum + (1|plot)),
                  family = negbinomial(),
                  data = npb2)

# save model for future reference
saveRDS(fruits.bm1, file = here("models/fruits.bm1.rds"))
saveRDS(fruits.bm2, file = here("models/fruits.bm2.rds"))


#fruits.bm1 <- readRDS(file = here("models/fruits.bm1.rds"))

# load model if not run
# fruits.bm1 <- readRDS(file = here("models/fruits.bm1.rds"))

# generate marginal_effects plot from model fit
me1 <- marginal_effects(fruits.bm1, effects = 'n.mined.leaves.sum', plot = FALSE, points = TRUE)
me2 <- marginal_effects(fruits.bm2, effects = 'n.mined.leaves.sum', plot = FALSE, points = TRUE)

# now re-construct plot with points:
seed_plot1 <- ggplot() +
  geom_line(data = me2[[1]], aes(x = n.mined.leaves.sum, y = estimate__), col = "gray40", lwd = 1) +
  geom_line(data = me2[[1]], aes(x = n.mined.leaves.sum, y = lower__), col = "gray80") +
  geom_line(data = me2[[1]], aes(x = n.mined.leaves.sum, y = upper__), col = "gray80") +
  geom_point(data = npb2, aes(x = n.mined.leaves.sum, y = n.fruits.sum), col = "gray20", alpha = 0.5) +
  scale_y_continuous(limits = c(0,255), breaks = seq(0,250,50)) +
  scale_x_continuous(limits = c(0,155), breaks = seq(0,150,25)) +
  xlab("mined leaves per patch") +
  ylab("fruits set per patch")

ggsave(filename = here("figs/new_seeds_1.pdf"), width = 2.5, height = 2)
saveRDS(seed_plot1, file = here("figs/seed_plot1.rds"))

# now produce what will become the marginal density plot:
seed_hist1 <- ggplot(npb2, aes(x = n.fruits.sum)) + geom_density(col = "gray60") + xlab("") + ylab("") +
  scale_x_continuous(limits = c(0,255), breaks = seq(0,250,50)) +
  # theme(axis.line = element_blank(),
  #       axis.text = element_blank()) +
  coord_flip()

saveRDS(seed_hist1, file = here("figs/seed_hist1.rds"))


## secondard analysis looking at hormone tx on relationship between fruit set and herbivore damage:
# plot relationship broken down by tx versus control:
# make CTR condition all one plot:
npb2$condition2 <- npb2$condition
npb2$tx[npb2$condition=='MOCK'] <- 'zm'
npb2$tx <- factor(npb2$tx, levels = c('zm','tx'))
seed_plot2 <- ggplot() +
  facet_wrap( ~ condition) +
  #geom_line(data = me1[[1]], aes(x = n.mined.leaves.sum, y = estimate__), col = "gray40", lwd = 1) +
  #geom_line(data = me1[[1]], aes(x = n.mined.leaves.sum, y = lower__), col = "gray80") +
  #geom_line(data = me1[[1]], aes(x = n.mined.leaves.sum, y = upper__), col = "gray80") +
  geom_point(data = npb2, aes(x = n.mined.leaves.sum, y = n.fruits.sum, col = tx), alpha = 0.5) +
  scale_y_continuous(limits = c(0,255), breaks = seq(0,250,50)) +
  scale_x_continuous(limits = c(0,155), breaks = seq(0,150,25)) +
  xlab("mined leaves per patch") +
  ylab("fruits set per patch") +
  scale_color_manual(values = c('gray40','dodgerblue'))
seed_plot2

# need to test for condition effect here...
fruits.bm3 <- brm(bf(n.fruits.sum ~ mean.stem.height + n.mined.leaves.sum + condition + (1|plot)),
                  family = negbinomial(),
                  data = npb2,
                  cores = 4,
                  iter = 8000)

summary(fruits.bm3)

LOO(fruits.bm2,fruits.bm1)


#### Output of model results in latex format ####
fr1.c <- tidy(fruits.bm1)

# multiple herbivore damage term by 10:
fr1.c[grep('*mined*',fr1.c$term),-1] <- fr1.c[grep('*mined*',fr1.c$term),-1] * 10

fr2 <- cbind(term = fr1.c[,1], round(fr1.c[,-c(1)], 2)) # round and create output data.frame
fr2$coef <- with(fr2, paste0(estimate,' [', lower, ';', upper, ']')) # collapse estimate into single string

fr3 <- fr2[,c('term','coef')]
fr3 <- rbind(fr3[grep('^b_',fr3$term),],
             fr3[grep('^sd_',fr3$term),],
             fr3[grep('^shape',fr3$term),])

# output:
fr3$term <- c('Intercept','height (cm)','mined leaves (x 10)','sigma(plot)','NB shape')

# write lines to output:
con <- file(here("tables/FruitSet_model_table_v1.tex"))
open(con, 'wr')
writeLines(kable(fr3, "latex", caption = "Model estimates for patch-level fruit set as function of herbivore damage", booktabs = T, escape = TRUE) %>%
             kable_styling(position = "center"),
           con = con, sep = ''
)
close(con)

