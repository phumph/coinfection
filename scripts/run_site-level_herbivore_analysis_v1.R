#! /usr/bin/Rscript
# Rscript to generate site-level analysis of herbivore patterns
# last updated 2018-SEP-18 by PTH

library(here)
source(here("scripts/phy_header.R"))
source(here("scripts/phy_functions.R"))

# load dataset to supply the newdata for predictive simulations
FOCAL.DAT.NAME <- here("data/EL.field.hormone.exp.master.txt")
focal.dat <- read.table(FOCAL.DAT.NAME,T,"\t")
names(focal.dat)[3] <- 'condition'

focal.dat$subplot <- with(focal.dat, paste0(plot,'_',tx))
focal.dat$plant <- with(focal.dat, paste0(sub.plot, stem.num))
focal.dat$condition <- relevel(focal.dat$condition, ref = 'CTR')


# do some exploratory plotting of aggregated herbivory levels across sub-plots:
subplot_data <- dplyr::group_by(focal.dat, subplot, condition, tx, plot) %>% summarise(n.mined.leaves.sum = sum(n.mined.leaves),
                                                                                       n.leaves.sum = sum(n.leaves),
                                                                                       height.avg = mean(height, na.rm = T))
subplot_data$group = paste0(subplot_data$subplot,'_',subplot_data$condition)
subplot_data$plot <- factor(subplot_data$plot)

ggplot(subplot_data, aes(x = tx, y = n.mined.leaves.sum, group = plot)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ condition)

# run model to estimate within-treatment effects:

el.ja.m1 <- brm(bf(n.mined.leaves ~ n.leaves + height + tx + (tx|plot) + (1|subplot)),
                  data = focal.dat[focal.dat$condition=='JA',],
                  family = negbinomial,
                  cores = 4)

# try model where counts are summed. Should get closer to normal approximation.
# this model seems over-parameterized.
el.ja.m2 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot) + (1|subplot)),
                data = subplot_data[subplot_data$condition=='JA',],
                family = gaussian,
                cores = 4)

el.ja.m3 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot)),
                data = subplot_data[subplot_data$condition=='JA',],
                family = gaussian,
                cores = 4,
                control = list(adapt_delta = 0.99))
plot(el.ja.m3)
el.ja.m4 <- brm(bf(n.mined.leaves.sum ~ n.leaves.sum + height.avg + tx + (tx|plot)),
                data = subplot_data[subplot_data$condition=='CTR',],
                family = gaussian,
                cores = 4,
                control = list(adapt_delta = 0.99))

# ggplot(subplot_data, aes(x = n.leaves.sum, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = n.mined.leaves.sum, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = n.leaves.sum, y = n.mined.leaves.sum, col = condition)) + geom_point()
#
#
# ggplot(subplot_data, aes(x = condition, y = n.leaves.sum, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = condition, y = height.avg, col = condition)) + geom_point()
# ggplot(subplot_data, aes(x = condition, y = n.mined.leaves.sum, col = condition)) + geom_point()


