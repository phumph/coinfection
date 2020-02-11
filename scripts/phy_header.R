### Insect herbivory reshapes a native leaf microbiome
### Humphrey & Whiteman
### last updated: PTH 15-JUN-2018

## Header file for all R scripts

## phy_header.R

# install.packages("here")
# install.packages("devtools")
# install.packages("ggplot2")
## devtools::install_github('hadley/ggplot2', build_vignettes = TRUE)
## devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
## install.packages("shiny")
## devtools::install_github("ropensci/plotly")
# install.packages("ggExtra")
# install.packages("ggpubr")
# install.packages("brms", type="source")
# install.packages('BH')
# install.packages('RcppEigen')
# install.packages('Rcpp')
# install.packages("broom")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("viridis")
# devtools::install_github("marcosci/cividis")
# devtools::install_github("ropensci/taxize") # need to install version >0.9.2
# install.packages("phangorn")
# install.packages("Matrix")
# install.packages("glmmTMB")
# install.packages("bbmle")
# install.packages("kableExtra", type="source")

## setup ggplot themes:
theme_ph1 <- function (base_size = 11, base_family = "sans")
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      #panel.border = element_rect(fill = NA, colour = "black"),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = 0.25),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_rect(fill = "white", colour = NA),
      complete = TRUE)
}

theme_phy1 <- function(base_size = 9, base_family = "sans") {
  theme_bw(
    base_family = base_family,
    base_size = base_size
  ) +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.4),
      #axis.ticks = element_line(size = 0.3),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = rel(0.9)),
      strip.placement = "outside",
      # strip.background = element_rect(fill = "gray95", color = NA),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.text = element_text(size = 9),
      legend.text.align = 0,
      legend.key = element_blank()
    )
}

## load libraries
library(ggExtra)
library(ggpubr)
library(dplyr)
library(reshape2)
library(brms)
library(Rcpp)
library(RcppEigen)
library(rlang)
library(ggplot2); theme_set(theme_phy1())
library(viridis)
library(cividis)
library(phangorn)
library(kableExtra)
library(broom)
