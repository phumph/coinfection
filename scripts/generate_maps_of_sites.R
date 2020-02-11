#!/usr/bin/env Rscript

# generate_maps_of_sites.R
# last updated 2018-OCT-24 by PTH

# tips from http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html

source("./phy_header.R")

# load mapping library
#install.packages(c("maps", "mapdata"))
#devtools::install_github("dkahle/ggmap", ref = "tidyup")

if (!requireNamespace("maps", quietly = TRUE))
    install.packages("maps")

if (!requireNamespace("mapdata", quietly = TRUE))
    install.packages("mapdata")

if (!requireNamespace("ggmap", quietly = TRUE))
    devtools::install_github("dkahle/ggmap", ref = "tidyup")

library(ggmap)

# load GPS data for EL:
gps_el <- read.table(file.path("../data/gps_EL2012.txt"),T,'\t')
gps_np <- read.table(file.path("../data/gps_NP2013.txt"),T,'\t')

# load Google Maps API key:
GAPI <- 'AIzaSyDtrCgJGbpxjttAntTnl81YX3oj4qCCn6g' # git yer own!!!
register_google(key = GAPI)

# load 'em up
library(ggmap)
library(maps)
library(mapdata)

# start by identifying base map of Colorado:
states <- map_data("state")
western_us <- dplyr::filter(states, region %in% c('california','idaho','oregon','washington','nevada','arizona','utah','colorado','new mexico','wyoming','montana'))
western_us_map <- ggplot(data = western_us) + theme_nothing() +
  geom_polygon(aes(x = long, y = lat, group = group), color = "gray40", fill = "gray90") +
  coord_fixed(1.3)

# try to make map of each site
all_gps <- rbind(gps_el,gps_np)

all_bbox <- make_bbox(lat = lat, lon = long, data = all_gps)
lon_offset <- 0.035
lat_offset <- 0.035

all_bbox[[1]] <- all_bbox[[1]]-lon_offset
all_bbox[[3]] <- all_bbox[[3]]+lon_offset

all_bbox[[2]] <- all_bbox[[2]]-lat_offset
all_bbox[[4]] <- all_bbox[[4]]+lat_offset

sites   <- get_map(location = all_bbox, source = "google", maptype = "terrain")

# make map:
region_map <- ggmap(sites) +
  geom_point(data = all_gps, mapping = aes(x = long, y = lat), col = 'darkorange2')

# now make NP map:
el_bbox <- make_bbox(lat = lat, lon = long, data = gps_el)
np_bbox <- make_bbox(lat = lat, lon = long, data = gps_np)

# generate small offsets:
offset_bbox <- function(bbox, offset){
  bbox[[1]] <- bbox[[1]]-offset
  bbox[[3]] <- bbox[[3]]+offset
  bbox[[2]] <- bbox[[2]]-offset
  bbox[[4]] <- bbox[[4]]+offset
  return(bbox)
}

el_bbox2 <- offset_bbox(el_bbox, offset = 0.0005)
np_bbox2 <- offset_bbox(np_bbox, offset = 0.003)

# grab maps:
el_site <- get_map(location = el_bbox2, source = "google", maptype = "satellite")
np_site <- get_map(location = np_bbox2, source = "google", maptype = "satellite")

# now make maps:

# define elevation breaks:
b <- c(10000,10500,11000,11500)

ELM <- ggmap(el_site) +
  geom_point(data = gps_el[sample(c(1:40),replace=F),], mapping = aes(x = long, y = lat, fill = elevation), pch = 21, col = 'black') +
  scale_fill_viridis(limits = c(10400,11900), alpha = 0.5, begin = 0.4,
                     #colours=c("navyblue", "darkmagenta", "darkorange1"),
                     breaks=b, labels=format(b))

NPM <- ggmap(np_site) +
  #geom_point(data = gps_np, mapping = aes(x = long, y = lat), pch = 21, fill = 'darkorange2', col = 'black')
  geom_point(data = gps_np, mapping = aes(x = long, y = lat, fill = elevation), pch = 21, col = 'black') +
  scale_fill_viridis(limits = c(10400,11900), alpha = 0.5, begin = 0.4,
                       #colours=c("navyblue", "darkmagenta", "darkorange1"),
                       breaks=b, labels=format(b))

MAPbottom <- ggarrange(plotlist = list(ELM,NPM), ncol = 2, align = 'hv', common.legend = T)

ggsave(MAPbottom, filename = file.path("../figs/map_bottom.pdf"), width = 6.5, height = 3)

MAPtop <- ggarrange(plotlist = list(western_us_map, region_map), ncol = 2, align = 'hv')
ggsave(MAPtop, filename = file.path("../figs/map_top.pdf"), width = 6.5, height = 3)

# end
