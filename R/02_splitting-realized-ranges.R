## script to split realized ranges into equatorward and poleward components and then rasterize them
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
library(rmapshaper)
select <- dplyr::select
## devtools::install_github("r-spatial/lwgeom")
library(lwgeom)


####################################################################################
#####                       SPLITTING RANGE SHAPEFILES                        ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in shape files
realized_ranges <- st_read()


## loop through all realized ranges
i = 1
while (i < length(st_geometry(realized_ranges)) + 1) {
  
  print(paste("On species number:", i, sep = " "))
  
  ## get range and latitudinal midpoint
  ## for ranges that are very complex and take forevvvvver to split, simplify first:
  if (i %in% c(100, 173, 190, 191, 192, 196, 197, 198, 199, 205, 206, 209, 225, 229, 426, 524)) {
    range <- ms_simplify(realized_ranges[i, ], keep_shapes = TRUE)
  }
  else {
    range <- realized_ranges[i, ]
  }
  
  extent <- extent(range)
  latitudinal_midpoint <- (as.numeric(ymin(extent)) + as.numeric(ymax(extent))) / 2
  lmp_ls <- st_linestring(rbind(c(-180, latitudinal_midpoint), c(180, latitudinal_midpoint)))
  
  ## for northern and southern species, chop the range at the latitudinal midpoint:
  if (range$hemisphere != "EQUATOR") {
    cut_line <- lmp_ls
    lat_of_cut <- latitudinal_midpoint
    print(paste("Species range is in the", range$hemisphere, "hemisphere.", sep = " "))
  }
  ## for equator crossers, chop the range at the equator:
  else if (range$hemisphere == "EQUATOR") {
    cut_line <- equator
    lat_of_cut <- 0
    print("Species range crosses the equator.")
  }
  
  ## split the range: 
  split <- st_split(range, cut_line)
  
  ## separate into polygons that are below cut vs polygons that are above cut line:
  polys <- list(st_geometry(split))[1][[1]][[1]]
  above <- which(unlist(lapply(polys, 
                               function (x) st_coordinates(st_centroid(x))[2]))
                 > lat_of_cut) %>%
    polys[.]
  
  below <- which(unlist(lapply(polys, function (x) st_coordinates(st_centroid(x))[2])) 
                 < lat_of_cut) %>%
    polys[.]
  
  ## create simple features to house the split pieces of the ranges:
  sf_below <- st_sfc(lapply(below, function (x) st_polygon(x))) %>%
    st_as_sf(.) %>%
    mutate(species = range$species) %>%
    mutate(latitudinal_midpoint = latitudinal_midpoint) %>%
    mutate(poleward_or_equatorward = "equatorward") %>%
    mutate(hemisphere = range$hemisphere) %>%
    mutate(source = range$source)
  
  sf_above <- st_sfc(lapply(above, function (x) st_polygon(x))) %>%
    st_as_sf(.) %>%
    mutate(species = range$species) %>%
    mutate(latitudinal_midpoint = latitudinal_midpoint) %>%
    mutate(poleward_or_equatorward = "poleward") %>%
    mutate(hemisphere = range$hemisphere) %>%
    mutate(source = range$source)
  
  if (i == 1) {
    sf_cumulative <- rbind(sf_above, sf_below)
  }
  else {
    sf_cumulative <- rbind(sf_cumulative, sf_above, sf_below)
  }
  
  i = i + 1
}

## combine polygons of each species equatorward and poleward range together into a multipolygon:
sf_cumulative <- sf_cumulative %>%
  mutate(species_type = paste(species, poleward_or_equatorward)) %>%
  group_by(species_type) %>%
  aggregate(., list(.$species_type), function(x) x[1]) %>%
  select(-species_type, -Group.1)

##saveRDS(sf_cumulative, "data-processed/sf_cumulative.rds")
rds <- readRDS("data-processed/sf_cumulative.rds")
st_write(sf_cumulative, "data-processed/realized-ranges_split.shp")



## code to visualize each range-splitting action:
countries <- ne_countries(returnclass = "sf") %>%
  st_transform(., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(st_geometry(range), col = "lightpink") ## range boundary
plot(st_geometry(polys[[1]]), add = TRUE, col = "lightblue") ## plot polygons separately
plot(st_geometry(polys[[2]]), add = TRUE, col = "lightyellow")
plot(st_geometry(sf_above), add = TRUE, col = "orange2") ## plots all polygons above midpoint
plot(st_geometry(sf_below), add = TRUE, col = "yellow1") ## plots all polygons below midpoint
plot(st_geometry(countries), add = TRUE) 
plot(st_geometry(equator), add = TRUE, col = "red") ## equator
plot(st_geometry(lmp_ls), add = TRUE, col = "orange") ## latitudinal midpoint