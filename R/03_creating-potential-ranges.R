## script to create potential range shapefiles for each species based on their thermal tolerance limits
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(ncdf4)
library(rnaturalearth)
library(smoothr)
library(lwgeom)
library(rmapshaper)
library(RColorBrewer)
select <- dplyr::select



####################################################################################
#####                       READ IN TEMPERATURE RASTERS                       ######
####################################################################################
## read temperature data created in script 00
raster_terr_low <- stack("./data-processed/raster_terr_low.grd")
raster_terr_high <- stack("./data-processed/raster_terr_high.grd")
raster_marine_low <- stack("./data-processed/raster_marine_low.grd")
raster_marine_high <- stack("./data-processed/raster_marine_high.grd")
raster_intertidal_low <- stack("./data-processed/raster_intertidal_low.grd")
raster_intertidal_high <- stack("./data-processed/raster_intertidal_high.grd")

## read in realm mask layers:
t_mask <- raster("./data-processed/raster_terr_mask.grd")
m_mask <- raster("./data-processed/raster_marine_mask.grd")
i_mask <- raster("./data-processed/raster_intertidal_mask.grd")

####################################################################################
#####                   RASTERIZE REALIZED RANGES                             ######
####################################################################################
## read in realized ranges:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
  mutate(range_id = paste(species, source, sep = "_")) %>%
  filter(source != 'GARD')


## constrain realized range rasters by habitat allowed to be in the potential range
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

i = 1
while (i < nrow(realized_ranges)+1) {
  ## for ranges that are very complex and take forevvvvver to rasterize, simplify first:
  if (i %in% c(100, 173, 190, 191, 192, 196, 197, 198, 199, 205, 206, 209, 225, 229, 426, 524)) {
    range <- ms_simplify(realized_ranges[i, ], keep_shapes = TRUE)
  }
  else {
    range <- realized_ranges[i, ]
  }
  rr_raster <- rasterize(range, r, background = NA, getCover = TRUE)
  rr_raster[rr_raster > 0] <- 1 ## set all values that overlap realized range in raster to 1
  rr_raster[rr_raster != 1] <- NA ## set all other values to NA
  
  ## constrain by habitat:
  if(range$realm == "Terrestrial" | range$realm == "Freshwater") {
    rr_raster <- mask(rr_raster, t_mask) ## set cells in mask that do not overlap realized range to NA
  }
  else if(range$realm == "Marine") {
    rr_raster <- mask(rr_raster, m_mask) 
  }
  else {
    rr_raster <- mask(rr_raster, i_mask) 
  }
  
  ## add to list of rasters
  if (i == 1) {
    rasterized_rrs <- rr_raster
  }
  else {
    rasterized_rrs <- addLayer(rasterized_rrs, rr_raster, updatevalue = NA)
  }
  
  print(paste("Finished range number:", i))
  i = i + 1
}

names(rasterized_rrs) <- realized_ranges$range_id

##saveRDS(rasterized_rrs, "data-processed/rasterized_rrs.rds")
rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")


####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))

## read in species traits:
traits <- read.csv("./data-processed/wrangled-traits.csv")

## subset lims to ones with traits for now:
thermal_limits <- thermal_limits[which(thermal_limits$genus_species %in% traits$genus_species),]

#######################################################
#####   TERRESTRIAL AND FRESHWATER SPECIES:      ######
#######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial" | realm == "Freshwater")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_terr_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPEICES WITH BOTH THERMAL LIMITS:
####################################=
combined_reg <- filter_by_tolerance(both_upper = both_upper, 
                                    both_lower = both_lower, 
                                    raster_high = raster_terr_high, 
                                    raster_low = raster_terr_low)

combined_dormancy <- filter_by_tolerance_dormancy(both_upper = both_upper, 
                                                       both_lower = both_lower, 
                                                       raster_high = raster_terr_high, 
                                                       raster_low = raster_terr_low)

combined_elev <- filter_by_tolerance_elev(both_upper = both_upper, 
                                          both_lower = both_lower, 
                                          raster_high = raster_terr_high, 
                                          raster_low = raster_terr_low)

combined_elev_x_dormancy <- filter_by_tolerance_elev_x_dormancy(both_upper = both_upper, 
                                          both_lower = both_lower, 
                                          raster_high = raster_terr_high, 
                                          raster_low = raster_terr_low)

# c_high <- combined_elev_x_dormancy[[1]]
# c_low <- combined_elev_x_dormancy[[2]]
# rr_high <- combined_elev_x_dormancy[[3]]
# rr_low <- combined_elev_x_dormancy[[4]]
# 
# plot(c_high[[1]])
# plot(c_low[[1]])
# plot(rr_high[[1]])
# plot(rr_low[[1]])

## restrict range to contiguous habitat that begins at the species realized range:
prs_terrestrial <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined_reg, 
                                           type = 'reg')
prs_terrestrial_dormancy <- create_potential_ranges(clumped_temps = clumped_temps, 
                                           realized_ranges = realized_ranges, combined = combined_dormancy,
                                           type = 'dormancy')
prs_terrestrial_elev <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                realized_ranges = realized_ranges, combined = combined_elev,
                                                type = 'elevation')
prs_terrestrial_elev_x_dormancy <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                realized_ranges = realized_ranges, combined = combined_elev_x_dormancy,
                                                type = 'elev_x_dormancy')


##saveRDS(prs_terrestrial, "./data-processed/prs_terrestrial.rds")
prs_terrestrial <- readRDS("./data-processed/prs_terrestrial.rds")
##saveRDS(prs_terrestrial_dormancy, "./data-processed/prs_terrestrial_dormancy_6mo.rds")
prs_terrestrial_dormancy <- readRDS("./data-processed/prs_terrestrial_dormancy_6mo.rds")
##saveRDS(prs_terrestrial_elev, "./data-processed/prs_terrestrial_elev.rds")
prs_terrestrial_elev <- readRDS("./data-processed/prs_terrestrial_elev.rds")
##saveRDS(prs_terrestrial_elev_x_dormancy, "./data-processed/prs_terrestrial_elev_x_dormancy.rds")
prs_terrestrial_elev_x_dormancy <- readRDS("./data-processed/prs_terrestrial_elev_x_dormancy.rds")

## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_reg_terr <- extract_thermal_niche(type = "reg", prs = prs_terrestrial, realm = 'terrestrial')
niche_dormancy_terr <- extract_thermal_niche(type = "dormancy", prs = prs_terrestrial_dormancy, 
                                             realm = 'terrestrial')
niche_elev_terr <- extract_thermal_niche(type = "elevation", prs = prs_terrestrial_elev, realm = 'terrestrial')
niche_elev_x_dormancy_terr <- extract_thermal_niche(type = "elev_x_dormancy", prs = prs_terrestrial_elev_x_dormancy,
                                                    realm = 'terrestrial')

## save:
write.csv(niche_reg_terr, "./data-processed/thermal-niche/niche_reg_terr.csv", row.names = FALSE)
write.csv(niche_dormancy_terr, "./data-processed/thermal-niche/niche_dormancy_terr.csv", row.names = FALSE)
write.csv(niche_elev_terr, "./data-processed/thermal-niche/niche_elev_terr.csv", row.names = FALSE)
write.csv(niche_elev_x_dormancy_terr, "./data-processed/thermal-niche/niche_elev_x_dormancy_terr.csv", 
          row.names = FALSE)



## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper, 
#                                            raster_high = raster_terr_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower, 
#                                           raster_low = raster_terr_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_terrestrial_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                                realized_ranges = realized_ranges,
#                                                                high_filtered = high_filtered,  
#                                                                low_filtered = low_filtered)



#######################################################
#####                MARINE SPECIES:             ######
#######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Marine")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Marine")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_marine_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
combined_reg <- filter_by_tolerance(both_upper = both_upper,
                                     both_lower = both_lower,
                                     raster_high = raster_marine_high, 
                                     raster_low = raster_marine_low)


## restrict range to contiguous habitat that begins at the species realized range:
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, combined = combined_reg,
                                      type = 'reg')

##saveRDS(prs_marine, "./data-processed/prs_marine.rds")
prs_marine <- readRDS("./data-processed/prs_marine.rds")

## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_reg_marine <- extract_thermal_niche(type = "reg", prs_marine, realm = 'marine')

write.csv(niche_reg_marine, "./data-processed/thermal-niche/niche_reg_marine.csv", row.names = FALSE)


## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
#                                            raster_high = raster_marine_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
#                                           raster_low = raster_marine_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_marine_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                           realized_ranges = realized_ranges, 
#                                                           high_filtered = high_filtered,  
#                                                           low_filtered = low_filtered)


######################################################
#####            INTERTIDAL SPECIES:            ######
######################################################
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Intertidal")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Intertidal")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]
only_upper <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps$geometry)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
combined_reg <- filter_by_tolerance(both_upper = both_upper,
                                     both_lower = both_lower,
                                     raster_high = raster_intertidal_high,
                                     raster_low = raster_intertidal_low)


## restrict range to contiguous habitat that begins at the species realized range:
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, combined = combined_reg,
                                          type = 'reg')

##saveRDS(prs_intertidal, "./data-processed/prs_intertidal.rds")
prs_intertidal <- readRDS("./data-processed/prs_intertidal.rds")

## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_reg_intertidal <- extract_thermal_niche(type = "reg", prs_intertidal, realm = 'intertidal')

write.csv(niche_reg_intertidal, "./data-processed/thermal-niche/niche_reg_intertidal.csv", row.names = FALSE)


## SPECIES WITH ONLY ONE THERMAL LIMIT:
# high_filtered <- filter_by_tolerance_upper(only_upper = only_upper,
#                                            raster_high = raster_intertidal_high)
# plot(high_filtered)
# low_filtered <- filter_by_tolerance_lower(only_lower = only_lower,
#                                           raster_low = raster_intertidal_low)
# plot(low_filtered)
# 
# ## restrict range to contiguous habitat that begins at the species realized range:
# prs_intertidal_one_limit <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                               realized_ranges = realized_ranges, 
#                                                               high_filtered = high_filtered,  
#                                                               low_filtered = low_filtered)


########################################################
###     COLLATING POTENTIAL RANGES ACROSS REALMS     ###
########################################################
potential_ranges_both_limits <- stack(prs_terrestrial, prs_marine, prs_intertidal)
potential_ranges_both_limits_dormancy <- stack(prs_terrestrial_dormancy, prs_marine, prs_intertidal)
potential_ranges_both_limits_elev <- stack(prs_terrestrial_elev, prs_marine, prs_intertidal)
potential_ranges_both_limits_elev_x_dormancy <- stack(prs_terrestrial_elev_x_dormancy, prs_marine, prs_intertidal)
##saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges.rds")
##saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges_notcutatequator.rds")
##saveRDS(potential_ranges_both_limits_dormancy, "data-processed/potential_ranges_notcutatequator_dormancy.rds")
##saveRDS(potential_ranges_both_limits_elev, "data-processed/potential_ranges_notcutatequator_elev.rds")
##saveRDS(potential_ranges_both_limits_elev_x_dormancy, "data-processed/potential_ranges_notcutatequator_elev_x_dormancy.rds")

#potential_ranges_one_limit <- stack(prs_terrestrial_one_limit, prs_marine_one_limit, prs_intertidal_one_limit)
##saveRDS(potential_ranges_one_limit, "data-processed/potential_ranges_one_limit_coldwarm.rds")


## COLLATING TEMPS IN THERMAL NICHE ACROSS REALMS
##################################################
niche_reg <- rbind(niche_reg_terr, niche_reg_marine, niche_reg_intertidal)
niche_elev <- rbind(niche_elev_terr, niche_reg_marine, niche_reg_intertidal)
niche_dormancy <- rbind(niche_dormancy_terr, niche_reg_marine, niche_reg_intertidal)
niche_elev_x_dormancy <- rbind(niche_elev_x_dormancy_terr, niche_reg_marine, 
                               niche_reg_intertidal)

write.csv(niche_reg, "./data-processed/thermal-niche/niche_reg.csv", row.names = FALSE)
write.csv(niche_elev, "./data-processed/thermal-niche/niche_elev.csv", row.names = FALSE)
write.csv(niche_dormancy, "./data-processed/thermal-niche/niche_dormancy.csv", row.names = FALSE)
write.csv(niche_elev_x_dormancy, "./data-processed/thermal-niche/niche_elev_x_dormancy.csv", row.names = FALSE)

######################################################
##                  FUNCTIONS:                      ##
######################################################

## functions to filter rasterized temperature data (raster_high, raster_low) by thermal tolerance limits (either upper, lower, or both) 
filter_by_tolerance <- function(both_upper, 
                                     both_lower, 
                                     raster_high, 
                                     raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low

  species = 1
  while (species < nrow(both_upper) + 1) {
    h <- raster_high$seasonal_high_temp
    l <- raster_low$seasonal_low_temp
    
    rr_high <- addLayer(rr_high, h) 
    rr_low <- addLayer(rr_low, l) 
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
    below_ctmax <- both_upper$thermal_limit[species] - h
    h[below_ctmax < 0] <- NA
    above_ctmin <- both_lower$thermal_limit[species] - l
    l[above_ctmin > 0] <- NA
    
    high <- addLayer(high, h) 
    low <- addLayer(low, l) 
  
    species = species + 1
  }
  
  names(high) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
  names(low) <- paste(both_upper$Genus, both_upper$Species, sep = "_")
  
  ## combine to find cells where seasonal high temp is less than CTmax AND seasonal low temp is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low))
}

filter_by_tolerance_dormancy <- function(both_upper, 
                                                   both_lower, 
                                                   raster_high, 
                                                   raster_low) {
  
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  species = 1
  high_names <- c()
  low_names <- c()
  tracker <- c()
  while (species < nrow(both_upper) + 1) {
    species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
      h <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
      l <- min(stack(raster_low$cold_dormancy_6mo,  raster_low$hot_dormancy_6mo), na.rm = TRUE)
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "both")
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      h <- raster_high$cold_dormancy_6mo
      l <-  raster_low$cold_dormancy_6mo
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "cold")
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      h <- raster_high$hot_dormancy_6mo
      l <- raster_low$hot_dormancy_6mo
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
  
      tracker <- append(tracker, rep("hot"))
    }
    else {
      h <- raster_high$seasonal_high_temp
      l <- raster_low$seasonal_low_temp
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "neither")
    }
    
    species = species + 1
  }

  names(high) <- paste(both_upper$Genus, both_upper$Species, tracker, sep = "_")
  names(low) <- paste(both_upper$Genus, both_upper$Species, tracker, sep = "_")
  ##plot(high)
  ##plot(low)
  
  ## combine to find cells where seasonal high temp is less than CTmax AND seasonal low temp is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
              
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low))
}

filter_by_tolerance_elev <- function(both_upper, 
                                     both_lower, 
                                     raster_high, 
                                     raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  species = 1
  while (species < nrow(both_upper) + 1) {
    h <- raster_high$seasonal_high_temp
    l <-  raster_low$low_at_min_elev
    
    rr_high <- addLayer(rr_high, h) 
    rr_low <- addLayer(rr_low, l) 
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
    below_ctmax <- both_upper$thermal_limit[species] - h
    h[below_ctmax < 0] <- NA
    above_ctmin <- both_lower$thermal_limit[species] - l
    l[above_ctmin > 0] <- NA
    
    high <- addLayer(high, h) 
    low <- addLayer(low, l) 
    
    species = species + 1
  }
  
  names(high) <- paste(both_upper$Genus, both_upper$Species, 'elev', sep = "_")
  names(low) <- paste(both_upper$Genus, both_upper$Species, 'elev', sep = "_")
  
  ## combine to find cells where seasonal high temp is less than CTmax AND seasonal low temp is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low))
}

filter_by_tolerance_elev_x_dormancy <- function(both_upper, 
                                     both_lower, 
                                     raster_high, 
                                     raster_low) {
  
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  species = 1
  high_names <- c()
  low_names <- c()
  tracker <- c()
  while (species < nrow(both_upper) + 1) {
    species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
      h <- max(stack(raster_high$hot_dormancy_6mo,  raster_high$cold_dormancy_6mo), na.rm = TRUE)
      l <- min(stack(raster_low$low_elev_x_cold_dormancy, raster_low$low_elev_x_hot_dormancy), na.rm = TRUE)
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "both")
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      h <- raster_high$cold_dormancy_6mo
      l <-  raster_low$low_elev_x_cold_dormancy
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "cold")
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      h <- raster_high$hot_dormancy_6mo
      l <-  raster_low$low_elev_x_hot_dormancy
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "cold")
    }
    else {
      h <- raster_high$seasonal_high_temp
      l <-  raster_low$low_at_min_elev
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- both_upper$thermal_limit[species] - h
      h[below_ctmax < 0] <- NA
      above_ctmin <- both_lower$thermal_limit[species] - l
      l[above_ctmin > 0] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      tracker <- append(tracker, "cold")
    }
    
    species = species + 1
  }

  names(high) <- paste(both_upper$Genus, both_upper$Species, tracker, 'elev', sep = "_")
  names(low) <- paste(both_upper$Genus, both_upper$Species, tracker, 'elev', sep = "_")
  ##plot(high)
  ##plot(low)
  
  ## combine to find cells where seasonal high temp is less than CTmax AND seasonal low temp is greater than CTmin 
  combined_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                          crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  combined_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                         crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  i = 1  
  while (i < nlayers(high) + 1) {
    combined_high <- addLayer(combined_high, mask(high[[i]], low[[i]]), updatevalue = NA)
    combined_low <- addLayer(combined_low, mask(low[[i]], high[[i]]), updatevalue = NA)
    
    i = i + 1
  }
  names(combined_high) <- names(high)
  names(combined_low) <- names(low)
  names(rr_high) <- names(high)
  names(rr_low) <- names(low)
  #plot(combined_high[[1]])
  #plot(combined_low[[1]])
  
  return (list(combined_high, combined_low, rr_high, rr_low))
}

filter_by_tolerance_upper <- function(only_upper, 
                                      raster_high) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  ## for upper limits, only filter use seasonal high temps 
  species = 1
  while (species < nrow(only_upper) + 1) {
    high <- addLayer(high, raster_high[[1]] - only_upper$thermal_limit[species]) 
    
    species = species + 1
  }
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  high[high > 0] <- NA
  names(high) <- c(paste(only_upper$Genus, only_upper$Species, sep = "_"))
  ##plot(high)
  
  return(high)
}

filter_by_tolerance_lower <- function(only_lower, 
                                      raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## for lower limits, only filter use seasonal low temps 
  species = 1
  while (species < nrow(only_lower) + 1) {
    low <- addLayer(low, raster_low[[1]] - only_lower$thermal_limit[species]) 
    
    species = species + 1
  }
  
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  low[low < 0] <- NA
  if (!is.na(minValue(low)[1])) {
    names(low) <- c(paste(only_lower$Genus, only_lower$Species, sep = "_"))
  }
  ##plot(low)
  
  return(low)
}


## creates potential range rasterstack containing a layer for each realized range that represents the species potential range containing only contiguous habitat in clumps that their realized range touches 
## returns 3 versions of potential ranges (boolean values, high temperature values and low temperature values) 
create_potential_ranges <- function (clumped_temps, 
                                     combined, 
                                     realized_ranges, 
                                     type) {
  
  combined_high <- combined[[1]]
  combined_low <- combined[[2]]
  
  rr_high <- combined[[3]]
  rr_low <- combined[[4]]
  
  rr_high_new <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                        crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  rr_low_new <- rr_high_new
  
  pr_restricted_bool <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  pr_restricted_high <- pr_restricted_bool
  pr_restricted_low <- pr_restricted_bool
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  names <- c()
  i = 1
  while (i < length(names(combined_high)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- combined_high[[i]] 
    potential_range <- potential_raster %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf()
    
    ## get realized range and realized range raster
    if(type == 'elevation') {
      split <- str_split_fixed(names(potential_raster), '_', n = 3) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
    }
    else if (type == 'elev_x_dormancy') {
      split <- str_split_fixed(names(potential_raster), '_', n = 4) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
    }
    else if(type == 'dormancy') {
      split <- str_split_fixed(names(potential_raster), "_", n = 3)
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    else {
      species <- names(potential_raster) %>%
        str_replace_all("_", " ") 
    }
  
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") 
    
    ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      ## overlay realized range with potential range
      ## restrict potential range to clumps that overlap with the realized range
      intersects <- st_intersects(clumped_temps, realized_range[num,], sparse = FALSE)[,]
      rr <- filter(clumped_temps, intersects == TRUE)
      
      intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
      
      if (nrow(rr) == 1) {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE),]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover = TRUE)
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      else {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE)[,2],]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      
      ## create boolean version of potential range 
      pr_restricted_bool <- addLayer(pr_restricted_bool, pr_restricted, updatevalue = NA)

      ## create version of prs with temperatures as raster values for thermal niche metrics 
      pr_restricted_high <- addLayer(pr_restricted_high, 
                                     mask(combined_high[[i]], pr_restricted, updatevalue = NA), 
                                     updatevalue = NA)
      pr_restricted_low <- addLayer(pr_restricted_low, 
                                     mask(combined_low[[i]], pr_restricted, updatevalue = NA), 
                                     updatevalue = NA)
      
      ## add layer to new rr:
      rr_high_new <- addLayer(rr_high_new, rr_high[[i]])
      rr_low_new <- addLayer(rr_low_new, rr_low[[i]])
      
      names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
      
      num = num + 1
    }
    
    i = i + 1
  }
  names(pr_restricted_bool) <- names
  names(pr_restricted_high) <- names
  names(pr_restricted_low) <- names
  names(rr_high_new) <- names
  names(rr_low_new) <- names
  #plot(pr_restricted_bool)
  
  return(list(pr_restricted_bool, pr_restricted_high, pr_restricted_low, rr_high_new, rr_low_new))
}


## returns a potential range raster containing a layer for all species in high_filtered and low_filtered that represents their potential equatorward and poleward range respectively containing only contiguous habitat in clumps that their realized range touches 
create_potential_ranges_one_limit <- function (clumped_temps, 
                                               realized_ranges, 
                                               high_filtered,  
                                               low_filtered) {
  ## WARM RANGES:
  pr_restricted_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  names <- c()
  i = 1
  if (nlayers(high_filtered) > 1) {
    while (i < length(names(high_filtered)) + 1) {
      ## get unrestricted potential range 
      potential_raster <- high_filtered[[i]] 
      potential_range <- potential_raster %>%
        clump(., directions = 8) %>%
        rasterToPolygons(., dissolve = TRUE) %>%
        st_as_sf()
      
      ## get realized range
      species <- names(potential_raster) %>%
        str_replace_all("_", " ")
      
      realized_range <- realized_ranges[which(as.character(realized_ranges$species)
                                                    %in% species),] 
      st_crs(realized_range) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
      
      ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
      num = 1
      while (num < nrow(realized_range)+1) {
        ## overlay realized range with potential range
        ## restrict potential range to clumps that overlap with the realized range
        intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
        rr <- filter(clumped_temps, intersects == TRUE)
        
        intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
        
        if (is.null(intersects_potential)) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- r
        }
        else if (nrow(rr) == 1) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        else {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE)[,2],]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        
        ## if species range does not cross the equator, cut off potential range at the equator based on where realized range latitudinal midpoint is (below or above equator)
        if (realized_range$hemisphere[num] == "EQUATOR") {
          p <- Polygon(matrix(c(-180,-90,-180,90,180,90,180,-90,-180,-90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
        }
        else if (realized_range$hemisphere[num] == "N") {
          p <- Polygon(matrix(c(-180, 90,-180,0,180,0,180, 90,-180, 90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1"))) 
        }
        else {
          p <- Polygon(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),
                              ncol=2, byrow=TRUE))
          rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
        }
        
        rect_raster <- rasterize(rect, r, getCover=TRUE)
        rect_raster[rect_raster==0] <- NA
        
        plot(pr_restricted, col = "orange", main = paste(realized_range$species[num], "warm", 
                                                         realized_range$hemsphere[num]))
        
        pr_restricted <- mask(pr_restricted, rect_raster)
        
        plot(pr_restricted, add = TRUE)
        plot(rect, add = TRUE)
        ##plot(pr_multipolygons,  add=TRUE, col = "purple")
        ##plot(potential_range,  add=TRUE, col = "yellow")
        
        dev.copy(png, filename = paste("figures/one-limits/", realized_range$species[num], "_warm_",
                                       realized_range$source[num], ".png", sep = ""), 
                 width = 614, height = 472);
        dev.off()
        
        pr_restricted_high <- addLayer(pr_restricted_high, pr_restricted, updatevalue = NA)
        names <- append(names, paste(species, realized_range$source[num], "warm", sep = "_"))
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_high) <- names
    ##plot(pr_restricted_high)
  }
  
  
  ## COLD RANGES:
  pr_restricted_low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  names <- c()
  i = 1
  if (!is.na(minValue(low_filtered)[1])) {
    while (i < length(names(low_filtered)) + 1) {
      ## get unrestricted potential range 
      potential_raster <- low_filtered[[i]] 
      potential_range <- potential_raster %>%
        clump(., directions = 8) %>%
        rasterToPolygons(., dissolve = TRUE) %>%
        st_as_sf()
      
      ## get realized range
      species <- names(potential_raster) %>%
        str_replace_all("_", " ")
      
      realized_range <- realized_ranges[which(as.character(realized_ranges$species)
                                                    %in% species),] 
      st_crs(realized_range) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
      
      ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
      num = 1
      while (num < nrow(realized_range)+1) {
        ## overlay realized range with potential range
        ## restrict potential range to clumps that overlap with the realized range
        intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
        rr <- filter(clumped_temps, intersects == TRUE)
        
        intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
        
        if (nrow(rr) == 1) {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE),]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        else {
          pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                    arr.ind=TRUE)[,2],]
          pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
          pr_rasterized[pr_rasterized == 0] <- NA
          pr_restricted <- mask(pr_rasterized, pr_multipolygons)
        }
        
        
        ## mark hemisphere since later equator species will need poleward range filling calculated twice
          plot(pr_restricted, col = "orange", main = paste(realized_range$species[num], "cold", 
                                                           realized_range$hemisphere[num]),
               legend = FALSE)
         
          ##plot(pr_multipolygons,  add=TRUE, col = "purple")
          ##plot(potential_range,  add=TRUE, col = "yellow")
          
          dev.copy(png, filename = paste("figures/one-limits/", realized_range$species[num], "_cold_",
                                         realized_range$source[num], ".png", sep = ""), 
                   width = 614, height = 472);
          dev.off()
          
          pr_restricted_low <- addLayer(pr_restricted_low, pr_restricted, updatevalue = NA)
          names <- append(names, paste(species, realized_range$source[num], "cold", 
                                       sep = "_"))

        
        num = num + 1
      }
      
      i = i + 1
    }
    names(pr_restricted_low) <- names
    plot(pr_restricted_low)
  }
  
  if(!is.na(minValue(high_filtered)[1]) & !is.na(minValue(low_filtered)[1])){
    pr_restricted_all <- addLayer(pr_restricted_high, pr_restricted_low)
  }
  else if(!is.na(minValue(low_filtered)[1]) & is.na(minValue(high_filtered)[1])) {
    pr_restricted_all <- pr_restricted_low
  }
  else {
    pr_restricted_all <- pr_restricted_high
  }
  
  return(pr_restricted_all)
}


extract_thermal_niche <- function(type, prs, realm) {
  
  rrs <- readRDS("data-processed/rasterized_rrs.rds")
  names(rrs) <- str_replace_all(names(rrs), '\\.', '_')
  
  traits <- read.csv("./data-processed/wrangled-traits.csv")
  
  prs_high <- prs[[2]]
  prs_low <- prs[[3]]
  
  rrs_high <- prs[[4]]
  rrs_low <- prs[[5]]
  
  range = 1
  while (range < nlayers(prs_high) + 1) {
    
    ## get species info
    ## get realized range and realized range raster
    if(type == 'elevation') {
      split <- str_split_fixed(names(prs_high)[[range]], '_', n = 4) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
      range_id <- paste(split[1,1], split[1,2], split[1,4], sep = "_")
    }
    else if (type == 'elev_x_dormancy') {
      split <- str_split_fixed(names(prs_high)[[range]], '_', n = 5) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
      range_id <- paste(split[1,1], split[1,2], split[1,5], sep = "_")
    }
    else if(type == 'dormancy') {
      split <- str_split_fixed(names(prs_high)[[range]], "_", n = 4)
      species <- paste(split[1,1], split[1,2], sep = " ")
      range_id <- paste(split[1,1], split[1,2], split[1,4], sep = "_")
    }
    else {
      split <- str_split_fixed(names(prs_high)[[range]], '_', n = 3) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
      range_id <- names(prs_high)[[range]]
    }

    species_traits <- traits[which(traits$genus_species == species),]
    
    cold_dormancy <- species_traits$cold_season_dormancy_ == "Yes" & !is.na(species_traits$cold_season_dormancy_)
    hot_dormancy <- species_traits$hot_season_dormancy_ == "Yes" & !is.na(species_traits$hot_season_dormancy_)
    dormancy <- cold_dormancy || hot_dormancy
    dormancy_type <- ifelse(!dormancy, "none", ifelse(cold_dormancy & hot_dormancy, "both", 
                                                      ifelse(cold_dormancy, "cold", "hot")))
    ## make sure pr is not empty:
    pr_high_empty <- length(which(is.na(values(prs_high[[range]])) == FALSE)) <= 0
    pr_low_empty <- length(which(is.na(values(prs_low[[range]])) == FALSE)) <= 0
    
    ## get temps within potential range 
    if (!pr_high_empty) {
      high_vals <- data.frame(high_or_low = "high", 
                              temps = values(prs_high[[range]])[which(is.na(values(prs_high[[range]])) == FALSE)],
                              type = paste("potential", type, sep = "_"),
                              dormancy = dormancy_type)
    }
    else {
      high_vals <- data.frame(high_or_low = "high", 
                              temps = NA,
                              type = paste("potential", type, sep = "_"),
                              dormancy = dormancy_type)
    }
    if (!pr_low_empty) {
      low_vals <- data.frame(high_or_low = "low", 
                             temps = values(prs_low[[range]])[which(is.na(values(prs_low[[range]])) == FALSE)],
                             type = paste("potential", type, sep = "_"),
                             dormancy = dormancy_type)
    }
    else {
      low_vals <- data.frame(high_or_low = "low", 
                             temps = NA,
                             type = paste("potential", type, sep = "_"),
                             dormancy = dormancy_type)
    }
  
    ## get temps within corresponding rasterized realized range
    rr <- rrs[[which(names(rrs) == range_id)]]
    
    rrs_high[[range]] <- mask(rrs_high[[range]], rr, updatevalue = NA)
    rrs_low[[range]] <- mask(rrs_low[[range]], rr, updatevalue = NA)
    
    ## make sure rr is not empty:
    rr_high_empty <- length(which(is.na(values(rrs_high[[range]])) == FALSE)) <= 0
    rr_low_empty <- length(which(is.na(values(rrs_low[[range]])) == FALSE)) <= 0
    
    ## get temps within realized range 
    if (!rr_high_empty) {
      high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                              temps = values(rrs_high[[range]])[which(is.na(values(rrs_high[[range]])) == FALSE)],
                              type = paste("realized", type, sep = "_"),
                              dormancy = dormancy_type))
    }
    else {
      high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                              temps = NA,
                              type = paste("realized", type, sep = "_"),
                              dormancy = dormancy_type))
    }
    if (!rr_low_empty) {
      low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                             temps = values(rrs_low[[range]])[which(is.na(values(rrs_low[[range]])) == FALSE)],
                             type = paste("realized", type, sep = "_"),
                             dormancy = dormancy_type))
    }
    else {
      low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                             temps = NA,
                             type = paste("realized", type, sep = "_"),
                             dormancy = dormancy_type))
    }
    
    ## bind dfs together:
    if(range == 1) {
      thermal_niche <- rbind(high_vals, low_vals) %>%
        mutate(range = range_id) %>%
        select(range, everything()) 
    }
    else {
      thermal_niche <- rbind(high_vals, low_vals) %>%
        mutate(range = range_id) %>%
        select(range, everything()) %>%
        rbind(thermal_niche, .)
    }
    
    range = range + 1
  }
 
  thermal_niche$realm <- realm
  
  return(thermal_niche)
}


##############################################
#####               PLOTS               ######
##############################################
## function that plots potential ranges and corresponding realized ranges into a folder 
plot_ranges_overlap  <- function (folder, potential_ranges) {
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  prs <- as.data.frame(potential_ranges, xy=TRUE)
  colnames(prs)[1:2] <- c("longitude", "latitude")
  
  land <- as.data.frame(raster("./data-processed/raster_terr_mask.grd"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("./data-processed/raster_marine_mask.grd"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("./data-processed/raster_intertidal_mask.grd"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  traits <- read.csv("./data-processed/wrangled-traits.csv")
  
  ## ggplots of all:
  i = 1
  while (i < ncol(prs) - 1) {
    range <- colnames(prs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("_", ".")
    
    if(TRUE %in% str_detect(names(potential_ranges), "both")) {
      split <- str_split_fixed(range, "\\_", n = 4)
      source <- split[1,4]
    } 
    else {
      split <- str_split_fixed(range, "\\_", n = 3)
      source <- split[1,3]
    }
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(traits$genus_species == paste(split[1,1], split[1,2], sep = " "))]
    
    if(realm == "Terrestrial" | realm == "Freshwater") {
      backdrop = land
    }
    else if (realm == "Marine") {
      backdrop = ocean
    }
    else {
      backdrop = intertidal
    }
    
    rr_index <- which(str_detect(colnames(rrs), species) & str_detect(colnames(rrs), source))
    
    r <- rrs[,c(1:2, rr_index)] 
    colnames(r)[3] <- "rr"
    r <- left_join(r, prs[,c(1:2, which(colnames(prs) == range))]) 
    colnames(r)[4] <- "pr"
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill=as.factor(rr))) + 
      scale_fill_manual(values = c("yellow"), aesthetics = 'fill', labels = ) +
      annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.5,
               fill = r$pr) +
      annotate(geom="raster", x=land$longitude, y=land$latitude, alpha=.1,
               fill = backdrop$mask) +
      labs(title = range,
           y = "Latitude",
           x = "Longitude") +
      scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
      scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
      theme(legend.position = "none")
    
    
    ## write to file:
    ggsave(r_gg, path = paste("figures/range-plots/", folder, "/", sep = ""), 
           filename = paste(range, ".png", sep = ""), 
           height = 6, width = 10, units = "in", device = "png")
    
    i = i + 1
  }
  
}

## plot the ranges:
potential_ranges_both_limits <- readRDS("data-processed/potential_ranges_notcutatequator.rds")
potential_ranges_both_limits_dormancy <- readRDS("data-processed/potential_ranges_notcutatequator_dormancy.rds")
potential_ranges_both_limits_elev <- readRDS("data-processed/potential_ranges_notcutatequator_elev.rds")

plot_ranges_overlap(folder = "no-modification", potential_ranges = potential_ranges_both_limits)
plot_ranges_overlap(folder = "dormancy", potential_ranges = potential_ranges_both_limits_dormancy)
plot_ranges_overlap(folder = "elev", potential_ranges = potential_ranges_both_limits_elev)




##############################################
#####               GARBAGE             ######
##############################################
## writes plots of realized range and restricted potential range to files for species with both thermal limits 
plot_ranges <- function (clumped_temps, 
                         combined, 
                         realized_ranges) {
  i = 1
  while (i < length(names(combined)) + 1) {
    ## get unrestricted potential range 
    potential_range <- combined[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    if(str_detect(species, "dormancy") ==  TRUE) {
      split <- str_split_fixed(species, " ", n = 3)
      
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/plot_ranges/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/plot_ranges/", names(potential_range), "_potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
}


## writes plots of realized range and restricted potential poleward or equatorward range to files for species with only one thermal limit
plot_ranges_one_limit <- function (clumped_temps, 
                                   realized_ranges, 
                                   high_filtered, 
                                   low_filtered) {
  i = 1
  while (i < length(names(high_filtered)) + 1) {
    ## get unrestricted potential range 
    potential_range <- high_filtered[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_equatorward-potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
  
  i = 1
  while (i < length(names(low_filtered)) + 1 & length(names(low_filtered)) > 1) {
    ## get unrestricted potential range 
    potential_range <- low_filtered[[i]] 
    potential_raster <- potential_range %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf() 
    
    ## get realized range:
    species <- names(potential_range) %>%
      str_replace_all("_", " ")
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ##plot(st_geometry(realized_range), add = TRUE, col = "red")
    
    ## overlay realized range with potential range and restrict potential range to clumps that overlap with the realized range
    intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
    sub <- filter(clumped_temps, intersects == TRUE)
    plot(sub)
    
    intersects_potential <- st_intersects(sub, potential_raster, sparse = FALSE)[,]
    if (nrow(sub) == 1) {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE),]
    }
    else {
      sub_potential <- potential_raster[which(intersects_potential == TRUE, arr.ind=TRUE)[,2],]
    }
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - realized range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(ms_simplify(realized_range, keep_shapes = TRUE)), add = TRUE, col = "red")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_realized-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    plot(st_geometry(clumped_temps), main = paste(names(potential_range), 
                                                  " - potential range", sep = ""))
    plot(sub, add = TRUE, col = "blue")
    plot(st_geometry(sub_potential), add = TRUE, col = "yellow")
    
    dev.copy(png, filename = paste("figures/selecting-contiguous-patches/", names(potential_range), "_poleward-potential-range.png", sep = ""), width = 1000, height = 500, res = 200);
    dev.off()
    
    i = i + 1
  }
}







## this function was used to make figures with complicated titles for species with dormancy
## simplified to match other sets of ranges 
filter_by_tolerance_dormancy <- function(both_upper, 
                                         both_lower, 
                                         raster_high, 
                                         raster_low) {
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  species = 1
  high_names <- c()
  low_names <- c()
  tracker <- c()
  while (species < nrow(both_upper) + 1) {
    species_traits <- traits[which(traits$genus_species == both_upper$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" & species_traits$hot_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1:2]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1:2]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                             c(0,6), "months_both_dormancy", sep = "_"))
      low_names <- append(low_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                           c(0,6), "months_both_dormancy", sep = "_"))
      tracker <- append(tracker, rep("both", 2))
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1:2]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, rep(paste(species_traits$Genus, species_traits$Species, 
                                                 "no_hot_dormancy", sep = "_"), 2))
      low_names <- append(low_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                           c(0,6), "months_cold_dormancy", sep = "_"))
      tracker <- append(tracker, rep("cold", 2))
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      high <- addLayer(high, raster_high[[1:2]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species])
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species])
      high_names <- append(high_names, paste(paste(species_traits$Genus, species_traits$Species, sep = "_"),
                                             c(0,6), "months_hot_dormancy", sep = "_"))
      low_names <- append(low_names, rep(paste(species_traits$Genus, species_traits$Species, "no_cold_dormancy",
                                               sep = "_"),2))
      tracker <- append(tracker, rep("hot", 2))
    }
    else {
      high <- addLayer(high, raster_high[[1]] - both_upper$thermal_limit[species]) 
      low <- addLayer(low, raster_low[[1]] - both_lower$thermal_limit[species]) 
      high_names <- append(high_names, paste(species_traits$Genus, species_traits$Species, 
                                             "no_dormancy", sep = "_"))
      low_names <- append(low_names, paste(species_traits$Genus, species_traits$Species, "no_dormancy",
                                           sep = "_"))
      tracker <- append(tracker, "neither")
    }
    
    species = species + 1
  }
  
  ## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
  high[high > 0] <- NA
  low[low < 0] <- NA
  names(high) <- high_names
  names(low) <- low_names
  ##plot(high)
  ##plot(high)
  
  ## combine to find cells where seasonal high temp is less than CTmax and seasonal low temp is greater than CTmin 
  combined <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  c_names <- c()
  i = 1  
  while (i < nlayers(high) + 1) {
    combined <- addLayer(combined, mask(high[[i]], low[[i]]), updatevalue = NA)
    c_names <- append(c_names, ifelse(tracker[i] == "both", high_names[i],
                                      ifelse(tracker[i] == "cold", low_names[i], 
                                             ifelse(tracker[i] == "hot", high_names[i],
                                                    ifelse(tracker[i] == "neither", low_names[i],
                                                           "oop")))))
    i = i + 1
  }
  names(combined) <- c_names
  ##plot(combined)
  
  return (combined)
}
create_potential_ranges <- function (clumped_temps, 
                                     combined, 
                                     realized_ranges) {
  
  pr_restricted_all <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  names <- c()
  i = 1
  while (i < length(names(combined)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- combined[[i]] 
    potential_range <- potential_raster %>%
      clump(., directions = 8) %>%
      rasterToPolygons(., dissolve = TRUE) %>%
      st_as_sf()
    
    ## get realized range
    species <- names(potential_raster) %>%
      str_replace_all("_", " ") 
    
    if(str_detect(species, "dormancy") ==  TRUE) {
      split <- str_split_fixed(species, " ", n = 3)
      
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
    
    ## if species has multiple realized ranges (IUCN and GBIF), loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      ## overlay realized range with potential range
      ## restrict potential range to clumps that overlap with the realized range
      intersects <- st_intersects(clumped_temps, realized_range, sparse = FALSE)[,]
      rr <- filter(clumped_temps, intersects == TRUE)
      
      intersects_potential <- st_intersects(rr, potential_range, sparse = FALSE)[,]
      
      if (nrow(rr) == 1) {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE),]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover = TRUE)
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      else {
        pr_multipolygons <- potential_range[which(intersects_potential == TRUE, 
                                                  arr.ind=TRUE)[,2],]
        pr_rasterized <- rasterize(pr_multipolygons, r, getCover=TRUE) 
        pr_restricted <- mask(pr_rasterized, pr_multipolygons)
      }
      
      ## if realized range does not cross equator, cut off potential range at the equator:
      ## DEC 14: commented out and rerun
      #lat_mp <- realized_range$lat_mp[num]
      
      # if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp > 0) {
      #   p <- Polygon(matrix(c(-180,0,-180,90,180,90,180,0,-180,0),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      # }
      # else if (!realized_range$hemisphere[num] == "EQUATOR" & lat_mp < 0) {
      #   p <- Polygon(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1"))) 
      # }
      # else {
      #   p <- Polygon(matrix(c(-180, 90,-180,-90,180,-90,180, 90,-180, 90),
      #                       ncol=2, byrow=TRUE))
      #   rect <- SpatialPolygons(list(Polygons(list(p), "p1")))
      # }
      # 
      # rect_raster <- rasterize(rect, r, getCover=TRUE)
      # rect_raster[rect_raster==0] <- NA
      
      #plot(pr_restricted, col = "orange")
      
      #pr_restricted <- mask(pr_restricted, rect_raster)
      
      ##plot(pr_restricted, col = "red", add = TRUE)
      
      ##plot(pr_restricted)
      ##plot(pr_multipolygons,  add=TRUE, col = "purple")
      ##plot(potential_range,  add=TRUE, col = "yellow")
      
      pr_restricted_all <- addLayer(pr_restricted_all, pr_restricted, updatevalue = NA)
      names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
      num = num + 1
    }
    
    i = i + 1
  }
  names(pr_restricted_all) <- names
  plot(pr_restricted_all)
  
  return(pr_restricted_all)
}