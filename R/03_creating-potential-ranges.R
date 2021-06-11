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
## read in function to rasterize realized ranges: 
source("R/shp2rast.R") ## function is shp2rast(rast, shp)



####################################################################################
#####                       READ IN TEMPERATURE RASTERS                       ######
####################################################################################
## read temperature data created in script 00
raster_terr_low <- stack("data-processed/raster_terr_low.grd")
raster_terr_high <- stack("data-processed/raster_terr_high.grd")
raster_marine_low <- stack("data-processed/raster_marine_low.grd")
raster_marine_high <- stack("data-processed/raster_marine_high.grd")
raster_intertidal_low <- stack("data-processed/raster_intertidal_low.grd")
raster_intertidal_high <- stack("data-processed/raster_intertidal_high.grd")
Te_sun_min <- stack("data-processed/Te_sun_min.grd")
Te_sun_max <- stack("data-processed/Te_sun_max.grd")
Te_shade_min <- stack("data-processed/Te_shade_min.grd")
Te_shade_max <- stack("data-processed/Te_shade_max.grd")

## read in realm mask layers:
t_mask <- raster("data-processed/raster_terr_mask.grd")
t_mask_nichemapr <- raster("data-processed/raster_terr_mask_nichemapr.grd")
m_mask <- raster("data-processed/raster_marine_mask.grd")
i_mask <- raster("data-processed/raster_intertidal_mask.grd")

## read in depth and elevation data:
elevs <- stack("data-processed/elevs_minmax.grd") 
depths <-  stack("data-processed/bathymetric-layers/depths_minmax.grd") 
shelf_mask <- raster("data-processed/bathymetric-layers/raster_200mdepth_mask.grd")

####################################################################################
#####                       READ IN THERMAL AND TRAIT DATA                    ######
####################################################################################
## read in thermal limits:
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges_taxized.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))

## read in acclimation data:
acc_data <- read.csv("data-processed/acclimation-data.csv")

## read in species traits:
traits <- read.csv("data-processed/rangetherm-traits_all-spp_filled.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "_"))



####################################################################################
#####                   RASTERIZE REALIZED RANGES                             ######
####################################################################################
## read in realized ranges:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
  mutate(range_id = paste(species, source, sep = "_")) 

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1, vals = 1)

i=1
while (i < nrow(realized_ranges)+1) {
  
  ## get speices' realized range:
  range <- realized_ranges[i, ]
  
  ## rasterize realized range:
  sp_range <- as_Spatial(range) 
  rr_raster <- shp2rast(shp = sp_range, rast = r) 
  rr_og <- rr_raster # store unrestricted version of rasterized range for comparison later
  
  ## get the species' traits:
  sp <- str_replace_all(traits$genus_species, "_", " ")
  sp_traits <- traits[which(sp == range$species), ]
  
  ## constrain realized range rasters by:
  ## - raster cells that were permitted to be in the potential range (in the realm)
  ## - elevational range for terrestrial species
  ## - depth distribution for marine and intertidal species
  if(range$realm == "Terrestrial") {
    rr_raster <- mask(rr_raster, t_mask) ## get rid of cells in realized range outside of the realm
    
    ## get elevational range limits 
    elev_max <- sp_traits$upper_elevation_limit_meters
    elev_min <- sp_traits$lower_elevation_limit_meters
    
    ## if both limits are not NA, restrict by both limits
    if (!is.na(elev_min) & !is.na(elev_max)) {
      rr_raster[elevs$elev_max < elev_min] <- NA 
      rr_raster[elevs$elev_min > elev_max] <- NA 
      ## exclude cells with max elev less than species min elevation
      ## exclude cells with min elev greater than species max elevation
    }
    else if (!is.na(elev_min)) {
      rr_raster[elevs$elev_max < elev_min] <- NA 
      ## exclude cells with max elev less than species min elevation
    }
    else if (!is.na(elev_max)) {
      rr_raster[elevs$elev_min > elev_max] <- NA 
      ## exclude cells with min elev greater than species max elevation
    }
  }
  else if(range$realm %in% c("Marine", "Intertidal")) {
    
    ## get rid of cells in realized range outside of the realm
    if (range$realm == "Marine") {
      rr_raster <- mask(rr_raster, m_mask) 
    }
    else if (range$realm == "Intertidal") {
      rr_raster <- mask(rr_raster, m_mask) 
    }
    
    ## get depth distribution limits 
    depth_upper <- sp_traits$upper_depth_limit
    depth_lower <- sp_traits$lower_depth_limit
    
    ## if both limits are not NA, restrict by both depth limits
    if (!is.na(depth_lower) & !is.na(depth_upper)) {
      rr_raster[depths$depth_min > depth_upper] <- NA 
      rr_raster[depths$depth_max < depth_lower] <- NA 
      ## exclude cells with deepest depth above species upper limit  
      ## exclude cells with shallowest depth less than species lower limit 
    }
    else if (!is.na(depth_upper)) {
      rr_raster[depths$depth_min > depth_upper] <- NA 
      ## exclude cells with deepest depth above species upper limit  
    }
    else if (!is.na(depth_lower)) {
      rr_raster[depths$depth_max < depth_lower] <- NA 
      ## exclude cells with shallowest depth less than species lower limit 
    }
    ## if lower depth limit is unknown, restrict to continental shelf if it is a coastal spp
    else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
      is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
      
      if(is_coastal) {
        rr_raster[is.na(shelf_mask)] <- NA
      }
    }
  }
  
  ## add to list of rasters
  if (i == 1) {
    rasterized_rrs <- rr_raster
    og_rasterized_rrs <- rr_og
  }
  else {
    rasterized_rrs <- addLayer(rasterized_rrs, rr_raster, updatevalue = NA)
    og_rasterized_rrs <- addLayer(og_rasterized_rrs, rr_og, updatevalue = NA)
  }
  
  print(paste("Finished range number:", i))
  i = i + 1
}

names(rasterized_rrs) <- realized_ranges$range_id
names(og_rasterized_rrs) <- realized_ranges$range_id

##saveRDS(rasterized_rrs, "data-processed/rasterized_rrs.rds")
rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
##saveRDS(og_rasterized_rrs, "data-processed/og_rasterized_rrs.rds")
og_rasterized_rrs <- readRDS("data-processed/og_rasterized_rrs.rds")

rasterized_rrs_nichemapr <- mask(rasterized_rrs, t_mask_nichemapr)
##saveRDS(rasterized_rrs_nichemapr, "data-processed/rasterized_rrs_nichemapr.rds")
rasterized_rrs_nichemapr <- readRDS("data-processed/rasterized_rrs_nichemapr.rds")

####################################################################################
#####                   CREATING POTENTIAL RANGE SHAPEFILES                   ######
####################################################################################
## partition thermal limits into upper and lower limits
upper_limits <- thermal_limits %>%
  filter(type == "max")
lower_limits <- thermal_limits %>%
  filter(type == "min") 

## keep limits of species with both an upper and lower thermal limits  
both_upper_all <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower_all <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

## keep limits of species with only an upper or lower thermal limit
only_upper_all <- upper_limits[!upper_limits$genus_species %in% lower_limits$genus_species,]
only_lower_all <- lower_limits[!lower_limits$genus_species %in% upper_limits$genus_species,]



#######################################################
#####            TERRESTRIAL  SPECIES:           ######
#######################################################
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Terrestrial")
both_lower <- filter(both_lower_all, realm == "Terrestrial")
only_upper <- filter(only_upper_all, realm == "Terrestrial")
only_lower <-  filter(only_lower_all, realm == "Terrestrial")

## prepare zoogeographic realms for use in restricting potential ranges:
zoo <- read_sf("data-raw/CMEC regions & realms/newRealms.shp")

# simplify the object to make it faster to use
zoo <- zoo %>% 
  st_simplify(dTolerance = 0.01) %>% 
  group_by(Realm) %>% 
  summarise()
plot(zoo)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
## filter out habitat with temperatures outside of the fundamental thermal niche 
combined_elev_x_dormancy <- filter_by_tolerance_elev_x_dormancy(both_upper = both_upper,
                                                                both_lower = both_lower,
                                                                raster_high = raster_terr_high, 
                                                                raster_low = raster_terr_low)
combined_Te_worst_tr <- filter_by_tolerance_Te(both_upper = both_upper,
                                               both_lower = both_lower,
                                               Te_sun = Te_sun_max, 
                                               Te_shade = Te_shade_min)

combined_Te_best_tr <- filter_by_tolerance_Te(both_upper = both_upper,
                                              both_lower = both_lower,
                                              Te_sun = Te_shade_max, 
                                              Te_shade = Te_shade_min)


combined_acclimated <- filter_by_tolerance_acclimated(both_upper = both_upper,
                                                      both_lower = both_lower,
                                                      raster_high = Te_shade_max, 
                                                      raster_low = Te_shade_min,
                                                      acc_data = acc_data, 
                                                      realm = "Terrestrial")

# combined_reg <- filter_by_tolerance(both_upper = both_upper, 
#                                     both_lower = both_lower, 
#                                     raster_high = raster_terr_high, 
#                                     raster_low = raster_terr_low)
# 
# combined_dormancy <- filter_by_tolerance_dormancy(both_upper = both_upper,
#                                                   both_lower = both_lower, 
#                                                   raster_high = raster_terr_high, 
#                                                   raster_low = raster_terr_low)
# 
# combined_elev <- filter_by_tolerance_elev(both_upper = both_upper, 
#                                           both_lower = both_lower, 
#                                           raster_high = raster_terr_high, 
#                                           raster_low = raster_terr_low)

## restrict range to:
## - zoogeographic realms that overlap species' realized range
## - habitat within the species elevational range 
prs_terrestrial_elev_x_dormancy <- create_potential_ranges(clumped_temps = zoo, 
                                                           realized_ranges = realized_ranges,
                                                           combined = combined_elev_x_dormancy,
                                                           type = 'elev_x_dormancy')
prs_terrestrial_Te_worst_tr <- create_potential_ranges(clumped_temps = zoo, 
                                                       realized_ranges = realized_ranges,
                                                       combined = combined_Te_worst_tr,
                                                       type = 'Te_worst')
prs_terrestrial_Te_best_tr <- create_potential_ranges(clumped_temps = zoo, 
                                                      realized_ranges = realized_ranges,
                                                      combined = combined_Te_best_tr,
                                                      type = 'Te_best')
prs_terrestrial_acclimated <- create_potential_ranges(clumped_temps = zoo, 
                                                      realized_ranges = realized_ranges,
                                                      combined = combined_acclimated,
                                                      type = 'acclimated')

# prs_terrestrial <- create_potential_ranges(clumped_temps = zoo, 
#                                            realized_ranges = realized_ranges, 
#                                            combined = combined_reg, 
#                                            type = 'reg')
# prs_terrestrial_dormancy <- create_potential_ranges(clumped_temps = zoo, 
#                                            realized_ranges = realized_ranges,
#       f                                     combined = combined_dormancy,
#                                            type = 'dormancy')
# prs_terrestrial_elev <- create_potential_ranges(clumped_temps = zoo, 
#                                                 realized_ranges = realized_ranges, 
#                                                 combined = combined_elev,
#                                                 type = 'elevation')


saveRDS(prs_terrestrial_elev_x_dormancy, "data-processed/prs_terrestrial_elev_x_dormancy.rds")
saveRDS(prs_terrestrial_Te_worst_tr, "data-processed/prs_terrestrial_Te_worst_tr.rds")
saveRDS(prs_terrestrial_Te_best_tr, "data-processed/prs_terrestrial_Te_best_tr.rds")
saveRDS(prs_terrestrial_acclimated, "data-processed/prs_terrestrial_acclimated.rds")
# saveRDS(prs_terrestrial, "data-processed/prs_terrestrial.rds")
# saveRDS(prs_terrestrial_dormancy, "data-processed/prs_terrestrial_dormancy_6mo.rds")
# saveRDS(prs_terrestrial_elev, "data-processed/prs_terrestrial_elev.rds")

prs_terrestrial_elev_x_dormancy <- readRDS("data-processed/prs_terrestrial_elev_x_dormancy.rds")
prs_terrestrial_Te_worst_tr <- readRDS("data-processed/prs_terrestrial_Te_worst_tr.rds")
prs_terrestrial_Te_best_tr <- readRDS("data-processed/prs_terrestrial_Te_best_tr.rds")
prs_terrestrial_acclimated <- readRDS("data-processed/prs_terrestrial_acclimated.rds")
# prs_terrestrial <- readRDS("data-processed/prs_terrestrial.rds")
# prs_terrestrial_dormancy <- readRDS("data-processed/prs_terrestrial_dormancy_6mo.rds")
# prs_terrestrial_elev <- readRDS("data-processed/prs_terrestrial_elev.rds")


## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_elev_x_dormancy_terr <- extract_thermal_niche(type = "elev_x_dormancy", 
                                                    prs = prs_terrestrial_elev_x_dormancy,
                                                    realm = 'terrestrial')
niche_Te_worst_terr <- extract_thermal_niche(type = "Te_worst", 
                                             prs = prs_terrestrial_Te_worst_tr,
                                             realm = 'terrestrial')
niche_Te_best_terr <- extract_thermal_niche(type = "Te_best", 
                                            prs = prs_terrestrial_Te_best_tr,
                                            realm = 'terrestrial')
# niche_reg_terr <- extract_thermal_niche(type = "reg", prs = prs_terrestrial, realm = 'terrestrial')
# niche_dormancy_terr <- extract_thermal_niche(type = "dormancy", prs = prs_terrestrial_dormancy, 
#                                              realm = 'terrestrial')
# niche_elev_terr <- extract_thermal_niche(type = "elevation", prs = prs_terrestrial_elev, realm = 'terrestrial')

## save:
write.csv(niche_elev_x_dormancy_terr, "data-processed/thermal-niche/niche_elev_x_dormancy_terr.csv", 
          row.names = FALSE)
write.csv(niche_Te_best_terr, "data-processed/thermal-niche/niche_Te_best_terr.csv", 
          row.names = FALSE)
write.csv(niche_Te_worst_terr, "data-processed/thermal-niche/niche_Te_worst_terr.csv", 
          row.names = FALSE)
# write.csv(niche_reg_terr, "data-processed/thermal-niche/niche_reg_terr.csv", row.names = FALSE)
# write.csv(niche_dormancy_terr, "data-processed/thermal-niche/niche_dormancy_terr.csv", row.names = FALSE)
# write.csv(niche_elev_terr, "data-processed/thermal-niche/niche_elev_terr.csv", row.names = FALSE)

niche_elev_x_dormancy_terr <- read.csv("data-processed/thermal-niche/niche_elev_x_dormancy_terr.csv")
niche_Te_terr <- read.csv("data-processed/thermal-niche/niche_Te_terr.csv")
# niche_reg_terr <- read.csv("data-processed/thermal-niche/niche_reg_terr.csv")
# niche_dormancy_terr <- read.csv("data-processed/thermal-niche/niche_dormancy_terr.csv")
# niche_elev_terr <- read.csv("data-processed/thermal-niche/niche_elev_terr.csv")



## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################
## filter out habitat with temperatures that exceed the species' single fundamental thermal niche limit
# combined_upper_e_x_d <- filter_by_tolerance_one_limit_e_x_d(lims = only_upper, 
#                                                             type = "upper", 
#                                                             raster_high = raster_terr_high, 
#                                                             raster_low = raster_terr_low)
# 
# combined_lower_e_x_d <- filter_by_tolerance_one_limit_e_x_d(lims = only_lower, 
#                                                             type = "lower",
#                                                             raster_high = raster_terr_high, 
#                                                             raster_low = raster_terr_low)

## restrict range to:
## - zoogeographic realms that overlap species' realized range
## - habitat within the species elevational range 
# prs_terrestrial_upper_e_x_d <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                            realized_ranges = realized_ranges, 
#                                            combined = combined_upper_e_x_d, 
#                                            type = 'elev_x_dormancy')
# 
# prs_terrestrial_lower_e_x_d <- create_potential_ranges_one_limit(clumped_temps = clumped_temps, 
#                                                        realized_ranges = realized_ranges, 
#                                                        combined = combined_lower_e_x_d, 
#                                                        type = 'elev_x_dormancy')


# ##saveRDS(prs_terrestrial_upper_e_x_d, "data-processed/prs_terrestrial_upper_e_x_d.rds")
# prs_terrestrial_upper_e_x_d <- readRDS("data-processed/prs_terrestrial_upper_e_x_d.rds")
# ##saveRDS(prs_terrestrial_lower_e_x_d, "data-processed/prs_terrestrial_lower_e_x_d.rds")
# prs_terrestrial_lower_e_x_d <- readRDS("data-processed/prs_terrestrial_lower_e_x_d.rds")


#######################################################
#####                MARINE SPECIES:             ######
#######################################################
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Marine")
both_lower <- filter(both_lower_all, realm == "Marine")
only_upper <- filter(only_upper_all, realm == "Marine")
only_lower <-  filter(only_lower_all, realm == "Marine")

## clump contiguous habitat together:
clumped_temps <- raster_marine_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
## filter out habitat with temperatures outside of the fundamental thermal niche 
combined_reg <- filter_by_tolerance(both_upper = both_upper,
                                    both_lower = both_lower,
                                    raster_high = raster_marine_high, 
                                    raster_low = raster_marine_low)
combined_bo <- filter_by_tolerance_biooracle(both_upper = both_upper,
                                             both_lower = both_lower,
                                             raster_high = raster_marine_high, 
                                             raster_low = raster_marine_low)

combined_acclimated <- filter_by_tolerance_acclimated(both_upper = both_upper,
                                                      both_lower = both_lower,
                                                      raster_high = raster_marine_high, 
                                                      raster_low = raster_marine_low,
                                                      acc_data = acc_data, 
                                                      realm = "Marine")


## restrict range to:
## - contiguous marine habitat that overlaps species' realized range
## - habitat within the species' depth distribution
prs_marine <- create_potential_ranges(clumped_temps = clumped_temps, 
                                      realized_ranges = realized_ranges, combined = combined_reg,
                                      type = 'reg')
prs_marine_bo <- create_potential_ranges(clumped_temps = clumped_temps, 
                                         realized_ranges = realized_ranges, combined = combined_bo,
                                         type = 'bio-oracle')
prs_marine_acclimated <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                 realized_ranges = realized_ranges, combined = combined_acclimated,
                                                 type = 'acclimated')

# saveRDS(prs_marine, "data-processed/prs_marine.rds")
# saveRDS(prs_marine_bo, "data-processed/prs_marine_bo.rds")
# saveRDS(prs_marine_acclimated, "data-processed/prs_marine_acclimated.rds")
prs_marine <- readRDS("data-processed/prs_marine.rds")
prs_marine_bo <- readRDS("data-processed/prs_marine_bo.rds")
prs_marine_acclimated <- readRDS("data-processed/prs_marine_acclimated.rds")

## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_reg_marine <- extract_thermal_niche(type = "reg", prs_marine, realm = 'marine')
niche_bo_marine <- extract_thermal_niche(type = "bio-oracle", prs_marine_bo, realm = 'marine')

write.csv(niche_reg_marine, "data-processed/thermal-niche/niche_reg_marine.csv", row.names = FALSE)
write.csv(niche_bo_marine, "data-processed/thermal-niche/niche_bo_marine.csv", row.names = FALSE)

niche_reg_marine <- read.csv("data-processed/thermal-niche/niche_reg_marine.csv")
niche_bo_marine <- read.csv("data-processed/thermal-niche/niche_bo_marine.csv")

## SPECIES WITH ONLY ONE THERMAL LIMIT:
#######################################


######################################################
#####            INTERTIDAL SPECIES:            ######
######################################################
## get rid of thermal limits of species in other realms:
both_upper <- filter(both_upper_all, realm == "Intertidal")
both_lower <- filter(both_lower_all, realm == "Intertidal")
only_upper <- filter(only_upper_all, realm == "Intertidal")
only_lower <-  filter(only_lower_all, realm == "Intertidal")

## clump contiguous habitat together:
clumped_temps <- raster_intertidal_high[[1]] %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps$geometry)

## SPECIES WITH BOTH THERMAL LIMITS:
####################################
## filter out habitat with temperatures outside of the fundamental thermal niche 
combined_reg <- filter_by_tolerance(both_upper = both_upper,
                                    both_lower = both_lower,
                                    raster_high = raster_intertidal_high,
                                    raster_low = raster_intertidal_low)
combined_bo <- filter_by_tolerance_biooracle(both_upper = both_upper,
                                             both_lower = both_lower,
                                             raster_high = raster_intertidal_high, 
                                             raster_low = raster_intertidal_low)

combined_acclimated <- filter_by_tolerance_acclimated(both_upper = both_upper,
                                                      both_lower = both_lower,
                                                      raster_high = raster_intertidal_high, 
                                                      raster_low = raster_intertidal_low,
                                                      acc_data = acc_data, 
                                                      realm = "Intertidal")


## restrict range to:
## - contiguous intertidal habitat (continental shelf beside a land mass) that overlaps species' realized range
## - habitat within the species' depth distribution
prs_intertidal <- create_potential_ranges(clumped_temps = clumped_temps, 
                                          realized_ranges = realized_ranges, combined = combined_reg,
                                          type = 'reg')
prs_intertidal_bo <- create_potential_ranges(clumped_temps = clumped_temps, 
                                             realized_ranges = realized_ranges, combined = combined_bo,
                                             type = 'bio-oracle')
prs_intertidal_acclimated <- create_potential_ranges(clumped_temps = clumped_temps, 
                                                     realized_ranges = realized_ranges, combined = combined_acclimated,
                                                     type = 'acclimated')

# saveRDS(prs_intertidal, "data-processed/prs_intertidal.rds")
# saveRDS(prs_intertidal_bo, "data-processed/prs_intertidal_bo.rds")
# saveRDS(prs_intertidal_acclimated, "data-processed/prs_intertidal_acclimated.rds")
prs_intertidal <- readRDS("data-processed/prs_intertidal.rds")
prs_intertidal_bo <- readRDS("data-processed/prs_intertidal_bo.rds")
prs_intertidal_acclimated <- readRDS("data-processed/prs_intertidal_acclimated.rds")

## EXTRACTING TEMPS IN THERMAL NICHE
#######################################
niche_reg_intertidal <- extract_thermal_niche(type = "reg", prs_intertidal, realm = 'intertidal')
niche_bo_intertidal <- extract_thermal_niche(type = "bio-oracle", prs_intertidal_bo, realm = 'intertidal')

write.csv(niche_reg_intertidal, "data-processed/thermal-niche/niche_reg_intertidal.csv", row.names = FALSE)
write.csv(niche_bo_intertidal, "data-processed/thermal-niche/niche_bo_intertidal.csv", row.names = FALSE)

niche_reg_intertidal <- read.csv("data-processed/thermal-niche/niche_reg_intertidal.csv")
niche_bo_intertidal <- read.csv("data-processed/thermal-niche/niche_bo_intertidal.csv")

## SPECIES WITH ONLY ONE THERMAL LIMIT:


########################################################
###     COLLATING POTENTIAL RANGES ACROSS REALMS     ###
########################################################
potential_ranges_both_limits_elev_x_dormancy <- stack(prs_terrestrial_elev_x_dormancy[[1]],
                                                      prs_marine[[1]],
                                                      prs_intertidal[[1]])
potential_ranges_both_limits_Te_best_tr <- stack(prs_terrestrial_Te_best_tr[[1]],
                                                 prs_marine[[1]],
                                                 prs_intertidal[[1]])
potential_ranges_both_limits_Te_worst_tr <- stack(prs_terrestrial_Te_worst_tr[[1]],
                                                  prs_marine[[1]],
                                                  prs_intertidal[[1]])
potential_ranges_both_limits_Te_exposed <- stack(prs_terrestrial_Te_exposed[[1]],
                                                 prs_marine[[1]],
                                                 prs_intertidal[[1]])

potential_ranges_both_limits_acclimated <- stack(prs_terrestrial_acclimated[[1]],
                                                 prs_marine_acclimated[[1]],
                                                 prs_intertidal_acclimated[[1]])

# potential_ranges_both_limits_elev_x_dormancy_bo <- stack(prs_terrestrial_elev_x_dormancy[[1]],
#                                                          prs_marine_bo[[1]],
#                                                          prs_intertidal_bo[[1]])
# potential_ranges_both_limits <- stack(prs_terrestrial[[1]], prs_marine[[1]], prs_intertidal[[1]])
# potential_ranges_both_limits_dormancy <- stack(prs_terrestrial_dormancy[[1]], prs_marine[[1]],
#                                                prs_intertidal[[1]])
# potential_ranges_both_limits_elev <- stack(prs_terrestrial_elev[[1]], prs_marine[[1]], prs_intertidal[[1]])


saveRDS(potential_ranges_both_limits_elev_x_dormancy, "data-processed/potential_ranges_notcutatequator_elev_x_dormancy.rds")
saveRDS(potential_ranges_both_limits_Te, "data-processed/potential_ranges_notcutatequator_Te.rds")
# saveRDS(potential_ranges_both_limits_elev_x_dormancy_bo, "data-processed/potential_ranges_notcutatequator_elev_x_dormancy_bo.rds")
# saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges.rds")
# saveRDS(potential_ranges_both_limits, "data-processed/potential_ranges_notcutatequator.rds")
# saveRDS(potential_ranges_both_limits_dormancy, "data-processed/potential_ranges_notcutatequator_dormancy.rds")
# saveRDS(potential_ranges_both_limits_elev, "data-processed/potential_ranges_notcutatequator_elev.rds")



## COLLATING TEMPS IN THERMAL NICHE ACROSS REALMS
##################################################
## do some rearranging to make niches across realms have the same number of types
niche_elev_x_dormancy_intertidal <- niche_reg_intertidal %>%
  mutate(type = ifelse(type == "realized_reg", "realized_elev_x_dormancy", 'potential_elev_x_dormancy'))
niche_Te_best_intertidal <- niche_reg_intertidal %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_best", 'potential_Te_best'))
niche_Te_worst_intertidal <- niche_reg_intertidal %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_worst", 'potential_Te_worst'))
# niche_elev_intertidal <- niche_reg_intertidal %>%
#   mutate(type = ifelse(type == "realized_reg", "realized_elevation", 'potential_elevation'))
# niche_dormancy_intertidal <- niche_reg_intertidal %>%
#   mutate(type = ifelse(type == "realized_reg", "realized_dormancy", 'potential_dormancy'))

niche_elev_x_dormancy_marine <- niche_reg_marine %>%
  mutate(type = ifelse(type == "realized_reg", "realized_elev_x_dormancy", 'potential_elev_x_dormancy'))
niche_Te_best_marine <- niche_reg_marine %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_best", 'potential_Te_best'))
niche_Te_worst_marine <- niche_reg_marine %>%
  mutate(type = ifelse(type == "realized_reg", "realized_Te_worst", 'potential_Te_worst'))
# niche_elev_marine <- niche_reg_marine %>%
#   mutate(type = ifelse(type == "realized_reg", "realized_elevation", 'potential_elevation'))
# niche_dormancy_marine <- niche_reg_marine %>%
#   mutate(type = ifelse(type == "realized_reg", "realized_dormancy", 'potential_dormancy'))


niche_elev_x_dormancy <- rbind(niche_elev_x_dormancy_terr, niche_elev_x_dormancy_marine, 
                               niche_elev_x_dormancy_intertidal)
niche_Te_best <- rbind(niche_Te_best_terr, niche_Te_best_marine, 
                       niche_Te_best_intertidal)
niche_Te_worst <- rbind(niche_Te_worst_terr, niche_Te_worst_marine, 
                        niche_Te_worst_intertidal)
# niche_reg <- rbind(niche_reg_terr, niche_reg_marine, niche_reg_intertidal)
# niche_elev <- rbind(niche_elev_terr, niche_elev_marine, niche_elev_intertidal)
# niche_dormancy <- rbind(niche_dormancy_terr, niche_dormancy_marine, niche_dormancy_intertidal)

write.csv(niche_elev_x_dormancy, "data-processed/thermal-niche/niche_elev_x_dormancy.csv", row.names = FALSE)
write.csv(niche_Te_best, "data-processed/thermal-niche/niche_Te_best.csv", row.names = FALSE)
write.csv(niche_Te_worst, "data-processed/thermal-niche/niche_Te_worst.csv", row.names = FALSE)
# write.csv(niche_reg, "data-processed/thermal-niche/niche_reg.csv", row.names = FALSE)
# write.csv(niche_elev, "data-processed/thermal-niche/niche_elev.csv", row.names = FALSE)
# write.csv(niche_dormancy, "data-processed/thermal-niche/niche_dormancy.csv", row.names = FALSE)



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

filter_by_tolerance_Te <- function(both_upper, 
                                   both_lower, 
                                   Te_sun, 
                                   Te_shade) {
  ## create an individual raster layer of difference between thermal limit and Te
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  ## subset to species that we have a Te for:
  upper <- filter(both_upper, both_upper$genus_species %in% names(Te_sun))
  lower <- filter(both_lower, both_lower$genus_species %in% names(Te_sun))
  
  species = 1
  while (species < nrow(upper) + 1) {
    
    ## get species name and extract its Te layers:
    sp <- upper$genus_species[species] 
    
    h <- Te_sun[[which(names(Te_sun) == sp)]]
    l <- Te_shade[[which(names(Te_shade) == sp)]]
    
    rr_high <- addLayer(rr_high, h) 
    rr_low <- addLayer(rr_low, l) 
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - hottest Te > 0 and where Tmin - coldest Te< 0)
    below_ctmax <- upper$thermal_limit[species] - h
    h[below_ctmax < 0] <- NA
    above_ctmin <- lower$thermal_limit[species] - l
    l[above_ctmin > 0] <- NA
    
    high <- addLayer(high, h) 
    low <- addLayer(low, l) 
    
    species = species + 1
  }
  
  names(high) <- upper$genus_species
  names(low) <- lower$genus_species
  
  ## combine to find cells where Te_sun is less than CTmax AND Te_shade is greater than CTmin 
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
      
      tracker <- append(tracker, "both")
    }
    else if (species_traits$cold_season_dormancy_ == "Yes") {
      h <- raster_high$cold_dormancy_6mo
      l <-  raster_low$cold_dormancy_6mo
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      tracker <- append(tracker, "cold")
    }
    else if (species_traits$hot_season_dormancy_ == "Yes") {
      h <- raster_high$hot_dormancy_6mo
      l <- raster_low$hot_dormancy_6mo
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      tracker <- append(tracker, rep("hot"))
    }
    else {
      h <- raster_high$seasonal_high_temp
      l <- raster_low$seasonal_low_temp
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      tracker <- append(tracker, "neither")
    }
    
    ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
    below_ctmax <- both_upper$thermal_limit[species] - h
    h[below_ctmax < 0] <- NA
    above_ctmin <- both_lower$thermal_limit[species] - l
    l[above_ctmin > 0] <- NA
    
    high <- addLayer(high, h) 
    low <- addLayer(low, l) 
    
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
    
    if (species_traits$cold_season_dormancy_ == "Yes" & 
        species_traits$hot_season_dormancy_ == "Yes") {
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
      
      tracker <- append(tracker, "hot")
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
      
      tracker <- append(tracker, "neither")
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

filter_by_tolerance_biooracle <- function(both_upper, 
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
    l <- raster_low$ltmin_at_max_depth  # use bio-oracle temp set 
    
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

filter_by_tolerance_one_limit <- function(lims, 
                                          raster_high, 
                                          raster_low,
                                          type) {
  
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  temp <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_temp <- temp
  
  species = 1
  while (species < nrow(lims) + 1) {
    
    if (type == 'upper') {
      t <- raster_high$seasonal_high_temp 
    }
    else if (type == 'lower') {
      t <- raster_low$seasonal_low_temp 
    }
    
    rr_temp <- addLayer(rr_temp, t) 
    
    ## exclude raster cells above upper thermal tolerance limit or below lower thermal tolerance limit 
    if (type == 'upper') {
      below_ctmax <- lims$thermal_limit[species] - t
      t[below_ctmax < 0] <- NA
    }
    else if (type == 'lower') {
      above_ctmin <- lims$thermal_limit[species] - t
      t[above_ctmin > 0] <- NA
    }
    
    temp <- addLayer(temp, t) 
    
    species = species + 1
  }
  
  names(temp) <- paste(lims$Genus, lims$Species, sep = "_")
  
  return (list(temp, rr_temp))
}

filter_by_tolerance_one_limit_e_x_d <- function(lims, 
                                                raster_high, 
                                                raster_low,
                                                type) {
  
  ## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
  temp <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_temp <- temp
  
  species = 1
  names <- c()
  tracker <- c()
  while (species < nrow(lims) + 1) {
    species_traits <- traits[which(traits$genus_species == lims$genus_species[species]),]
    
    if (species_traits$cold_season_dormancy_ == "Yes" && type == "lower") {
      t <-  raster_low$low_elev_x_cold_dormancy
      
      rr_temp <- addLayer(rr_temp, t) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      above_ctmin <- lims$thermal_limit[species] - t
      t[above_ctmin > 0] <- NA
      
      temp <- addLayer(temp, t) 
      
      tracker <- append(tracker, "cold")
    }
    else if (species_traits$hot_season_dormancy_ == "Yes" && type == "upper") {
      t <- raster_high$hot_dormancy_6mo
      
      rr_temp <- addLayer(rr_temp, t) 
      
      ## exclude raster cells outside of the thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      below_ctmax <- lims$thermal_limit[species] - t
      t[below_ctmax < 0] <- NA
      
      temp <- addLayer(temp, t) 
      
      tracker <- append(tracker, "hot")
    }
    else {
      rr_temp <- addLayer(rr_temp, t) 
      
      if (type == 'upper') {
        t <- raster_high$seasonal_high_temp
        below_ctmax <- lims$thermal_limit[species] - t
        t[below_ctmax < 0] <- NA
      }
      else if (type == 'lower') {
        t <- raster_low$low_at_min_elev
        above_ctmin <- lims$thermal_limit[species] - t
        t[above_ctmin > 0] <- NA
      }
      
      temp <- addLayer(temp, t) 
      
      tracker <- append(tracker, "neither")
    }
    
    species = species + 1
  }
  
  names(temp) <- paste(lims$Genus, lims$Species, tracker, 'elev', sep = "_")
  names(rr_temp) <- names(temp)
  
  return (list(temp, rr_temp))
}

filter_by_tolerance_acclimated <- function(both_upper,
                                           both_lower,
                                           raster_high, 
                                           raster_low,
                                           acc_data, realm) {
  
  ## create an individual raster layer of difference between thermal limit and Te
  high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  low <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  
  ## create set of temps to use to extract realized thermal niches later 
  rr_high <- high
  rr_low <- low
  
  tlims_upper_all <- high
  tlims_lower_all <- low
  
  ## read in acclimation temperatures:
  if (realm == "Terrestrial") {
    acc_cold <- stack("data-processed/terr_acc_cold.grd")
    acc_hot <- stack("data-processed/terr_acc_hot.grd")
    
    ## subset to species that we have a Te for:
    upper <- filter(both_upper, both_upper$genus_species %in% names(raster_high))
    lower <- filter(both_lower, both_lower$genus_species %in% names(raster_low))
    
    species = 1
    while (species < nrow(upper) + 1) {
      
      ## get species name and extract its Te layers:
      sp <- upper$genus_species[species] 
      
      h <- raster_high[[which(names(raster_high) == sp)]]
      l <- raster_low[[which(names(raster_low) == sp)]]
      
      rr_high <- addLayer(rr_high, h) 
      rr_low <- addLayer(rr_low, l) 
      
      ## create raster layers of thermal limits if species was acclimated to the cell
      tlims_upper <- acc_hot[[which(names(acc_hot) == sp)]]
      tlims_lower <- acc_cold[[which(names(acc_cold) == sp)]]
      
      equ <- acc_data[which(acc_data$genus_species == sp),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      if (nrow(hot_equ) != 0) {
        if (!is.na(hot_equ$acclimation_temperature)){
          slope = hot_equ$ARR_equ_slope
          int = hot_equ$ARR_equ_int
          
          tlims_upper <- slope*tlims_upper + int
        }
        else {
          tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
        }
      }
      else {
        tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
      }
      if (nrow(cold_equ) != 0) {
        if (!is.na(cold_equ$acclimation_temperature)){
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims_lower <- slope*tlims_lower + int
        }
        else {
          tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
        }
      }
      else {
        tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[which(lower$genus_species == sp)]
      }
      
      ## exclude raster cells outside of the thermal tolerance (where acclimated Tmax - hottest Te > 0 and where acclimated Tmin - coldest Te< 0)
      h[tlims_upper < h]  <- NA
      l[tlims_lower > l] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      species = species + 1
    }
  }
  else {
    acc_cold <- raster_low[[4]]
    acc_hot <- raster_high[[4]]
    upper <- both_upper
    lower <- both_lower
    
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
        
        tracker <- append(tracker, "both")
      }
      else if (species_traits$cold_season_dormancy_ == "Yes") {
        h <- raster_high$cold_dormancy_6mo
        l <-  raster_low$cold_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
        
        tracker <- append(tracker, "cold")
      }
      else if (species_traits$hot_season_dormancy_ == "Yes") {
        h <- raster_high$hot_dormancy_6mo
        l <- raster_low$hot_dormancy_6mo
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
        
        tracker <- append(tracker, "hot")
      }
      else {
        h <- raster_high$seasonal_high_temp
        l <- raster_low$seasonal_low_temp
        
        rr_high <- addLayer(rr_high, h) 
        rr_low <- addLayer(rr_low, l) 
        
        tracker <- append(tracker, "neither")
      }
      
      ## create raster layer of thermal limits if species was acclimated to the cell
      equ <- acc_data[which(acc_data$genus_species == species_traits$genus_species),]
      hot_equ <- equ[which(equ$type == "max"),]
      cold_equ <- equ[which(equ$type == "min"),]
      
      if (nrow(hot_equ) != 0) {
        
        slope = hot_equ$ARR_equ_slope
        int = hot_equ$ARR_equ_int
        
        tlims_upper <- slope*acc_hot + int
        
        if (nrow(cold_equ) != 0) {
          slope = cold_equ$ARR_equ_slope
          int = cold_equ$ARR_equ_int
          
          tlims_lower <- slope*acc_cold + int
        }
      }
      else if (nrow(cold_equ) != 0) {
        slope = cold_equ$ARR_equ_slope
        int = cold_equ$ARR_equ_int
        
        tlims_lower <- slope*acc_cold + int
        
        tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
      }
      else {
        tlims_upper[!is.na(tlims_upper)] <- upper$thermal_limit[species]
        tlims_lower[!is.na(tlims_lower)] <- lower$thermal_limit[species]
      }
      
      ## exclude raster cells outside of the acclimated thermal tolerance (where Tmax - seasonal_high > 0 and where Tmin - seasonal_low< 0)
      h[tlims_upper < h]  <- NA
      l[tlims_lower > l] <- NA
      
      high <- addLayer(high, h) 
      low <- addLayer(low, l) 
      
      ## get rr and crop:
      rr <- rasterized_rrs[[which(str_split_fixed(names(rasterized_rrs), "_", n=2)[,1] 
                                  == str_replace_all(sp, "_", "."))]]
      tlims_upper <- mask(tlims_upper, rr)
      tlims_lower <- mask(tlims_lower, rr)
      
      tlims_upper_all <- addLayer(tlims_upper_all, tlims_upper) 
      tlims_lower_all <- addLayer(tlims_lower_all, tlims_lower) 
      
      species = species + 1
    }
    
    names(high) <- paste(upper$Genus, upper$Species, tracker, sep = "_")
    names(low) <- paste(upper$Genus, upper$Species, tracker, sep = "_")
    names(tlims_upper_all) <- paste(upper$Genus, upper$Species, sep = "_")
    names(tlims_lower_all) <- paste(upper$Genus, upper$Species, sep = "_")
    ##plot(high)
    ##plot(low)
  }
  
  
  ## combine to find cells where Te_sun is less than CTmax AND Te_shade is greater than CTmin 
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
  
  return (list(combined_high, combined_low, rr_high, rr_low, tlims_upper_all, tlims_lower_all))
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
  
  if (type == "acclimated") {
    tlims_upper_all <- combined[[5]]
    tlims_lower_all <- combined[[6]]
    tlims_upper <- rr_high_new
    tlims_lower <- rr_high_new
  }
  
  pr_restricted_bool <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  pr_restricted_high <- pr_restricted_bool
  pr_restricted_low <- pr_restricted_bool
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1, vals = 1)
  
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  names <- c()
  i = 1
  while (i < length(names(combined_high)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- combined_high[[i]] 
    
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
    else if (type == "acclimated") {
      if (str_detect(names(potential_raster), "neither")) {
        split <- str_split_fixed(names(potential_raster), "_", n = 3)
        species <- paste(split[1,1], split[1,2], sep = " ")
      }
      else {
        species <- str_replace_all(names(potential_raster), "_", " ")
      }
    }
    else {
      species <- str_replace_all(names(potential_raster), "_", " ")
    }
    
    realized_range <- realized_ranges[which(as.character(realized_ranges$species) %in% species),] %>%
      st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0") 
    
    ## get traits 
    sp_traits <- traits[which(str_replace_all(traits$genus_species, "_", " ") == species), ]
    
    ## if terrestrial, restrict by elevational range limits:
    if (sp_traits$Realm == "Terrestrial") {
      ## get elevational range limits 
      elev_max <- sp_traits$upper_elevation_limit_meters
      elev_min <- sp_traits$lower_elevation_limit_meters
      
      ## if both limits are not NA, restrict by both limits
      if (!is.na(elev_min) & !is.na(elev_max)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with max elev less than species min elevation
        ## exclude cells with min elev greater than species max elevation
      }
      else if (!is.na(elev_min)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        ## exclude cells with max elev less than species min elevation
      }
      else if (!is.na(elev_max)) {
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with min elev greater than species max elevation
      }
    }
    ## if marine, restrict by depth distribution limits:
    else if (sp_traits$Realm %in% c("Marine", "Intertidal")) {
      ## get depth distribution limits 
      depth_upper <- sp_traits$upper_depth_limit
      depth_lower <- sp_traits$lower_depth_limit
      
      ## if both limits are not NA, restrict by both depth limits
      if (!is.na(depth_lower) & !is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with deepest depth above species upper limit  
        ## exclude cells with shallowest depth less than species lower limit 
      }
      else if (!is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        ## exclude cells with deepest depth above species upper limit  
      }
      else if (!is.na(depth_lower)) {
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with shallowest depth less than species lower limit 
      }
      ## if lower depth limit is unknown, restrict to continental shelf if it is a coastal spp
      else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
        is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
        
        if(is_coastal) {
          potential_raster[is.na(shelf_mask)] <- NA
        }
      }
    }
    
    ## if species has multiple realized ranges, loop through them:
    num = 1
    while (num < nrow(realized_range)+1) {
      
      ## find areas of habitat that realized range intersects
      intersects <- st_intersects(clumped_temps, realized_range[num,], sparse = FALSE)[,]
      intersects <- filter(clumped_temps, intersects == TRUE) 
      intersects$geometry <- st_union(intersects) # combine all habitat into one for rasterizing
      
      ## rasterize them:
      r_intersects <- shp2rast(r, as_Spatial(intersects[1,]))
      
      if (type == "acclimated") {
        ## save acclimatized thermal limits across area allowed to be in potential range
        tlims_upper <- addLayer(tlims_upper, mask(tlims_upper_all[[i]], r_intersects))
        tlims_lower <- addLayer(tlims_lower, mask(tlims_lower_all[[i]], r_intersects))
      }
      
      pr_restricted <- mask(potential_raster, intersects) ## get rid of potential range outside of areas that realized range intersects
      pr_restricted[!is.na(pr_restricted)] <- 1 ## make it boolean
      
      ## create boolean version of potential range 
      pr_restricted_bool <- addLayer(pr_restricted_bool, pr_restricted, updatevalue = NA)
      
      ## create version of prs with temperature values for extracting thermal niche later 
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
  names(tlims_upper) <- names(tlims_lower) <- names
  #plot(pr_restricted_bool)
  
  return(list(pr_restricted_bool, pr_restricted_high, pr_restricted_low, rr_high_new, rr_low_new,
              tlims_upper, tlims_lower))
}

create_potential_ranges_one_limit <- function (clumped_temps, 
                                               combined, 
                                               realized_ranges, 
                                               type) {
  
  combined_temp <- combined[[1]]
  rr_temp <- combined[[2]]
  
  rr_new <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  pr_restricted_bool <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  pr_restricted <- pr_restricted_bool
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  names <- c()
  i = 1
  while (i < length(names(combined_temp)) + 1) {
    ## get unrestricted potential range 
    potential_raster <- combined_temp[[i]] 
    
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
    
    ## get traits 
    sp_traits <- traits[which(traits$genus_species == species), ]
    
    ## if terrestrial, restrict by elevational range limits:
    if (sp_traits$Realm == "Terrestrial") {
      elev_max <- sp_traits$upper_elevation_limit_meters
      elev_min <- sp_traits$lower_elevation_limit_meters
      
      ## if not NA, restrict by elevation
      if (!is.na(elev_max)) {
        potential_raster[elevs$elev_max < elev_min] <- NA 
        ## exclude cells with max elev less than species min
      }
      else if (!is.na(elev_min)) {
        potential_raster[elevs$elev_min > elev_max] <- NA 
        ## exclude cells with min elev greater than species max 
      }
    }
    ## if marine, restrict by depth distribution limits:
    else if (sp_traits$Realm == "Marine") {
      depth_upper <- sp_traits$upper_depth_limit
      depth_lower <- sp_traits$lower_depth_limit
      
      ## if not NA, restrict by depth
      if (!is.na(depth_upper)) {
        potential_raster[depths$depth_min > depth_upper] <- NA 
        ## exclude cells with deepest depth above species upper limit  
      }
      else if (!is.na(depth_lower)) {
        potential_raster[depths$depth_max < depth_lower] <- NA 
        ## exclude cells with shallowest depth less than species lower limit 
      }
      ## if max depth is na, restrict to continental shelf if it is a coastal spp
      else if (is.na(depth_lower) & !is.na(sp_traits$coastal_or_oceanic)) {
        is_coastal <- sp_traits$coastal_or_oceanic == "coastal"
        
        if(is_coastal) {
          potential_raster[is.na(shelf_mask)] <- NA
        }
        
      }
    }
    
    ## if species has multiple realized ranges, loop through them:
    num = 1
    r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
                res = 1, vals = 1)
    while (num < nrow(realized_range)+1) {
      
      intersects <- st_intersects(clumped_temps, realized_range[num,], sparse = FALSE)[,]
      intersects <- filter(clumped_temps, intersects == TRUE) ## find areas of habitat that realized range intersects
      intersects$geometry <- st_union(intersects) # combine all habitat into one for rasterizing
      
      ## rasterize them:
      r_intersects <- shp2rast(r, as_Spatial(intersects[1,]))
      
      pr_restricted <- mask(potential_raster, intersects) ## get rid of potential range outside of areas that realized range intersects
      pr_restricted[!is.na(pr_restricted)] <- 1 ## make it boolean
      
      ## create boolean version of potential range 
      pr_restricted_bool <- addLayer(pr_restricted_bool, pr_restricted, updatevalue = NA)
      
      ## create version of prs with temperatures as raster values for thermal niche metrics 
      pr_restricted <- addLayer(pr_restricted, 
                                mask(combined_temp[[i]], pr_restricted, updatevalue = NA), 
                                updatevalue = NA)
      
      ## add layer to new rr:
      rr_new <- addLayer(rr_new, rr_temp[[i]])
      rr_new <- addLayer(rr_new, rr_temp[[i]])
      
      names <- append(names, paste(names(potential_raster), realized_range$source[num], sep = "_"))
      
      num = num + 1
    }
    
    i = i + 1
  }
  names(pr_restricted_bool) <- names
  names(pr_restricted) <- names
  names(rr_new) <- names
  #plot(pr_restricted_bool)
  
  return(list(pr_restricted_bool, pr_restricted, rr_new))
}


extract_thermal_niche <- function(type, prs, realm) {
  
  rrs <- readRDS("data-processed/rasterized_rrs.rds")
  names(rrs) <- str_replace_all(names(rrs), '\\.', '_')
  
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
    range_id <- str_replace_all(range_id, "\\.", "_")
    
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
plot_ranges_overlap  <- function (folder, potential_ranges, type) {
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  
  land <- as.data.frame(raster("data-processed/raster_terr_mask.grd"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("data-processed/raster_marine_mask.grd"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("data-processed/raster_intertidal_mask.grd"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  if (type %in% c("Te_best_tr", "Te_worst_tr", "acclimated")) {
    land <- as.data.frame(raster("data-processed/raster_terr_mask_nichemapr.grd"), xy=TRUE)
    colnames(land)[1:3] <- c("longitude", "latitude", "mask")
    
    rasterized_rrs <- readRDS("data-processed/rasterized_rrs_nichemapr.rds")
  }
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  prs <- as.data.frame(potential_ranges, xy=TRUE)
  colnames(prs)[1:2] <- c("longitude", "latitude")
  
  
  ## ggplots of all:
  i = 1
  while (i < ncol(prs) - 1) {
    range <- colnames(prs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("_", ".")
    
    if(type == "elev_x_dormancy") {
      split <- str_split_fixed(range, "\\_", n = 5)
      source <- split[1,5]
    } 
    else if (type == "dormancy") {
      split <- str_split_fixed(range, "\\_", n = 4)
      source <- split[1,4]
    }
    else if (type == "acclimated") {
      if (str_detect(range, "neither")) {
        split <- str_split_fixed(range, "_", n = 4)
        species <- paste(split[1,1], split[1,2], sep = ".")
        source <- split[1,4]
      }
      else {
        split <- str_split_fixed(range, "_",  n = 3)
        source <- split[1,3]
      }
    }
    else {
      split <- str_split_fixed(range, "\\_", n = 3)
      source <- split[1,3]
    }
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(traits$genus_species == paste(split[1,1], split[1,2], sep = "_"))]
    
    if(realm == "Terrestrial") {
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
    r$rr <- as.factor(r$rr)
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill = rr)) + 
      scale_fill_manual(values = "yellow", na.value = "transparent") +
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
    ggsave(r_gg, path = paste("figures/range-plots/", folder, sep = ""), 
           filename = paste(range, ".png", sep = ""), 
           height = 6, width = 10, units = "in", device = "png")
    
    i = i + 1
  }
  
}

## plot the ranges:
plot_ranges_overlap(folder = "elev-x-dormancy", 
                    potential_ranges = potential_ranges_both_limits_elev_x_dormancy,
                    type = "elev_x_dormancy")
plot_ranges_overlap(folder = "Te-best-tr", potential_ranges = potential_ranges_both_limits_Te_best_tr,
                    type = "Te_best_tr")
plot_ranges_overlap(folder = "Te-worst-tr", potential_ranges = potential_ranges_both_limits_Te_worst_tr,
                    type = "Te_worst_tr")
plot_ranges_overlap(folder = "Te-exposed", potential_ranges = potential_ranges_both_limits_Te_exposed,           type = "Te_worst_tr")
plot_ranges_overlap(folder = "acclimated", potential_ranges = potential_ranges_both_limits_acclimated,           type = "acclimated")

# plot_ranges_overlap(folder = "no-modification", potential_ranges = potential_ranges_both_limits, type = "reg")
# plot_ranges_overlap(folder = "dormancy", potential_ranges = potential_ranges_both_limits_dormancy, type = "dormancy")
# plot_ranges_overlap(folder = "elev", potential_ranges = potential_ranges_both_limits_elev, type = "elev")
# plot_ranges_overlap(folder = "bio-oracle-elev-x-dormancy", potential_ranges = 
#                       potential_ranges_both_limits_elev_x_dormancy_bo,
#                     type = "elev_x_dormancy")


plot_rrs <- function () {
  rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
  og_rasterized_rrs <- readRDS("data-processed/og_rasterized_rrs.rds")
  
  rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
  colnames(rrs)[1:2] <- c("longitude", "latitude")
  
  og_rrs <- as.data.frame(og_rasterized_rrs, xy=TRUE)
  colnames(og_rrs)[1:2] <- c("longitude", "latitude")
  
  land <- as.data.frame(raster("data-processed/raster_terr_mask.grd"), xy=TRUE)
  colnames(land)[1:3] <- c("longitude", "latitude", "mask")
  ocean <- as.data.frame(raster("data-processed/raster_marine_mask.grd"), xy=TRUE)
  colnames(ocean)[1:3] <- c("longitude", "latitude", "mask")
  intertidal <- as.data.frame(raster("data-processed/raster_intertidal_mask.grd"), xy=TRUE)
  colnames(intertidal)[1:3] <- c("longitude", "latitude", "mask")
  
  ## ggplots of all:
  i = 1
  while (i < ncol(og_rrs) - 1) {
    range <- colnames(og_rrs)[i+2]
    
    ## get realized range:
    species <- range %>%
      str_replace_all("\\_", ".")
    
    split <- str_split_fixed(species, "\\.", n = 3)
    source <- split[1,3]
    
    ## get species name and realm
    species <- paste(split[1,1], split[1,2], sep = ".")
    realm <- traits$Realm[which(paste(traits$Genus, traits$Species, sep = "_") == 
                                  paste(split[1,1], split[1,2], sep = "_"))]
    
    if(realm == "Terrestrial" | realm == "Freshwater") {
      backdrop = land
    }
    else if (realm == "Marine") {
      backdrop = ocean
    }
    else {
      backdrop = intertidal
    }
    
    r <- rrs[,c(1:2, i+2)] 
    colnames(r)[3] <- "rr"
    r <- left_join(r, og_rrs[,c(1:2, i+2)]) 
    colnames(r)[4] <- "og_rr"
    
    r_gg <- r %>%
      ggplot(., aes(x = longitude, y = latitude)) +
      xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
      geom_raster(aes(fill=as.factor(og_rr))) + 
      scale_fill_manual(values = c("red"), aesthetics = 'fill', labels = ) +
      annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.5,
               fill = r$rr) +
      annotate(geom="raster", x=land$longitude, y=land$latitude, alpha=.1,
               fill = backdrop$mask) +
      labs(title = range,
           y = "Latitude",
           x = "Longitude") +
      scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
      scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
      theme(legend.position = "none")
    
    
    ## write to file:
    ggsave(r_gg, path = paste("figures/realized-range-plots/", sep = ""), 
           filename = paste(range, ".png", sep = ""), 
           height = 6, width = 10, units = "in", device = "png")
    
    
    i = i + 1
  }
  
}



## convincing myself Te ranges are mostly restricted by upper thermal limit 
## get a species' thermal tolerance 
spp1 <- thermal_limits %>%
  filter(genus_species == "Plestiodon_gilberti") 

ctmax <- spp1$thermal_limit[which(spp1$type == "max")]
ctmin <- spp1$thermal_limit[which(spp1$type == "min")]

## get its max Te sun and min Te shade rasters
sun <- Te_sun[[which(names(Te_sun) == "Plestiodon_gilberti")]]
shade <- Te_shade_min[[which(names(Te_shade_min) == "Plestiodon_gilberti")]]

## look at where max Te sun > upper thermal limit 
not_too_hot <- ctmax - sun
not_too_hot[not_too_hot < 0] <- -1
not_too_hot[not_too_hot >=0] <- 1
not_too_cold <- ctmin - shade
not_too_cold[not_too_cold >0] <- 2
not_too_cold[not_too_cold<=0] <- 1
not_too_cold[not_too_cold == 2] <- -1

plot(not_too_hot == -1)
plot(not_too_cold == -1)

plot(not_too_hot == 1 & not_too_cold == 1)


plot(sun < 41.5)

## look at where max Te shade < lower thermal limit 



acc_data %>%
  filter(duplicated(genus_species) & (realm %in% c("Marine", "Intertidal") | (acc_data$genus_species %in% names(Te_sun)))) %>%
  nrow(.) 



## plotting thermal lims:
tup <- data.frame(rasterToPoints(tlims_upper_terr))
colnames(tup) <- c("x","y",names(tlims_upper_terr))
tup <- janitor::remove_constant(tup)

tup <- tup %>%
  mutate(x = as.character(x), y = as.character(y))

tlow <- data.frame(rasterToPoints(tlims_lower_terr))
colnames(tlow) <- c("x","y",names(tlims_lower_terr))
tlow <- janitor::remove_constant(tlow)

tlow <- tlow %>%
  mutate(x = as.character(x), y = as.character(y))

mlow <- data.frame(rasterToPoints(tlims_lower_marine))
colnames(mlow) <- c("x","y",names(tlims_lower_marine))
mlow <- janitor::remove_constant(mlow)

mlow <- mlow %>%
  mutate(x = as.character(x), y = as.character(y))

mup <- data.frame(rasterToPoints(tlims_upper_marine))
colnames(mup) <- c("x","y",names(tlims_upper_marine))
mup <- janitor::remove_constant(mup)

mup <- mup %>%
  mutate(x = as.character(x), y = as.character(y))

iup <- data.frame(rasterToPoints(tlims_upper_int))
colnames(iup) <- c("x","y",names(tlims_upper_int))
iup <- janitor::remove_constant(iup)

iup <- iup %>%
  mutate(x = as.character(x), y = as.character(y))

ilow <- data.frame(rasterToPoints(tlims_lower_int))
colnames(ilow) <- c("x","y",names(tlims_lower_int))
ilow <- janitor::remove_constant(ilow)

ilow <- ilow %>%
  mutate(x = as.character(x), y = as.character(y))


tlims_upper_all <- left_join(tup, mup, iup, by = c("x", "y"))
tlims_lower_all <- left_join(tlow, mlow, ilow, by = c("x", "y"))

mean_data <- tlims_upper_all %>%
  summarise(across(1:ncol(.), ~ mean(.x, na.rm = TRUE))) %>%
  gather(key = "spp", value = "mean_acc_lim", c(3:ncol(.))) %>%
  left_join(.,filter(acc_data, type == "max"), by = c("spp" = "genus_species")) 

mean_data %>%
  ggplot(., aes(x = thermal_limit, y = mean_acc_lim)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Original upper thermal limit (C)", y = "Mean acclimatized upper thermal limit (C)")

mean_data <- tlims_lower_all %>%
  summarise(across(1:ncol(.), ~ mean(.x, na.rm = TRUE))) %>%
  gather(key = "spp", value = "mean_acc_lim", c(3:ncol(.))) %>%
  left_join(.,filter(acc_data, type == "min"), by = c("spp" = "genus_species")) 

mean_data %>%
  ggplot(., aes(x = thermal_limit, y = mean_acc_lim)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Original lower thermal limit (C)", y = "Mean acclimatized lower thermal limit (C)")


data <- data.frame(rasterToPoints(tlims_upper_all[[35]]))
colnames(data) <- c("x","y","acc_tlim")

ggplot(data, aes(x = acc_tlim)) + geom_histogram() + geom_vline(xintercept = 44.9, col = "red") +
  labs(x = "Acclimatized upper thermal limit (C)", y = "Frequency") + 
  annotate(geom = "text", x = 44.9 + 1.5, y = 1500, label = "Original\nlimit", col = "red", size = 3)

data <- data.frame(rasterToPoints(tlims_lower_all[[35]]))
colnames(data) <- c("x","y","acc_tlim")

ggplot(data, aes(x = acc_tlim)) + geom_histogram() + geom_vline(xintercept = 1, col = "blue") +
  labs(x = "Acclimatized lower thermal limit (C)", y = "Frequency") + 
  annotate(geom = "text", x = 1 + 1.5, y = 1500, label = " Original\nlimit", col = "blue", size = 3)


data_upper <- data.frame(rasterToPoints(tlims_upper_all))


tlims_upper_all %>%
  select(x, y, Eremias_argus) %>%
  mutate(Eremias_argus = as.numeric(Eremias_argus - 44.9), x = as.numeric(x), y = as.numeric(y)) %>%
  filter(!is.na(Eremias_argus)) %>%
  ggplot(., aes(x = x, y = y)) +
   coord_fixed(ratio = 1) +
  geom_raster(aes(fill = Eremias_argus)) +
  labs( y = "Latitude",
       x = "Longitude", fill = "Acclimatized limit\n-original upper limit")

tlims_lower_all %>%
  select(x, y, Eremias_argus) %>%
  mutate(Eremias_argus = as.numeric(1 - Eremias_argus), x = as.numeric(x), y = as.numeric(y)) %>%
  filter(!is.na(Eremias_argus)) %>%
  ggplot(., aes(x = x, y = y)) +
  coord_fixed(ratio = 1) +
  geom_raster(aes(fill = Eremias_argus)) +
  scale_fill_continuous(type = "gradient") +
  labs( y = "Latitude",
        x = "Longitude", fill = "Original lower limit\n- acclimatized limit")

tlims_upper_all %>%
  select(x, y, Eremias_argus) %>%
  ggplot(., aes(x = Eremias_argus)) + geom_histogram() + geom_vline(xintercept = 44.9, col = "red") +
  labs(x = "Acclimatized upper thermal limit (C)", y = "Frequency") + 
  annotate(geom = "text", x = 44.9 + 1.3, y = 75, label = "Original\nlimit", col = "red", size = 3)

tlims_lower_all %>%
  select(x, y, Eremias_argus) %>%
  ggplot(., aes(x = Eremias_argus)) + geom_histogram() + geom_vline(xintercept = 1, col = "blue") +
  labs(x = "Acclimatized lower thermal limit (C)", y = "Frequency") + 
  annotate(geom = "text", x = 1 + 1, y = 60, label = "Original\nlimit", col = "blue", size = 3)
