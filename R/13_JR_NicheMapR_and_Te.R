######################## Install nichemapR from github
# require(devtools)
# install_github("mrke/NicheMapR")
require(NicheMapR)
library(evobiR)
#get.global.climate(folder="/Volumes/SundayLab/temperature-data/terrestrial")
# download global climate database

######################## RUN NICHEMAPR ######################## 

## Parameters
loc <- c(x=-3.70256, y=40.4165) # set lon / lat
minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt) # run nichemapr
micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
micro_shade <- as.data.frame(micro$shadmet) # meteorological conditions in the shade
micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the shade

head(micro_sun)
# most important here is: DOY (day of the year), TIME (minutes), TALOC (Air temperature), RHLOC (relative humidity), VLOC (wind speed), and SOLR (solar radiation)

head(micro_soil_sun)
# soil DOMAINE LAUSANNE CAMPSITEtempearture database: here we want soil surface temperature (D0cm)

## Subsetting microclimate data
# NicheMapR computes hourly estimations of microclimatic conditions in the middle day of each month: January 15th (DOY: 15), February 15th (DOY: 46), March 15th (DOY: 74)...
# if we want to model Te during, e.g., the period March - Aug:
require(dplyr)
month <- unique(micro_sun$DOY)
micro_sun_march_aug <- micro_sun %>% filter(DOY %in% month[3:8])
micro_shade_march_aug <- micro_shade %>% filter(DOY %in% month[3:8])
micro_soil_sun_march_aug <- micro_soil_sun %>% filter(DOY %in% month[3:8])
micro_soil_shade_march_aug <- micro_soil_shade %>% filter(DOY %in% month[3:8])

# this is estimated air temperature in the sun every hour for the middle day of march, april, may... aug (each peack corresponds to one day - total 6 days/months)
plot(micro_sun_march_aug$TALOC)


######################## MODEL OPERATIVE TEMPERATURE ######################## 

# This function uses all these microclimatic data to estimate operative (=equilibrium) temperature of an animal in the sun or in the shade

Te_function <- function(S,    # Solar radiation (Wm-2)
                        Ta,   # Air temperature (?C)
                        Tg,   # Soil temperature (?C)
                        v,    # Wind speed (m/s)
                        RH,   # Relative humidity (%)
                        d,    # Body length (m)
                        r,    # Total resistance to water loss (s m-1)
                        alpha_lw=0.965, # Skin absorbance (long wave) (Buckley 2007)
                        alpha_s=0.9,    # Skin absorbance (short wave) (Buckley 2007)
                        eps=0.965,      # Skin IR emissivity (Buckley 2007)
                        Fa=0.5,         # View factor from sky (Algar et al. 2018)
                        Fg=0.5){        # View factor from the ground (Algar et al. 2018)
  
  ## Emission and absorption of long-wave radiation
  sigma = 5.67e-8 # W m-2 K-4, Stefan-Boltzmann constant
  eps_sky = 9.2e-6 * (Ta + 273)^2 # Clear sky emissivity of long wave radiation (Buckley 2007)
  La = eps_sky * sigma * (Ta + 273)^4 # long wave radiation from the sky
  
  eps_g = 0.965 # emisivity of the soil surface (Algar et al. 2018)
  Lg = eps_g * sigma * (Tg + 273)^4 # long wave radiation from the soil surface
  
  Rlw = alpha_lw * (Fa * La + Fg * Lg) # absorbed long-wave radiation
  
  ## Absorption of short-wave solar radiation
  # This is a simplification of Lauren's Buckley's (2007) model for the geometry of a lizard. 
  # Here we  assume that the the upper half of the body receives both direct and scattered solar radiation
  Rsol = alpha_s * Fa * S 
  
  ## Evaporative cooling
  ps_a = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2  # Air vapor pressure (Pa)
  rho_saturated = 2.2 * ps_a / (Ta+273) # Trasform into density (g m-3)
  RH_prop <- RH * 1e-2
  
  EWL = (rho_saturated - RH_prop * rho_saturated) / r  # Evaporative water loss (g s-1) (Spotila and Berman 1976)
  Qewl = 2257 * EWL # Evaporative cooling (W) = latent heat of vaporization of water (J g-1) x EWL (g s-1)
  
  ## Operative temperature
  cp = 29.3 # J mol-1 K-1, specific heat of the air
  Te = Ta + (Rsol + Rlw - eps * sigma * (Ta + 273)^4 - Qewl) / (4 * sigma * (Ta + 273)^3 + cp*(1.4 + 0.135*sqrt(v/d)))
  
  return(Te)
}

# Trait values
d <- 0.1 # body length (m)
r <- 6e5 # skin resistance to water loss (approx. 300 sm-1 for amphibians, and 6e5 sm-1 for lizards and other dry-skinned ectotherms)

# Operative temperature in the sun
S_sun <- micro_sun_march_aug$SOLR # solar radiation (Wm-2)
Ta_sun <- micro_sun_march_aug$TALOC # air temperature (?C)
Tg_sun <- micro_soil_sun_march_aug$D0cm  # Soil surface temperature (?C)
v_sun <- micro_sun_march_aug$VLOC   # Wind velocity (m/s)
RH_sun <- micro_sun_march_aug$RH # relative humidity (%)

Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)

# Operative temperature in the shade
S_shade <- micro_shade_march_aug$SOLR * 0.1 # solar radiation (Wm-2)
Ta_shade <- micro_shade_march_aug$TALOC # air temperature (?C)
Tg_shade <- micro_soil_shade_march_aug$D0cm  # Soil surface temperature (?C)
v_shade <- micro_shade_march_aug$VLOC   # Wind velocity (m/s)
RH_shade <- micro_shade_march_aug$RH # relative humidity (%)

Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)

plot(Te_sun) # Te in full sun (march to aug)
points(Te_shade, pch=20) # Te full shade (march to aug)


######################## MAPPING OPERATIVE TEMPERATURES ######################## 
require(raster)

## map of terrestrial areas:
# map <- raster(xmn = -180, xmx = 180, ncol = 2160, 
#               ymn = -90, ymx = 90, nrow = 1080, vals = 1) ## use map at 10' resolution
map <- raster("data-processed/raster_terr_mask.nc")
regions0 <- shapefile("data-raw/Shapes/newRealms.shp") # this is a shapefile to crop terrestrial regions
regions <- aggregate(rbind(regions0))
map <- mask(map, regions)

plot(map)

xy <- rasterToPoints(map)
cells <- which(values(map) == 1)
xy.values <- xy[,1:2]
nrow(xy.values) # n cells - 536820

# Trait and parameter values
d <- 0.1 # body length (m)
r <- 6e5 # skin resistance to water loss (sm-1)
activity_period <- 3:8 # e.g. march to august

minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

Te_data <- data.frame(xy.values, cells, "maxTe_sun"=NA, "maxTe_shade"=NA, "minTe_sun"=NA, "minTe_shade"=NA)
for(i in 1:nrow(Te_data)){
  start <- Sys.time()
  
  # run NicheMapR
  loc <- c(x=Te_data$x[i], y=Te_data$y[i]) # set lon / lat
  micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt, timeinterval = 365) # run nichemapr
  micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
  micro_shade <- as.data.frame(micro$shadmet) # and in the shade
  micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
  micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the sun
  
  # subset activity period
  month <- unique(micro_sun$DOY)
  micro_sun_activity <- micro_sun %>% filter(DOY %in% month[activity_period])
  micro_shade_activity <- micro_shade %>% filter(DOY %in% month[activity_period])
  micro_soil_sun_activity <- micro_soil_sun %>% filter(DOY %in% month[activity_period])
  micro_soil_shade_activity <- micro_soil_shade %>% filter(DOY %in% month[activity_period])
  
  # Operative temperature in the sun
  S_sun <- micro_sun_activity$SOLR # solar radiation (Wm-2)
  Ta_sun <- micro_sun_activity$TALOC # air temperature (?C)
  Tg_sun <- micro_soil_sun_activity$D0cm  # Soil surface temperature (?C)
  v_sun <- micro_sun_activity$VLOC   # Wind velocity (m/s)
  RH_sun <- micro_sun_activity$RH # relative humidity (%)
  
  Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
  
  # Operative temperature in the shade
  S_shade <- micro_shade_activity$SOLR * 0.1 # solar radiation (Wm-2)
  Ta_shade <- micro_shade_activity$TALOC # air temperature (?C)
  Tg_shade <- micro_soil_shade_activity$D0cm  # Soil surface temperature (?C)
  v_shade <- micro_shade_activity$VLOC   # Wind velocity (m/s)
  RH_shade <- micro_shade_activity$RH # relative humidity (%)
  
  Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
  
  # Extract maximum and minimum Te
  Te_data$maxTe_sun[i] <- max(Te_sun)
  Te_data$maxTe_shade[i] <- max(Te_shade)
  Te_data$minTe_sun[i] <- min(Te_sun)
  Te_data$minTe_shade[i] <- min(Te_shade)
  
  end <- Sys.time()
  elapsed <- end - start
  expected <- elapsed[[1]] * (nrow(xy.values)-i) / 60
  print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours")) # not super precise... but will give you an idea of the remaining time
}

## Making maps of max Te sun and min Te shade
map_maxTe_sun <- map_minTe_shade <- map
for(i in 1:nrow(Te_data)){
  cell <- Te_data$cells[i]
  map_maxTe_sun[cell] <- Te_data$maxTe_sun[i]
  map_minTe_shade[cell] <- Te_data$minTe_shade[i]
}

plot(map_maxTe_sun)
plot(map_minTe_shade)


######################## RUN NicheMapR and Te FOR MANY SPECIES ######################## 

# To model many species, you just need to run NicheMapR once and then use the microclimatic data to estimate Te of each species

# 1) Run NicheMapR across the world map (~4 hours)

minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

microclim_data <- list()
for (i in 1:nrow(xy.values)) { 
  start <- Sys.time()
  
  # run NicheMapR
  loc <- c(x=xy.values[i,1], y=xy.values[i,2]) # set lon / lat
  tryCatch({
    micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt, timeinterval = 365)
    
    micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
    micro_shade <- as.data.frame(micro$shadmet) # and in the shade
    micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
    micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the sun
    
    # Make a dataset with all varables we need
    microclim_data_cell <- data.frame("DOY" = micro_sun$DOY,
                                      "TIME" = micro_sun$TIME,
                                      
                                      "S_sun" = micro_sun$SOLR, 
                                      "Ta_sun" = micro_sun$TALOC, 
                                      "Tg_sun" = micro_soil_sun$D0cm,  
                                      "v_sun" = micro_sun$VLOC,   
                                      "RH_sun" = micro_sun$RH,
                                      
                                      "S_shade" = micro_shade$SOLR * 0.1, 
                                      "Ta_shade" = micro_shade$TALOC, 
                                      "Tg_shade" = micro_soil_shade$D0cm,  
                                      "v_shade" = micro_shade$VLOC,   
                                      "RH_shade" = micro_shade$RH
    )
    
    # Store microclimatic data of each cell
    microclim_data[[i]] <- microclim_data_cell
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
  
  end <- Sys.time()
  elapsed <- end - start
  expected <- elapsed[[1]] * (nrow(xy.values)-i) / 60
  print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours")) # not super precise... but will give you an idea of the remaining time
  
  print(paste("Done cell number: ", i, sep = ""))
}

saveRDS(microclim_data, "data-processed/NicheMapR-microclimate.rds")
microclim_data <- readRDS("data-processed/NicheMapR-microclimate.rds")


# each entry of the list "microclim_data" contains a dataframe with all the microclimatic variables we need to compute Te at each cell

# 2) Create microclimate sets for species that are dormant during the hot and cold season
cold_dormancy_microclim <- list()
hot_dormancy_microclim <- list()
for(i in 1:length(microclim_data)) {
  
  # open microclimatic data of this cell
  microclim_data_cell <- microclim_data[[i]]
  
  if (is.null(microclim_data_cell)) {
    hot_dormancy_microclim[[i]] <- NULL
    cold_dormancy_microclim[[i]] <- NULL
  }
  else {
    ## find which 6 consecutive months are hottest, which are coldest:
    monthly_highs <- microclim_data_cell %>%
      group_by(DOY) %>%
      do(mutate(., Ta_sun = max(.$Ta_sun, na.rm = TRUE))) %>%
      ungroup() %>%
      dplyr::select(DOY, Ta_sun) %>%
      unique(.) %>%
      rbind(., .[1:183,])
    
    monthly_lows <- microclim_data_cell %>%
      group_by(DOY) %>%
      do(mutate(., Ta_sun = min(.$Ta_sun, na.rm = TRUE))) %>%
      ungroup() %>%
      dplyr::select(DOY, Ta_sun) %>%
      unique(.) %>%
      rbind(., .[1:183,])
    
    sw_high <- SlidingWindow(monthly_highs$Ta_sun, window = 180, FUN = sum, step = 1)
    sw_low <- SlidingWindow(monthly_lows$Ta_sun, window = 180, FUN = sum, step = 1)
    
    ## get index of days to block out
    hot_start <- which(sw_high == max(sw_high))[1]
    cold_start <- which(sw_low == min(sw_low))[1]
    
    hot_end <- hot_start + 179
    cold_end <- cold_start + 179
    
    if(hot_end > 365) {
      hot_end = hot_end - 365
      hot_indecies <- c(hot_start:365, 1:hot_end)
    }
    else {
      hot_indecies <- c(hot_start:hot_end)
    }
    
    if(cold_end > 365) {
      cold_end = cold_end - 365
      cold_indecies <- c(cold_start:365, 1:cold_end)
    }
    else {
      cold_indecies <- c(cold_start:cold_end)
    }
    
    ## block out temperatures in 6 hottest and coldest consecutive months
    days <- unique(microclim_data_cell$DOY)
    
    h_d <- microclim_data_cell %>%
      filter(!DOY %in% days[hot_indecies]) 
    
    c_d <- microclim_data_cell %>%
      filter(!DOY %in% days[cold_indecies]) 
    
    ## add temps to list;
    hot_dormancy_microclim[[i]] <- h_d
    cold_dormancy_microclim[[i]] <- c_d
  }
  
  print(paste("Done cell number: ", i, sep = ""))
}

saveRDS(hot_dormancy_microclim, "data-processed/NicheMapR-microclimate_hot-dormancy.rds")
saveRDS(cold_dormancy_microclim, "data-processed/NicheMapR-microclimate_cold-dormancy.rds")
hot_dormancy_microclim <- readRDS("data-processed/NicheMapR-microclimate_hot-dormancy.rds")
cold_dormancy_microclim <- readRDS("data-processed/NicheMapR-microclimate_cold-dormancy.rds")


# 3) Compute Te of each species (~15 min)

## read in species traits and subset to terrestrial species
traits <- read.csv("data-processed/ectotherm-traits_all-spp.csv") %>%
  filter(Realm == "Terrestrial") %>%
  filter(!is.na(.$maximum_body_size_SVL_HBL_cm_)) # 344 spp (55/399 spp have NA body length)

## convert body length to m
traits$maximum_body_size_SVL_HBL_cm_ <- traits$maximum_body_size_SVL_HBL_cm_/100
Te_allspp <- list()
## for each species
for (x in 1:nrow(traits)) {
  
  sp_traits <- traits[x,]
  
  #### set species-specific parameters ####
  # body length (m)
  d <- sp_traits$maximum_body_size_SVL_HBL_cm_ 
  
  # skin resistance to water loss 
  if (sp_traits$Class == "Amphibia") {
    # approx. 300 sm-1 for amphibians
    r <- 300
  }
  else {
    # 6e5 sm-1 for lizards and other dry-skinned ectotherms
    r = 6e5
  }
  
  # if seasonal dormancy, set microclimatic data to be dormancy data
  if(sp_traits$cold_season_dormancy_ == "Yes") {
    microclim <- cold_dormancy_microclim
  }
  else if (sp_traits$hot_season_dormancy_ == "Yes") {
    microclim <- hot_dormancy_microclim
  }
  else {
    microclim <- microclim_data
  }
  
  Te_data <- data.frame(xy.values, cells, "maxTe_sun"=NA, "maxTe_shade"=NA, "minTe_sun"=NA,
                        "minTe_shade"=NA, "hot_acc_temp" = NA, "cold_acc_temp" = NA, 
                        "doy_hot" = NA, "doy_cold" = NA)
  
  for(i in 1:nrow(Te_data)){
    
    # open microclimatic data of this cell
    microclim_data_cell <- microclim[[i]]
    
    ## if no climate data for a cell, set all values of Te to NA
    if (is.null(microclim_data_cell)) {
      Te_data$maxTe_sun[i] <- NA
      Te_data$maxTe_shade[i] <- NA
      Te_data$minTe_sun[i] <- NA
      Te_data$minTe_shade[i] <- NA
    }
    else {
      # Operative temperature in the sun
      S_sun <- microclim_data_cell$S_sun # solar radiation (Wm-2)
      Ta_sun <- microclim_data_cell$Ta_sun # air temperature (?C)
      Tg_sun <- microclim_data_cell$Tg_sun  # Soil surface temperature (?C)
      v_sun <- microclim_data_cell$v_sun   # Wind velocity (m/s)
      RH_sun <- microclim_data_cell$RH_sun # relative humidity (%)
      
      Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
      
      # Operative temperature in the shade
      S_shade <- microclim_data_cell$S_shade # solar radiation (Wm-2)
      Ta_shade <- microclim_data_cell$Ta_shade # air temperature (?C)
      Tg_shade <- microclim_data_cell$Tg_shade  # Soil surface temperature (?C)
      v_shade <- microclim_data_cell$v_shade   # Wind velocity (m/s)
      RH_shade <- microclim_data_cell$RH_shade # relative humidity (%)
      
      Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
      
      # Extract maximum and minimum Te
      Te_data$maxTe_sun[i] <- max(Te_sun)
      Te_data$maxTe_shade[i] <- max(Te_shade)
      Te_data$minTe_sun[i] <- min(Te_sun)
      Te_data$minTe_shade[i] <- min(Te_shade)
      
      ## compute the mean daily air temperature the month before species experieces max and min Te
      doy_hot_extreme <- microclim_data_cell$DOY[which(Te_shade == max(Te_shade))]
      doy_cold_extreme <- microclim_data_cell$DOY[which(Te_shade == min(Te_shade))]
      
      acclim_data = microclim_data[[i]]
      days <- unique(acclim_data$DOY)
      
      if (which(days == doy_hot_extreme) == 1) {
        doys_before_hot = days[365-30]:days[365]
      }
      else {
        doys_before_hot <- days[which(days == doy_hot_extreme) - 30]:days[which(days == doy_hot_extreme) - 1]
      }
      if (which(days == doy_cold_extreme) == 1) {
        doys_before_hot = days[365-30]:days[365]
      }
      else {
        doys_before_cold <- days[which(days == doy_cold_extreme) - 30]:days[which(days == doy_cold_extreme) - 1]
      }
      
      Te_data$hot_acc_temp[i] <- filter(acclim_data, DOY %in% doys_before_hot_st) %>%
        summarize(mean(Ta_sun)) %>%
        as.numeric(as.character(.))
      Te_data$cold_acc_temp[i] <- filter(acclim_data, DOY %in% doys_before_cold) %>%
        summarize(mean(Ta_sun))%>%
        as.numeric(as.character(.))
      Te_data$doy_hot[i] = doy_before_hot
      Te_data$doy_cold[i] = doy_before_cold
      
      
    }
    
  }
  
  Te_allspp[[x]] <- Te_data
  names(Te_allspp)[x] <- as.character(sp_traits$genus_species)
  
  print(paste("Done species number: ", x, sep = ""))
}

saveRDS(Te_allspp, "data-processed/Te_allspp.rds")
Te_allspp <- readRDS("data-processed/Te_allspp.rds")

## Making maps of Te sun, Te shade and acclimated thermal limits for each species 
for (x in 1:length(Te_allspp)) {
  Te_data <- Te_allspp[[x]]
  
  map_minTe_sun <- map_maxTe_sun <- map_minTe_shade <- map_maxTe_shade <- map
  map_acc_cold <- map_acc_hot <- map
  for(i in 1:nrow(Te_data)){
    cell <- Te_data$cells[i]
    map_minTe_sun[cell] <- Te_data$minTe_sun[i]
    map_maxTe_sun[cell] <- Te_data$maxTe_sun[i]
    map_minTe_shade[cell] <- Te_data$minTe_shade[i]
    map_maxTe_shade[cell] <- Te_data$maxTe_shade[i]
    map_acc_cold[cell] <- Te_data$cold_acc_temp[i]
    map_acc_hot[cell] <- Te_data$hot_acc_temp[i]
   }
  
  #plot(map_maxTe_sun)
  #plot(map_minTe_shade)
  
  if (x == 1) {
    Te_sun_min <- map_minTe_sun
    Te_sun_max <- map_maxTe_sun
    Te_shade_min <- map_minTe_shade
    Te_shade_max <- map_maxTe_shade
    terr_acc_hot <- map_acc_hot
    terr_acc_cold <- map_acc_cold
  }
  else {
    Te_sun_min <- addLayer(Te_sun_min, map_minTe_sun)
    Te_sun_max <- addLayer(Te_sun_max, map_maxTe_sun)
    Te_shade_min <- addLayer(Te_shade_min, map_minTe_shade)
    Te_shade_max <- addLayer(Te_shade_max, map_maxTe_shade)
    terr_acc_hot <- addLayer(terr_acc_hot, map_acc_hot)
    terr_acc_cold <- addLayer(terr_acc_cold, map_acc_cold)
  }
  print(paste("Done species number: ", x, sep = ""))
}

names(Te_sun_min) <- names(Te_sun_max) <- names(Te_shade_min) <- names(Te_shade_max) <-
  names(terr_acc_cold) <- names(terr_acc_hot) <- names(Te_allspp)

## save: 
writeRaster(Te_sun_min, "data-processed/Te_sun_min.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_sun_max, "data-processed/Te_sun_max.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_shade_min, "data-processed/Te_shade_min.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_shade_max, "data-processed/Te_shade_max.grd", overwrite = TRUE, format = "raster")
writeRaster(terr_acc_cold, "data-processed/terr_acc_cold.grd", overwrite = TRUE, format = "raster")
writeRaster(terr_acc_hot, "data-processed/terr_acc_hot.grd", overwrite = TRUE, format = "raster")


## resample to our lower 1 degree by 1 degree resolution
## assign each cell the minimum Te value of the smaller pixels falling within 

## (CHECK FACTOR)
Te_sun_min <- aggregate(Te_sun_min, fact = 6, fun = min, na.rm = TRUE)
Te_sun_max <- aggregate(Te_sun_max, fact = 6, fun = min, na.rm = TRUE)
Te_shade_min <- aggregate(Te_shade_min, fact = 6, fun = min, na.rm = TRUE)
Te_shade_max <- aggregate(Te_shade_max, fact = 6, fun = min, na.rm = TRUE)

## save a mask layer:
writeRaster(nichemapr_mask, "data-processed/raster_terr_mask_nichemapr.grd", overwrite = TRUE)






## try a version that ignores dormancy:
## convert body length to m
traits$maximum_body_size_SVL_HBL_cm_ <- traits$maximum_body_size_SVL_HBL_cm_/100
Te_allspp_nodormancy <- list()
## for each species
for (x in 1:nrow(traits)) {
  
  sp_traits <- traits[x,]
  
  #### set species-specific parameters ####
  # body length (m)
  d <- sp_traits$maximum_body_size_SVL_HBL_cm_ 
  
  # skin resistance to water loss 
  if (sp_traits$Class == "Amphibia") {
    # approx. 300 sm-1 for amphibians
    r <- 300
  }
  else {
    # 6e5 sm-1 for lizards and other dry-skinned ectotherms
    r = 6e5
  }

  microclim <- microclim_data
  
  Te_data <- data.frame(xy.values, cells, "maxTe_sun"=NA, "maxTe_shade"=NA, "minTe_sun"=NA,
                        "minTe_shade"=NA)
  
  for(i in 1:nrow(Te_data)){
    
    # open microclimatic data of this cell
    microclim_data_cell <- microclim[[i]]
    
    ## if no climate data for a cell, set all values of Te to NA
    if (is.null(microclim_data_cell)) {
      Te_data$maxTe_sun[i] <- NA
      Te_data$maxTe_shade[i] <- NA
      Te_data$minTe_sun[i] <- NA
      Te_data$minTe_shade[i] <- NA
    }
    else {
      # Operative temperature in the sun
      S_sun <- microclim_data_cell$S_sun # solar radiation (Wm-2)
      Ta_sun <- microclim_data_cell$Ta_sun # air temperature (?C)
      Tg_sun <- microclim_data_cell$Tg_sun  # Soil surface temperature (?C)
      v_sun <- microclim_data_cell$v_sun   # Wind velocity (m/s)
      RH_sun <- microclim_data_cell$RH_sun # relative humidity (%)
      
      Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
      
      # Operative temperature in the shade
      S_shade <- microclim_data_cell$S_shade # solar radiation (Wm-2)
      Ta_shade <- microclim_data_cell$Ta_shade # air temperature (?C)
      Tg_shade <- microclim_data_cell$Tg_shade  # Soil surface temperature (?C)
      v_shade <- microclim_data_cell$v_shade   # Wind velocity (m/s)
      RH_shade <- microclim_data_cell$RH_shade # relative humidity (%)
      
      Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
      
      # Extract maximum and minimum Te
      Te_data$maxTe_sun[i] <- max(Te_sun)
      Te_data$maxTe_shade[i] <- max(Te_shade)
      Te_data$minTe_sun[i] <- min(Te_sun)
      Te_data$minTe_shade[i] <- min(Te_shade)
    }
    
  }
  
  Te_allspp_nodormancy[[x]] <- Te_data
  names(Te_allspp_nodormancy)[x] <- as.character(sp_traits$genus_species)
  
  print(paste("Done species number: ", x, sep = ""))
}

saveRDS(Te_allspp_nodormancy, "data-processed/Te_allspp_nodormancy.rds")
Te_allspp_nodormancy <- readRDS("data-processed/Te_allspp_nodormancy.rds")

## Making maps of max Te sun and min Te shade for each species 
for (x in 1:length(Te_allspp_nodormancy)) {
  Te_data <- Te_allspp_nodormancy[[x]]
  
  map_maxTe_sun <- map_minTe_shade <- map
  for(i in 1:nrow(Te_data)){
    cell <- Te_data$cells[i]
    map_maxTe_sun[cell] <- Te_data$maxTe_sun[i]
    map_minTe_shade[cell] <- Te_data$minTe_shade[i]
  }
  
  #plot(map_maxTe_sun)
  #plot(map_minTe_shade)
  
  if (x == 1) {
    Te_sun_nodormancy <- map_maxTe_sun
    Te_shade_nodormancy <- map_minTe_shade
  }
  else {
    Te_sun_nodormancy <- addLayer(Te_sun_nodormancy, map_maxTe_sun)
    Te_shade_nodormancy <- addLayer(Te_shade_nodormancy, map_minTe_shade)
  }
  print(paste("Done species number: ", x, sep = ""))
}

names(Te_sun_nodormancy) <- names(Te_allspp_nodormancy)
names(Te_shade_nodormancy) <- names(Te_allspp_nodormancy)

## save: 
writeRaster(Te_sun_nodormancy, "data-processed/Te_sun_nodormancy.grd", overwrite = TRUE, format = "raster")
writeRaster(Te_shade_nodormancy, "data-processed/Te_shade_nodormancy.grd", overwrite = TRUE, format = "raster")



###### How do mean Te_sun and mean Te_shade comare to highest and lowest climatologies? ###### 
Te_sun <- stack("data-processed/Te_sun_nodormancy.grd")
Te_shade <- stack("data-processed/Te_shade_nodormancy.grd")

mean_sun <- merge(Te_sun, fun = mean)
mean_shade <- merge(Te_shade, fun = mean)
plot(mean_sun)
plot(mean_shade)

clim_low <- stack("data-processed/raster_terr_low.grd")
clim_high <- stack("data-processed/raster_terr_high.grd")
plot(clim_low[[1]])
plot(clim_high[[1]])

library(tidyverse)
library(sf)

clim_low <- mask(clim_low[[1]], mean_shade)
clim_high <-  mask(clim_high[[1]], mean_sun)
mean_shade <- mask(mean_shade, clim_low[[1]])
mean_sun <- mask(mean_sun, clim_high[[1]])


anom_low <-  mean_shade -  clim_low
anom_high <- mean_sun - clim_high 
plot(anom_low)  
plot(anom_high) 

anom_low_xy <- as.data.frame(rasterToPoints(anom_low))
anom_high_xy <- as.data.frame(rasterToPoints(anom_high))

## places where Te is colder than climatology = negative values = blue
## places where Te is warmer than climatology = positive values = red
## mostly colder, large difference at high Northern latitudes
anom_low_xy %>%
  ggplot(., aes(x = x, y = y)) +
  geom_tile(aes(fill=layer)) + theme_classic() + 
  labs(fill = "min Te shade - min T air") +
  coord_equal() +
  scale_fill_gradient2(
    low = "skyblue",
    mid = "grey94",
    high = "darkred")

## places where climatology is colder than Te = negative values = blue
## places where climatology is warmer than Te = positive values = red
## only colder
anom_high_xy %>%
  ggplot(., aes(x = x, y = y)) +
  geom_tile(aes(fill=layer)) + theme_classic() + 
  coord_equal() +
  labs(fill = "max Te sun - max T air") +
  scale_fill_gradient2(
    low = "skyblue",
    mid = "grey94",
    high = "darkred") 

## explains patterns of the changes in potential ranges:
# 1. australian ranges become much more precise - operative shade temp and cold climatology in this area are comparable, but operative sun temp much closer to upper thermal limit 
# 2. potential ranges of species with realized ranges in europe and russia shift northward - operative temp in shade is much warmer than cold climatology in these areas 
# 3. potential ranges of species with realized ranges in North America get very small


mean_sun <- as.data.frame(rasterToPoints(mean_sun))
mean_shade <- as.data.frame(rasterToPoints(mean_shade))

## places where climatology is colder than Te = negative values = blue
## places where climatology is warmer than Te = positive values = red
## mostly colder, large difference at high Northern latitudes
mean_sun %>%
  ggplot(., aes(x = x, y = y)) +
  geom_tile(aes(fill=layer)) + theme_classic() + 
  labs(fill = "Mean max Te in sun") +
  coord_equal() +
  scale_fill_gradient2(
    low = "skyblue",
    mid = "grey",
    high = "darkred")

mean_shade %>%
  ggplot(., aes(x = x, y = y)) +
  geom_tile(aes(fill=layer)) + theme_classic() + 
  labs(fill = "Min Te in shade") +
  coord_equal() +
  scale_fill_gradient2(
    low = "skyblue",
    mid = "grey",
    high = "darkred")

## places where climatology is colder than Te = negative values = blue
## places where climatology is warmer than Te = positive values = red
## only colder
anom_high_xy %>%
  ggplot(., aes(x = x, y = y)) +
  geom_tile(aes(fill=layer)) + theme_classic() + 
  coord_equal() +
  scale_fill_gradient2(
    low = "skyblue",
    mid = "white",
    high = "darkred")



## plotting distribution of sun and shade temps
library(tidyverse)
te_data <- do.call(rbind, Te_allspp)

ggplot(data = te_data, aes(x = maxTe_sun)) + geom_histogram()
ggplot(data = te_data, aes(x = minTe_sun)) + geom_histogram()
ggplot(data = te_data, aes(x = maxTe_shade)) + geom_histogram()
ggplot(data = te_data, aes(x = minTe_shade)) + geom_histogram()


## look at only cordylus spp
# 111,112,113

c_c <- Te_allspp[[111]]
c_n <- Te_allspp[[112]]
c_o <- Te_allspp[[113]]

c_c_gathered <- c_c %>%
  gather(key = 'Te_type', value = "Te", c(maxTe_sun, maxTe_shade))

c_c_gathered %>%
  ggplot(data = ., aes(x = Te, fill = Te_type)) + geom_histogram(alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("darkred", "red")) + 
  geom_vline(xintercept = 32.1) + # tsel
  geom_vline(xintercept = 41.6) # ctmax

length(which(c_c$maxTe_sun > 41.6 & c_c$maxTe_shade < 32.1))
## 2525/14918 = 16.9%

# ctmin 7.5
ggplot(data = c_c_gathered, aes(x = Te, y = minTe_shade, colour = Te_type)) + geom_point() +
  scale_colour_manual(values = c("darkred", "red"))

length(which(c_c$maxTe_sun > 41.6 & c_c$maxTe_shade < 32.1 & c_c$minTe_shade > 7.5))
## 767/14918 = 5.15%
       