######################## Install nichemapR from github
# require(devtools)
# install_github("mrke/NicheMapR")
require(NicheMapR)
# get.global.climate(folder="path") # download global climate database

######################## RUN NICHEMAPR ######################## 

## Parameters
loc <- c(x=-3.70256, y=40.4165) # set lon / lat
minshade <- 0 # minimum shade level (full sun: 0% of shade)
maxshade <- 90 # maximum shade (full shade: 90%)
Usrhyt <- 0.01 # animal's height (1 cm)

loc <- c(x=143.5, y=46.5)
micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt) # run nichemapr
micro_sun <- as.data.frame(micro$metout) # meteorological conditions in the sun 
micro_shade <- as.data.frame(micro$shadmet) # meteorological conditions in the shade
micro_soil_sun <- as.data.frame(micro$soil) # soil temperature in the sun
micro_soil_shade <- as.data.frame(micro$shadsoil) # soil temperature in the shade

plot(map)
points(x=143.5, y=46.5)

head(micro_sun)
# most important here is: DOY (day of the year), TIME (minutes), TALOC (Air temperature), RHLOC (relative humidity), VLOC (wind speed), and SOLR (solar radiation)

head(micro_soil_sun)
# soil tempearture database: here we want soil surface temperature (D0cm)

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
                        Ta,   # Air temperature (ºC)
                        Tg,   # Soil temperature (ºC)
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
Ta_sun <- micro_sun_march_aug$TALOC # air temperature (ºC)
Tg_sun <- micro_soil_sun_march_aug$D0cm  # Soil surface temperature (ºC)
v_sun <- micro_sun_march_aug$VLOC   # Wind velocity (m/s)
RH_sun <- micro_sun_march_aug$RH # relative humidity (%)

Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)

# Operative temperature in the shade
S_shade <- micro_shade_march_aug$SOLR * 0.1 # solar radiation (Wm-2)
Ta_shade <- micro_shade_march_aug$TALOC # air temperature (ºC)
Tg_shade <- micro_soil_shade_march_aug$D0cm  # Soil surface temperature (ºC)
v_shade <- micro_shade_march_aug$VLOC   # Wind velocity (m/s)
RH_shade <- micro_shade_march_aug$RH # relative humidity (%)

Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)

plot(Te_sun) # Te in full sun (march to aug)
points(Te_shade, pch=20) # Te full shade (march to aug)


######################## MAPPING OPERATIVE TEMPERATURES ######################## 

require(raster)

# I will open a WorldClim 1ºx1º map as a reference and then substitute the value of each cell by our estimated Te (don't need to do this if you have your own maps!)

r <- getData("worldclim",var="bio",res=5)
bio1 <- r$bio1
map <- aggregate(bio1, fact=1/res(bio1)[1]) # Set resolution (deg)
regions0 <- shapefile("Shapes/newRealms.shp") # this is a shapefile to crop terrestrial regions
regions <- aggregate(rbind(regions0))
map <- mask(map, regions)

plot(map)

xy <- xyFromCell(map, 1:ncell(map))
cells <- cellFromPolygon(map, regions)[[1]]
xy.values <- xy[cells,]
nrow(xy.values) # n cells

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
  micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt) # run nichemapr
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
  Ta_sun <- micro_sun_activity$TALOC # air temperature (ºC)
  Tg_sun <- micro_soil_sun_activity$D0cm  # Soil surface temperature (ºC)
  v_sun <- micro_sun_activity$VLOC   # Wind velocity (m/s)
  RH_sun <- micro_sun_activity$RH # relative humidity (%)
  
  Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
  
  # Operative temperature in the shade
  S_shade <- micro_shade_activity$SOLR * 0.1 # solar radiation (Wm-2)
  Ta_shade <- micro_shade_activity$TALOC # air temperature (ºC)
  Tg_shade <- micro_soil_shade_activity$D0cm  # Soil surface temperature (ºC)
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
for(i in 1:nrow(xy.values)){
  start <- Sys.time()
  
  # run NicheMapR
  loc <- c(x=xy.values[i,1], y=xy.values[i,2]) # set lon / lat
  tryCatch({
      micro <- micro_global(loc = loc, minshade = minshade, maxshade = maxshade, Usrhyt = Usrhyt)
      
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
  print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours"))
}

# each entry of the list "microclim_data" contains a dataframe with all the microclimatic variables we need to compute Te at each cell

# 2) Compute Te of each species (~15 min)

d <- 0.1 # body length (m)
r <- 6e5 # skin resistance to water loss (sm-1)
activity_period <- 3:8 # months of activity

Te_data <- data.frame(xy.values, cells, "maxTe_sun"=NA, "maxTe_shade"=NA, "minTe_sun"=NA, "minTe_shade"=NA)
for(i in 1:nrow(Te_data)){
  start <- Sys.time()
  
  if(!is.null(microclim_data[[i]])){
    # open microclimatic data of this cell
    microclim_data_cell <- microclim_data[[i]]
    
    # subset activity period
    month <- unique(microclim_data_cell$DOY)
    microclim_data_cell_activity <- microclim_data_cell %>% filter(DOY %in% month[activity_period])
    
    # Operative temperature in the sun
    S_sun <- microclim_data_cell_activity$S_sun # solar radiation (Wm-2)
    Ta_sun <- microclim_data_cell_activity$Ta_sun # air temperature (ºC)
    Tg_sun <- microclim_data_cell_activity$Tg_sun  # Soil surface temperature (ºC)
    v_sun <- microclim_data_cell_activity$v_sun   # Wind velocity (m/s)
    RH_sun <- microclim_data_cell_activity$RH_sun # relative humidity (%)
    
    Te_sun <- Te_function(S_sun, Ta_sun, Tg_sun, v_sun, RH_sun, d, r)
    
    # Operative temperature in the shade
    S_shade <- microclim_data_cell_activity$S_shade # solar radiation (Wm-2)
    Ta_shade <- microclim_data_cell_activity$Ta_shade # air temperature (ºC)
    Tg_shade <- microclim_data_cell_activity$Tg_shade  # Soil surface temperature (ºC)
    v_shade <- microclim_data_cell_activity$v_shade   # Wind velocity (m/s)
    RH_shade <- microclim_data_cell_activity$RH_shade # relative humidity (%)
    
    Te_shade <- Te_function(S_shade, Ta_shade, Tg_shade, v_shade, RH_shade, d, r)
    
    # Extract maximum and minimum Te
    Te_data$maxTe_sun[i] <- max(Te_sun)
    Te_data$maxTe_shade[i] <- max(Te_shade)
    Te_data$minTe_sun[i] <- min(Te_sun)
    Te_data$minTe_shade[i] <- min(Te_shade)
  }
  
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







