## wrangling bathymetric data and elevation data to use in potential range restrictions
library(tidyverse)
library(raster)
library(sf)
rasterOptions(maxmemory = 1000000000000)


#############################################################################
#####                      RESAMPLING ELEVATION                        ######
#############################################################################

##############################################################
## resample minimum and maximum elevation in each grid cell ##
##############################################################

## read in EarthEnv data of minimum elevations at 10km resolution
minelev <- raster("/Volumes/SundayLab/elevation-data/elevation_10KMmi_GMTEDmi.tif")
minelev[is.infinite(minelev)] <- NA
##plot(minelev, asp = 1)
## read in max
maxelev <- raster("/Volumes/SundayLab/elevation-data/elevation_10KMma_GMTEDma.tiff")
maxelev[is.infinite(maxelev)] <- NA
##plot(maxelev, asp = 1)

## aggregate to 1deg x 1deg resolution, selecting minimum value of grid cells
minelev <- aggregate(minelev, fact = 12, fun = min, na.rm = TRUE)
##plot(minelev, asp = 1)
maxelev <- aggregate(maxelev, fact = 12, fun = max, na.rm = TRUE)
##plot(maxelev, asp = 1)

## crop raster to include only terrestrial areas using the terrestrial mask
t_mask <- raster("./data-processed/raster_terr_mask.grd")
minelev <- extend(minelev, t_mask) ## extend to same extent as temperature data 
minelev <- mask(minelev, t_mask) ## masks temps outside of the terrestrial mask 
maxelev <- extend(maxelev, t_mask) ## extend to same extent as temperature data 
maxelev <- mask(maxelev, t_mask) ## masks temps outside of the terrestrial mask 

## ## stack and write out to file:
elevs <- stack(minelev, maxelev)
names(elevs) <- c("elev_min", "elev_max")
writeRaster(elevs, "./data-processed/elevs_minmax.grd", 
            overwrite = TRUE, format = "raster")


#############################################################################
#####                      RESAMPLING BATHYMETRY                       ######
#############################################################################
## read in bathymetry data and layer informing where land vs ocean occurs:
bath <- raster("/Volumes/SundayLab/bathymetric-data/gebco_2020_netcdf/GEBCO_2020.nc")
tid <- brick("/Volumes/SundayLab/bathymetric-data/gebco_2020_tid_netcdf/GEBCO_2020_TID.nc")
#plot(bath, asp = 1)
#plot(tid, asp = 1)

## set cells on land (0 in mask) to have elevation/depth of 0:
stack <- stack(bath, tid)

fun <- function(stack) {stack[[1]][stack[[2]] == 0] <- 0; return(stack[[1]]) }
bath <- raster::calc(stack, fun)
#plot(bath, asp = 1)

## save really long step:
writeRaster(bath, "/Volumes/SundayLab/bathymetric-data/gebco_2020_netcdf/bath_intermediate.grd", 
            overwrite = TRUE, format = "raster")
bath <- raster("/Volumes/SundayLab/bathymetric-data/gebco_2020_netcdf/bath_intermediate.grd")

tid <- aggregate(tid, fact = 240, fun = max, na.rm = TRUE) ## by assigning max value to aggregate, if aggregate area has any bit of ocean in it (value > 0) then the cell will not be considered land
tid[!tid == 0] <- 1 ## create land mask where land = 0, ocean = 1
#plot(tid, asp = 1)
writeRaster(tid, "./data-processed/bathymetric-layers/TID_agg.grd", 
            overwrite = TRUE, format = "raster")

## aggregate data to same resolution as temperature data:
depth_max <- aggregate(bath, fact = 240, fun = max, na.rm = TRUE) ## assign cell a value of the maximum depth within the aggregate cells 
depth_max[tid == 0] <- NA
writeRaster(depth_max, "./data-processed/bathymetric-layers/bath_agg_min.grd", 
            overwrite = TRUE, format = "raster")

depth_min <- aggregate(bath, fact = 240, fun = min, na.rm = TRUE) ## assign cell a value of the minimum depth within the aggregate cells 
depth_min[tid == 0] <- NA
writeRaster(depth_min, "./data-processed/bathymetric-layers/bath_agg_min.grd", 
            overwrite = TRUE, format = "raster")

## stack and write out to file:
depths <- stack(depth_min, depth_max)
names(depths) <- c("depth_min", "depth_max")
writeRaster(depths, "./data-processed/bathymetric-layers/depths_minmax.grd", 
            overwrite = TRUE, format = "raster")


################################################ 
## Creating general depth restriction layers  ##
################################################ 

## create boolean raster informing which areas are:
##      - have max that is less than -200m elev
##      - have min that is greater than 0 
##      - are defined as being on land in tid
depth <- depth_min
depth[depth_max < -200] <- NA
depth[tid == 0] <- NA
depth[depth_min > 0] <- NA
#plot(depth, asp = 1)

## looks good, save it and a mask version of it!
depth_mask <- depth
depth_mask[!is.na(depth_mask)] <- 1 ## change all areas where elev > -200m to 1
#plot(depth_mask, asp = 1)
writeRaster(depth, "./data-processed/bathymetric-layers/raster_200mdepth.grd", 
            overwrite = TRUE, format = "raster") ## with depth values 
writeRaster(depth_mask, "./data-processed/bathymetric-layers/raster_200mdepth_mask.grd", 
            overwrite = TRUE, format = "raster") ## mask 



## look at depth where occurrence records and realized ranges of most marine species are
## classify depths into equal elevation contours 

## count # occurences/area within each 





## investigating overlap between ocean and land temperature data
## right now, intertidal data gives priority to sea temperatures 
## could change so that priority given to whichever temperature is more extreme? less extreme?
## another sensitivity parameter? oh boy
terr <- raster("./data-processed/raster_terr_high.grd")
marine <- raster("./data-processed/raster_marine_high.grd")

terr_ol <- terr
terr_ol[is.na(marine)] <- NA
## plot(terr_ol, asp = 1)

marine_ol <- marine
marine_ol[is.na(terr)] <- NA
## plot(marine_ol, asp = 1)

## see how many cells
length(marine_ol[!is.na(marine_ol)]) ## 3737 cells 
length(marine[!is.na(marine)]) ## 44219 cells
length(terr[!is.na(terr)]) ## 24311 cells

## see how the air temperatures compare to sst for these cells 
ol <- terr_ol - marine_ol
## plot(ol, asp = 1)