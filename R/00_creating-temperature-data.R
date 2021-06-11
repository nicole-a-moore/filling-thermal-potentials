library(ncdf4)
library(tidyverse)
library(raster)
library(evobiR)
library(sf)
#library(devtools)
#install_github("adamlilith/enmSdm")
library(enmSdm1)
select <- dplyr::select


####################################################################################
#####               TERRESTRIAL SEASONAL TEMPERATURE HIGH AND LOW             ######
####################################################################################
## max and min terrestrial temps:
filename <- paste("/Volumes/SundayLab/temperature-data/terrestrial/Complete_TMAX_Daily_LatLong1_1950.nc", sep = "")
ncfile <- nc_open(filename)

## create variables for things needed to use data
lat <- ncvar_get(ncfile, "latitude")
long <- ncvar_get(ncfile, "longitude")

## close the file
nc_close(ncfile)

## create arrays to hold mean temperatures from each 10 year dataset
mean_max <- array(dim = c(360, 180, 3650))
mean_max[,,] <- NaN
mean_min <- array(dim = c(360, 180, 3650))
mean_min[,,] <- NaN

rep = 1950
while (rep < 2020) {
  ## open max and minimum files 
  filename_max <- paste("/Volumes/SundayLab/temperature-data/terrestrial/Complete_TMAX_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_max <- nc_open(filename_max)
  filename_min <- paste("/Volumes/SundayLab/temperature-data/terrestrial/Complete_TMIN_Daily_LatLong1_", rep, ".nc", sep = "")
  ncfile_min <- nc_open(filename_min)
  
  ## create variables for data
  date <- ncvar_get(ncfile_max, "date_number")
  arr.anom_max <- ncvar_get(ncfile_max, "temperature")
  arr.clim_max <- ncvar_get(ncfile_max, "climatology")
  arr.anom_min <- ncvar_get(ncfile_min, "temperature")
  arr.clim_min <- ncvar_get(ncfile_min, "climatology")
  
  ## close files
  nc_close(ncfile_max)
  nc_close(ncfile_min)
  
  ## figure out which years in this time frame are leap years
  leap_years <- seq(from = rep, to = (rep+9), by = 1) %% 4 == 0
  
  ## figure out which index represents the end of each year 
  x <- c(1)
  i = 1
  while(i < length(leap_years)) {
    if(leap_years[i] == FALSE) {
      x <- append(x, x[i] + 365)
    }
    else {
      x <- append(x, x[i] + 366)
    }
    i = i+1
  }
  
  leap_years <- data.frame(leap_year = leap_years, index = x)
  
  ## create arrays to store temperatures in 
  arr.temps_max <- array(dim = c(nrow(arr.clim_max), ncol(arr.clim_max), 3650))
  arr.temps_min <- array(dim = c(nrow(arr.clim_min), ncol(arr.clim_min), 3650))
  
  ## loop through each element (unique pairs of row x column)
  row =  1
  while(row < nrow(arr.anom_max) + 1) {
    col = 1
    while(col < ncol(arr.anom_max) + 1) {
      ## retrieve climatology and anomaly in the cell 
      anom_max <- arr.anom_max[row, col, ]
      anom_min <- arr.anom_min[row, col, ]
      
      year <- 1
      this_index <- leap_years
      ## for leap years, remove element 60 within that year in anomaly data 
      while(year < nrow(leap_years) + 1) {
        if (this_index$leap_year[year] == TRUE) {
          index <- this_index$index[year]
          anom_max <- anom_max[-index+60]
          anom_min <- anom_min[-index+60]
          
          this_index$index <- this_index$index - 1
        }
        year = year + 1
      }
      ## extend climatology and add anomaly to climatology to get temperature 
      clim_max <- rep(arr.clim_max[row, col, ], times = 10)
      temps_max <- clim_max + anom_max 
      clim_min <- rep(arr.clim_min[row, col, ], times = 10)
      temps_min <- clim_min + anom_min 
      
      ## store temperatures for that cell 
      arr.temps_max[row, col, ] <- temps_max
      arr.temps_min[row, col, ] <- temps_min
      
      ## calculate average temp in each cell for each day over the period and store 
      x = 1
      iteration = (rep-1950)/10
      while (x < 366) {
        mean_max[row, col, x+365*iteration] <- mean(arr.temps_max[row, col, seq(x, 3650, 365)], 
                                                    na.rm = TRUE)
        mean_min[row, col, x+365*iteration] <- mean(arr.temps_min[row, col, seq(x, 3650, 365)], 
                                                    na.rm = TRUE)
        x = x +1
      }
      col = col + 1
    }
    row = row + 1
  }
  print(paste("On rep number:", rep))
  rep = rep + 10
}

mean_min <- readRDS("/Volumes/SundayLab/temperature-data/terrestrial/mean_min.rds")
mean_max <- readRDS("/Volumes/SundayLab/temperature-data/terrestrial/mean_max.rds")

## create temperature set matricies and lists to store data in  
## add columns for temps species with cold and hot dormancy
mean_max_new <- mean_min_new <- array(dim = c(360, 180, 365))
final_max <- final_min <- list()

element = 1
row = 1
while(row < nrow(mean_max) + 1) {
  col = 1
  while (col < ncol(mean_max) +1) {
    x = 1
    while (x < 366) {
      ## get mean min and max temp for each cell over all 70 years:
      mean_max_new[row, col, x] <- mean(mean_max[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      mean_min_new[row, col, x] <- mean(mean_min[row, col, seq(x, 3650, 365)], 
                                        na.rm = TRUE)
      x = x+1
    }
    ## get maximum mean daily temperature and minimum mean daily temperature for cell:
    max <- max(mean_max_new[row, col,], na.rm=TRUE)
    min <- min(mean_min_new[row, col,], na.rm=TRUE)
    
    ## get mean temp over the 30 days before and after to min and max temps:
    max_index2 <- first(which(mean_max_new[row, col,] == max))
    min_index2 <- first(which(mean_min_new[row, col,] == min))
    
    acc_max <- append(mean_max_new[row, col,], mean_max_new[row, col,1:183]) 
    acc_min <- append(mean_min_new[row, col,], mean_min_new[row, col,1:183]) 
    max_index1 <- ifelse(max_index2 - 30 < 0, max_index2-30+183, max_index2-30)
    min_index1 <- ifelse(min_index2 - 30 < 0, min_index2-30+183, min_index2-30)
    
    if (!is.infinite(max) & !is.infinite(min)) {
      acc_max <- mean(acc_max[max_index1:max_index2], na.rm = TRUE)
      acc_min <- mean(acc_min[min_index1:min_index2], na.rm = TRUE)
    } 
    else {
      acc_max = acc_min = NA
    }
  
    ## create df for temps and add column for month and month index
    daily_highs <- data.frame(day = c(1:365), temp = mean_max_new[row, col,], 
                              month = factor(c(rep("Jan", 31),
                                               rep("Feb", 28),
                                               rep("Mar", 31),
                                               rep("Apr", 30),
                                               rep("May", 31),
                                               rep("Jun", 30),
                                               rep("Jul", 31),
                                               rep("Aug", 31),
                                               rep("Sep", 30),
                                               rep("Oct", 31),
                                               rep("Nov", 30),
                                               rep("Dec", 31)), levels = c("Jan", "Feb", "Mar",
                                                                           "Apr", "May", "Jun",
                                                                           "Jul", "Aug", "Sep", 
                                                                           "Oct", "Nov", "Dec")), 
                              month_index = c(rep(1, 31),
                                              rep(2, 28),
                                              rep(3, 31),
                                              rep(4, 30),
                                              rep(5, 31),
                                              rep(6, 30),
                                              rep(7, 31),
                                              rep(8, 31),
                                              rep(9, 30),
                                              rep(10, 31),
                                              rep(11, 30),
                                              rep(12, 31)))
    daily_lows <- daily_highs %>%
      mutate(temp = mean_min_new[row, col,])
    
    ## get most extreme temperature high and low temp from each month, and create 'circular' data frame
    monthly_highs <- daily_highs %>%
      group_by(month) %>%
      do(mutate(., temp = max(.$temp, na.rm = TRUE))) %>%
      select(-day) %>%
      unique(.) %>%
      rbind(., .[1:6,])
    
    monthly_lows <- daily_lows %>%
      group_by(month) %>%
      do(mutate(., temp = min(.$temp, na.rm = TRUE))) %>%
      select(-day) %>%
      unique(.) %>%
      rbind(., .[1:6,])
    
    ## find 6 months with the highest mean extreme daily high temp
    sw_high <- SlidingWindow(monthly_highs$temp, window = 6, FUN = sum, step = 1)
    
    ## find 6 months with the lowest mean extreme daily low temp
    sw_low <- SlidingWindow(monthly_lows$temp, window = 6, FUN = sum, step = 1)
    
    ## get index of months to block out
    hot_start <- which(sw_high == max(sw_high))[1]
    cold_start <- which(sw_low == min(sw_low))[1]
    
    hot_end <- hot_start + 5
    cold_end <- cold_start + 5
    
    ## get month indecies to block out
    if(hot_end > 12) {
      hot_end = hot_end - 12
      hot_indecies <- c(hot_start:12, 1:hot_end)
    }
    else {
      hot_indecies <- c(hot_start:hot_end)
    }
    
    if(cold_end > 12) {
      cold_end = cold_end - 12
      cold_indecies <- c(cold_start:12, 1:cold_end)
    }
    else {
      cold_indecies <- c(cold_start:cold_end)
    }
    
    ## block out temperatures in 6 hottest consecutive months and find max
    h_d_max_high <- daily_highs %>%
      filter(!month_index %in% hot_indecies) %>%
      select(temp) %>%
      max()
    
    h_d_min_low <- daily_lows %>%
      filter(!month_index %in% hot_indecies) %>%
      select(temp) %>%
      min()
    
    ## block out temperatures in 6 coldest consecutive months and find min
    c_d_min_low <- daily_lows %>%
      filter(!month_index %in% cold_indecies) %>%
      select(temp) %>%
      min()
    
    c_d_max_high <- daily_highs %>%
      filter(!month_index %in% cold_indecies) %>%
      select(temp) %>%
      max()
    
    ## calculate mean temp before this extreme:
    max_index2 <- first(which(mean_max_new[row, col,] == h_d_max_high))
    min_index2 <- first(which(mean_min_new[row, col,] == c_d_min_low))
    
    acc_max_dormancy <- append(mean_max_new[row, col,], mean_max_new[row, col,1:26]) 
    acc_min_dormancy <- append(mean_min_new[row, col,], mean_min_new[row, col,1:26]) 
    max_index1 <- ifelse(max_index2 - 3 < 0, max_index2-3+26, max_index2-3)
    min_index1 <- ifelse(min_index2 - 3 < 0, min_index2-3+26, min_index2-3)
    
    if (!is.infinite(max) & !is.infinite(min)) {
      acc_max_dormancy <- mean(acc_max_dormancy[max_index1:max_index2], na.rm = TRUE)
      acc_min_dormancy <- mean(acc_min_dormancy[min_index1:min_index2], na.rm = TRUE)
    } 
    else {
      acc_max_dormancy = acc_min_dormancy = NA
    }
    
    
    final_max[[element]] <- c(lat[col], long[row], max, h_d_max_high, c_d_max_high, 
                              acc_max, acc_max_dormancy)
    final_min[[element]] <- c(lat[col], long[row],  min, c_d_min_low, h_d_min_low, 
                              acc_min,acc_min_dormancy)
    
    element = element + 1
    col = col + 1
  }
  row = row + 1
}

final_max <- data.frame(do.call(rbind, final_max), stringsAsFactors = FALSE) 
colnames(final_max) <- c("latitude", "longitude", "seasonal_high_temp",
                         "hot_dormancy_6mo", "cold_dormancy_6mo", "hot_acc_temp",
                         "hot_acc_temp_dormancy")
final_max <- filter(final_max, !is.infinite(seasonal_high_temp)) %>%
  select(longitude, latitude, seasonal_high_temp, hot_dormancy_6mo, cold_dormancy_6mo, 
         hot_acc_temp, hot_acc_temp_dormancy)

final_min <- data.frame(do.call(rbind, final_min), stringsAsFactors = FALSE) 
colnames(final_min) <- c("latitude", "longitude", "seasonal_low_temp",
                         "cold_dormancy_6mo", 'hot_dormancy_6mo', "cold_acc_temp", 
                         "cold_acc_temp_dormancy")
final_min <- filter(final_min, !is.infinite(seasonal_low_temp)) %>%
  select(longitude, latitude, seasonal_low_temp, cold_dormancy_6mo, hot_dormancy_6mo, cold_acc_temp, 
         cold_acc_temp_dormancy)


## save data:
write.csv(final_max, "data-processed/terrestrial_seasonal-max-temps_6mo-dormancy.csv", row.names = FALSE)
write.csv(final_min, "data-processed/terrestrial_seasonal-min-temps_6mo-dormancy.csv", row.names = FALSE)


####################################################################################
#####                 MARINE SEASONAL TEMPERATURE HIGH AND LOW                ######
####################################################################################
firstset <- "/Volumes/SundayLab/temperature-data/marine/sst.wkmean.1981-1989.nc"
secondset <- "/Volumes/SundayLab/temperature-data/marine/sst.wkmean.1990-present.nc"
landmask <- "/Volumes/SundayLab/temperature-data/marine/lsmask.nc" 

ncfile_first <- nc_open(firstset)
ncfile_second <- nc_open(secondset)
ncfile_mask <- nc_open(landmask)

## create variables for things needed to use data
lat <- ncvar_get(ncfile_first, "lat")
long <- ncvar_get(ncfile_first, "lon")
mask <- ncvar_get(ncfile_mask, "mask")

## close the files
nc_close(ncfile_first)
nc_close(ncfile_second)
nc_close(ncfile_mask)

lt_weekly_means <- array(dim = c(360, 180, 52))

firstset <-
  paste("/Volumes/SundayLab/temperature-data/marine/sst.wkmean.1981-1989.nc",
        sep = "")
secondset <-
  paste("/Volumes/SundayLab/temperature-data/marine/sst.wkmean.1990-present.nc",
        sep = "")
ncfile_first <- nc_open(firstset)
ncfile_second <- nc_open(secondset)

## week centred on Sunday
## first day is October 29, 1981
## last day is December 28, 1989
time_first <- ncvar_get(ncfile_first, "time")

## week centred on Wednesday
## first day is December 31, 1989
## last day is October 4, 2020
time_second <- ncvar_get(ncfile_second, "time")

## create variables for data
sst_first <- ncvar_get(ncfile_first, "sst")
sst_second <- ncvar_get(ncfile_second, "sst")

## close files
nc_close(ncfile_first)
nc_close(ncfile_second)

max_weekly_mean <- min_weekly_mean <- list()

## loop through each element (unique pairs of row x column)
element = row = 1
while (row < nrow(sst_first) + 1) {
  col = 1
  while (col < ncol(sst_first) + 1) {
    ## if it is a cell on land, skip it
    if (mask[row, col] == 0) {
      col = col + 1
    }
    else {
      ## retrieve temps in cell
      temps <- sst_first[row, col,]
      temps <- append(temps, sst_second[row, col,])
      
      ## calculate average temp in each cell for each week over the period and store
      x = 1
      while (x < 53) {
        lt_weekly_means[row, col, x] <- mean(temps[seq(x, 2033, 52)],
                                             na.rm = TRUE)
        
        x = x + 1
      }
      ## get maximum weekly mean temperature and minimum weekly mean temperature for cell
      max <- max(lt_weekly_means[row, col,], na.rm=TRUE)
      min <- min(lt_weekly_means[row, col,], na.rm=TRUE)
      
      ## get mean temp over the 30 days before and after to min and max temps:
      max_index2 <- first(which(lt_weekly_means[row, col,] == max))
      min_index2 <- first(which(lt_weekly_means[row, col,] == min))
      
      acc_max <- append(lt_weekly_means[row, col,], lt_weekly_means[row, col,1:26]) 
      acc_min <- append(lt_weekly_means[row, col,], lt_weekly_means[row, col,1:26]) 
      max_index1 <- ifelse(max_index2 - 3 < 0, max_index2-3+26, max_index2-3)
      min_index1 <- ifelse(min_index2 - 3 < 0, min_index2-3+26, min_index2-3)
      acc_max <- mean(acc_max[max_index1:max_index2], na.rm = TRUE)
      acc_min <- mean(acc_min[min_index1:min_index2], na.rm = TRUE)
      
      ## create df for temps, create 'circular' data frame
      weekly_temps <- data.frame(week = c(1:52), temp = lt_weekly_means[row, col,]) %>%
        rbind(., .[1:24,])
      
      ## find 6 months with the highest mean extreme daily high temp
      sw_high <- SlidingWindow(weekly_temps$temp, window = 24, FUN = sum, step = 4)
      
      ## find 6 months with the lowest mean extreme daily low temp
      sw_low <- SlidingWindow(weekly_temps$temp, window = 24, FUN = sum, step = 4)
      
      ## get index of months to block out
      hot_start <- which(sw_high == max(sw_high))[1]
      cold_start <- which(sw_low == min(sw_low))[1]
      
      hot_end <- hot_start + 23
      cold_end <- cold_start + 23
      
      ## get month indecies to block out
      if(hot_end > 52) {
        hot_end = hot_end - 52
        hot_indecies <- c(hot_start:52, 1:hot_end)
      }
      else {
        hot_indecies <- c(hot_start:hot_end)
      }
      
      if(cold_end > 52) {
        cold_end = cold_end - 52
        cold_indecies <- c(cold_start:52, 1:cold_end)
      }
      else {
        cold_indecies <- c(cold_start:cold_end)
      }
      
      ## block out temperatures in 6 hottest consecutive months and find max
      h_d_max_high <- weekly_temps %>%
        filter(!week %in% hot_indecies) %>%
        select(temp) %>%
        max()
      
      h_d_min_low <- weekly_temps %>%
        filter(!week %in% hot_indecies) %>%
        select(temp) %>%
        min()
      
      ## block out temperatures in 6 coldest consecutive months and find min
      c_d_min_low <- weekly_temps %>%
        filter(!week %in% cold_indecies) %>%
        select(temp) %>%
        min()
      
      c_d_max_high <- weekly_temps %>%
        filter(!week %in% cold_indecies) %>%
        select(temp) %>%
        max()
      
      ## calculate mean temp before this extreme:
      max_index2 <- first(which(lt_weekly_means[row, col,] == h_d_max_high))
      min_index2 <- first(which(lt_weekly_means[row, col,] == c_d_min_low))
      
      acc_max_dormancy <- append(lt_weekly_means[row, col,], lt_weekly_means[row, col,1:26]) 
      acc_min_dormancy <- append(lt_weekly_means[row, col,], lt_weekly_means[row, col,1:26]) 
      max_index1 <- ifelse(max_index2 - 3 < 0, max_index2-3+26, max_index2-3)
      min_index1 <- ifelse(min_index2 - 3 < 0, min_index2-3+26, min_index2-3)
      acc_max_dormancy <- mean(acc_max_dormancy[max_index1:max_index2], na.rm = TRUE)
      acc_min_dormancy <- mean(acc_min_dormancy[min_index1:min_index2], na.rm = TRUE)
      
      max_weekly_mean[[element]] <- c(lat[col], long[row], max, h_d_max_high, c_d_max_high, 
                                      acc_max, acc_max_dormancy)
      min_weekly_mean[[element]] <- c(lat[col], long[row],  min, c_d_min_low, h_d_min_low, 
                                      acc_min, acc_min_dormancy)
      
      element = element + 1
      col = col + 1
    }
  }
  row = row + 1
}

max_weekly_mean <- data.frame(do.call(rbind, max_weekly_mean), stringsAsFactors = FALSE) 
colnames(max_weekly_mean) <- c("latitude", "longitude", "seasonal_high_temp", 
                               "hot_dormancy_6mo", "cold_dormancy_6mo",
                               "hot_acc_temp", "hot_acc_temp_dormancy")
max_weekly_mean <- mutate(max_weekly_mean, longitude = ifelse(longitude < 180, longitude,
                                                              longitude - 360)) %>%
  select(longitude, latitude, seasonal_high_temp, hot_dormancy_6mo, cold_dormancy_6mo, 
         hot_acc_temp, hot_acc_temp_dormancy)

min_weekly_mean <- data.frame(do.call(rbind, min_weekly_mean), stringsAsFactors = FALSE) 
colnames(min_weekly_mean) <- c("latitude", "longitude", "seasonal_low_temp",
                               "cold_dormancy_6mo", "hot_dormancy_6mo", 
                               "cold_acc_temp", "cold_acc_temp_dormancy")
min_weekly_mean <- mutate(min_weekly_mean, longitude = ifelse(longitude < 180, longitude, 
                                                              longitude - 360)) %>%
  select(longitude, latitude, seasonal_low_temp, cold_dormancy_6mo, hot_dormancy_6mo, 
         cold_acc_temp, cold_acc_temp_dormancy)



## save these datasets as seasonal high and low temps
write.csv(max_weekly_mean, "data-processed/marine_seasonal-max-temps_6mo.csv", row.names = FALSE)
write.csv(min_weekly_mean, "data-processed/marine_seasonal-min-temps_6mo.csv", row.names = FALSE)

## Bio-ORACLE temperature sets:
tmin_dmax <- raster("/Volumes/SundayLab/temperature-data/marine/Present.Benthic.Max.Depth.Temperature.Lt.min.tif")
tmin_dmin <- raster("/Volumes/SundayLab/temperature-data/marine/Present.Benthic.Min.Depth.Temperature.Lt.min.tif")

ttest <- tmin_dmax <= tmin_dmin ## ok, so max depth here means deepest depth 
## want tmin at deepest (max) depth and tmax at shallowest (min) depth

tmax_dmin <- raster("/Volumes/SundayLab/temperature-data/marine/Present.Benthic.Min.Depth.Temperature.Lt.max.tif") 

## resample to 1 deg x 1 deg resolution:
tmin_dmax <- aggregate(tmin_dmax, fact = 12, fun = min, na.rm = TRUE)
tmax_dmin <- aggregate(tmax_dmin, fact = 12, fun = max, na.rm = TRUE)

## write out:
writeRaster(tmin_dmax, "data-processed/marine_min-at-max-depth.tif")
writeRaster(tmax_dmin, "data-processed/marine_max-at-min-depth.tif")


####################################################################################
#####            MARINE & TERRESTRIAL SEASONAL TEMPERATURE RASTERS            ######
####################################################################################
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## read in seasonal high and low temp data:
terr_seasonal_high <- read.csv("data-processed/terrestrial_seasonal-max-temps_6mo-dormancy.csv")
terr_seasonal_low <- read.csv("data-processed/terrestrial_seasonal-min-temps_6mo-dormancy.csv")

## rasterize:
raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3:7], fun=mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- c("seasonal_high_temp", "hot_dormancy_6mo", "cold_dormancy_6mo", 
                             "hot_acc_temp", "hot_acc_temp_dormancy")
##plot(raster_terr_high, asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3:7], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <-  c("seasonal_low_temp", "cold_dormancy_6mo", "hot_dormancy_6mo", 
                             "cold_acc_temp", "cold_acc_temp_dormancy")
##plot(raster_terr_low, asp = 1)

## restrict by zoogeographical realms:
cmec <- read_sf("data-raw/CMEC regions & realms/newRealms.shp")

# simplify the object to make it 'usable'
cmec <- cmec %>% 
  st_simplify(dTolerance = 0.01) %>% 
  group_by(Realm) %>% 
  summarise()

## rasterize using Adam's function:
rast <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
               crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
               res = 1, vals = 1)

source("R/shp2rast.R") 
for (i in 1:nrow(cmec)) {
  cur <- cmec[i,]
  shp <- as_Spatial(cur)
  
  if (i==1) {
    zoo_rast <- shp2rast(shp = shp, rast = rast) 
    zoo_rast[!is.na(zoo_rast)] = i
  }
  else {
    temp <- shp2rast(shp = shp, rast = rast)
    temp[!is.na(temp)] = i
    zoo_rast <- addLayer(zoo_rast, temp)
  }
}

mosaic <- zoo_rast[[1]]
for (i in 1:nlayers(zoo_rast)-1) {
  mosaic <- mosaic(mosaic, zoo_rast[[i+1]], fun = min)
}
plot(mosaic)

## get rid of cells in realm that there is no temperature data for:
#plot(!is.na(mosaic) & is.na(raster_terr_high[[1]])) # get rid of these 

mosaic <- mask(raster_terr_high[[1]], mosaic)

raster_terr_high_zoo <- mask(raster_terr_high, mosaic)
raster_terr_low_zoo <- mask(raster_terr_low, mosaic)

## write out:
writeRaster(raster_terr_low_zoo, "data-processed/raster_terr_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_terr_high_zoo, "data-processed/raster_terr_high.grd", 
            overwrite = TRUE, format = "raster")

## write out mask layer for use in restricting realized ranges:
raster_terr_mask <- raster_terr_high_zoo
raster_terr_mask[!is.na(raster_terr_mask)] = 1
##plot(raster_terr_mask, asp = 1)
writeRaster(raster_terr_mask, "data-processed/raster_terr_mask.grd", 
            overwrite = TRUE, format = "raster")

## write out new version of temp data, excluding cells outside of realm data:
terr_seasonal_high_zoo <- as.data.frame(rasterToPoints(raster_terr_high_zoo))
terr_seasonal_low_zoo <- as.data.frame(rasterToPoints(raster_terr_low_zoo))
names(terr_seasonal_high_zoo)[1:2] <- names(terr_seasonal_low_zoo)[1:2] <- c("longitude", "latitude")
write.csv(terr_seasonal_high_zoo, 
          "data-processed/terrestrial_seasonal-max-temps_6mo-dormancy_zoo.csv", row.names = FALSE)
write.csv(terr_seasonal_low_zoo, 
          "data-processed/terrestrial_seasonal-min-temps_6mo-dormancy_zoo.csv", row.names = FALSE)

## read in seasonal high and low temp data:
marine_seasonal_high <- read.csv("data-processed/marine_seasonal-max-temps_6mo.csv") 
marine_seasonal_low <- read.csv("data-processed/marine_seasonal-min-temps_6mo.csv")

## rasterize:
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3:7], fun=mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <-  names(raster_terr_high) 
##plot(raster_marine_high[[2]], asp = 1)

raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3:7], fun=mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <- names(raster_terr_low) 
##plot(raster_marine_low, asp = 1)

## write out mask layer for use in restricting realized ranges:
raster_marine_mask <- raster_marine_low
raster_marine_mask[!is.na(raster_marine_mask)] = 1
##plot(raster_marine_mask, asp = 1)
writeRaster(raster_marine_mask[[1]], "data-processed/raster_marine_mask.grd", 
            overwrite = TRUE, format = "raster")

## read in temps at min and max depths:
tmin_dmax <- raster("data-processed/marine_min-at-max-depth.tif")
tmax_dmin <- raster("data-processed/marine_max-at-min-depth.tif")
## mask by other marine data to stay consistent:
tmin_dmax <- mask(tmin_dmax, raster_marine_mask[[1]])
tmax_dmin <- mask(tmax_dmin, raster_marine_mask[[1]])
## add to temps:
raster_marine_high <- addLayer(raster_marine_high, tmax_dmin)
names(raster_marine_high)[6] <- "ltmax_at_min_depth"
raster_marine_low <- addLayer(raster_marine_low, tmin_dmax)
names(raster_marine_low)[6] <- "ltmin_at_max_depth"

## write out:
writeRaster(raster_marine_low, "data-processed/raster_marine_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_marine_high, "data-processed/raster_marine_high.grd", 
            overwrite = TRUE, format = "raster")


####################################################################################
#####             CREATING INTERTIDAL SEASONAL TEMPERATURE RASTERS            ######
####################################################################################

## restrict intertidal species to continental shelf around land masses:
bath <- raster("./data-processed/bathymetric-layers/raster_200mdepth_mask.grd")

## clump and turn into polygons 
cs_clumps <- bath %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(cs_clumps)

## buffer land by 1 cell:
t_poly <- raster_terr_high[[1]]
t_poly[!is.na(t_poly)] <- 1
t_poly <- t_poly %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
smooth <- smoothr::smooth(t_poly, method = "ksmooth", smoothness = 10) 
st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
## plot(st_geometry(smooth), axes = TRUE) 

buffer_land <- st_buffer(smooth, dist = 1)
plot(buffer_land)

## filter out clumps that don't intersect with the terrestrial mask:
intersects <- st_intersects(cs_clumps, buffer_land, sparse = FALSE)[,]
intersects <- filter(cs_clumps, intersects == TRUE)
plot(intersects)

## get temps within these areas:
raster_marine_low <- stack("data-processed/raster_marine_low.grd")
raster_marine_high <- stack("data-processed/raster_marine_high.grd")
raster_terr_low <- stack("data-processed/raster_terr_low.grd")
raster_terr_high <- stack("data-processed/raster_terr_high.grd")

# merge marine and terr temps, giving priority to SST in cells with both 
raster_intertidal_low <- merge(raster_marine_low[[1:5]], raster_terr_low[[1:5]]) 
names(raster_intertidal_low) <- names(raster_marine_low)[1:5]
raster_intertidal_high <- merge(raster_marine_high[[1:5]], raster_terr_high[[1:5]])
names(raster_intertidal_high) <- names(raster_marine_high)[1:5]

# mask by intertidal coninental shelf areas:
raster_intertidal_low <- mask(raster_intertidal_low, intersects)
raster_intertidal_high <- mask(raster_intertidal_high, intersects)
plot(raster_intertidal_low)
plot(raster_intertidal_high)

## write out:
writeRaster(raster_intertidal_low, "data-processed/raster_intertidal_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_intertidal_high, "data-processed/raster_intertidal_high.grd", 
            overwrite = TRUE, format = "raster")

## write out mask layer for use in restricting realized ranges:
raster_intertidal_mask <- raster_intertidal_low
raster_intertidal_mask[!is.na(raster_intertidal_mask)] = 1
##plot(raster_intertidal_mask[[1]], asp = 1)
writeRaster(raster_intertidal_mask[[1]], "data-processed/raster_intertidal_mask.grd", 
            overwrite = TRUE, format = "raster")


#############################################################################
#####               ENCORPORATING ELEVATION PERMISSIVITY               ######
#############################################################################

###########################################################
## find the minimum and mean elevation in each grid cell ##
###########################################################

## read in EarthEnv data of minimum elevations at 1km resolution
minelev <- raster("/Volumes/SundayLab/elevation-data/elevation_10KMmi_GMTEDmi.tif")
minelev[is.infinite(minelev)] <- NA
##plot(minelev, asp = 1)
## read in mean
meanelev <- raster("/Volumes/SundayLab/elevation-data/elevation_10KMmn_GMTEDmn.tif")
meanelev[is.infinite(meanelev)] <- NA
##plot(meanelev, asp = 1)

## find difference between mean and minimum elevation (mean - min) in each cell
elev_dif <- meanelev - minelev
##plot(elev_dif, asp = 1)

## aggregate to 1deg x 1deg resolution, selecting minimum value of grid cells
elev_dif <- aggregate(elev_dif, fact = 12, fun = mean, na.rm = TRUE)
##plot(elev_dif, asp = 1)

## crop raster to include only terrestrial areas using the terrestrial mask
t_mask <- raster("data-processed/raster_terr_mask.grd")
elev_dif <- extend(elev_dif, t_mask) ## extend to same extent as temperature data 
elev_dif <- mask(elev_dif, t_mask) ## masks temps outside of the terrestrial mask 
##plot(elev_dif, asp = 1)

## write raster to data frame:
elev_df <- data.frame(rasterToPoints(elev_dif))
colnames(elev_df) <- c("longitude", "latitude", "elev_dif")

write.csv(elev_df, "data-processed/xyz_elev_diff.csv", row.names = FALSE)
elev_df <- read.csv("data-processed/xyz_elev_diff.csv")

#############################################################################
## create new set of minimum temperatures at lowest elevation in each cell ##
#############################################################################
## read in terrestrial min temps 
terr_min <- read.csv("data-processed/terrestrial_seasonal-min-temps_6mo-dormancy_zoo.csv")

## merge data by geographic coordinates
terr_min <- left_join(terr_min, elev_df, by = c("latitude", "longitude"))

## apply lapse rate of +5.5deg C/km of elevation difference between mean and min elevation
terr_min <- terr_min %>%
  mutate(elev_dif = ifelse(is.na(elev_dif), 0, elev_dif)) %>%
  mutate(low_at_min_elev = ifelse(elev_dif == 0, seasonal_low_temp, 
                                  seasonal_low_temp + elev_dif*0.0055)) %>%
  mutate(low_elev_x_cold_dormancy = ifelse(elev_dif == 0, cold_dormancy_6mo, 
                                           cold_dormancy_6mo + elev_dif*0.0055)) %>%
  mutate(low_elev_x_hot_dormancy = ifelse(elev_dif == 0, hot_dormancy_6mo, 
                                          hot_dormancy_6mo + elev_dif*0.0055)) 


## plot the difference in temperatures between mean/min elev within cells:
terr_min %>%
  mutate(elev_correction_degC = low_at_min_elev - seasonal_low_temp) %>%
  ggplot(., aes(x = log(elev_correction_degC))) +
  geom_histogram()

terr_min %>%
  ggplot(., aes(x = seasonal_low_temp, y = low_at_min_elev)) +
  geom_point()

terr_min %>%
  mutate(elev_correction_degC = low_at_min_elev - seasonal_low_temp) %>%
  ggplot(., aes(x = latitude, y = elev_correction_degC)) +
  geom_point()

terr_min <- terr_min %>%
  select(-elev_dif)

## write to file
write.csv(terr_min, "data-processed/terrestrial_seasonal-min-temps_dormancy-and-elev_zoo.csv",
          row.names = FALSE)

## add layer to terrestrial low temperature raster
raster_terr_low <- rasterize(terr_min[, 1:2], r, terr_min[,3:8], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
#plot(raster_terr_low, asp = 1)

## write to raster file
writeRaster(raster_terr_low, "data-processed/raster_terr_low.grd", 
            overwrite = TRUE, format = "raster")




##############################################
#####               GARBAGE             ######
##############################################
## looking at dormancy max/min temps vs absolute max/min temps
ggplot(final_max, aes(x = seasonal_high_temp, y = hot_dormancy_6mo)) + geom_point() + 
  geom_abline(slope = 1, col = "red")

ggplot(final_min, aes(x = seasonal_low_temp, y = cold_dormancy_6mo)) + geom_point() + 
  geom_abline(slope = 1, col = "red")

length(which(final_min$seasonal_low_temp >= final_min$cold_dormancy_6mo ))
length(which(final_min$seasonal_low_temp >= final_min$cold_dormancy_6mo - 1))
length(which(final_min$seasonal_low_temp >= final_min$cold_dormancy_6mo - 5))

length(which(final_max$seasonal_high_temp == final_max$hot_dormancy_6mo))
length(which(final_max$seasonal_high_temp <= final_max$hot_dormancy_6mo + 1))
length(which(final_max$seasonal_high_temp <= final_max$hot_dormancy_6mo + 5))

ggplot(max_weekly_mean, aes(x = seasonal_high_temp, y = hot_dormancy_6mo)) + geom_point() + 
  geom_abline(slope = 1, col = "red")

ggplot(min_weekly_mean, aes(x = seasonal_low_temp, y = cold_dormancy_6mo)) + geom_point() + 
  geom_abline(slope = 1, col = "red")

length(which(min_weekly_mean$seasonal_low_temp >= min_weekly_mean$cold_dormancy_6mo ))
length(which(min_weekly_mean$seasonal_low_temp >= min_weekly_mean$cold_dormancy_6mo - 1))
length(which(min_weekly_mean$seasonal_low_temp >= min_weekly_mean$cold_dormancy_6mo - 5))

length(which(max_weekly_mean$seasonal_high_temp == max_weekly_mean$hot_dormancy_6mo))
length(which(max_weekly_mean$seasonal_high_temp <= max_weekly_mean$hot_dormancy_6mo + 1))
length(which(max_weekly_mean$seasonal_high_temp <= max_weekly_mean$hot_dormancy_6mo + 5))






## create intertidal temperature data:
## create polygon representing the edge of land:
land <- raster_terr_high 
land[is.infinite(land)] = NA 
land[land > 0 | land < 0] <- 1

polygon <- land %>%
  rasterToPolygons(., dissolve = TRUE) %>% 
  st_as_sf()
smooth <- smoothr::smooth(polygon, method = "ksmooth", smoothness = 10) 
st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
## plot(st_geometry(smooth), axes = TRUE) 

## create buffer around intertidal area into sea and onto land
buffer_sea <- st_buffer(smooth, dist = 1)
buffer_land <- st_buffer(smooth, dist = -1)
##plot(st_geometry(buffer_sea))
##plot(st_geometry(buffer_land))

## subset temperature data to include only temperatures in buffer 
intertidal_sea_high <- raster_marine_high %>%
  mask(., buffer_sea) 

intertidal_land_high <- raster_terr_high %>%
  mask(., buffer_land, inverse = TRUE) 

intertidal_sea_low <- raster_marine_low %>%
  mask(., buffer_sea) 

intertidal_land_low <- raster_terr_low %>%
  mask(., buffer_land, inverse = TRUE)

high_temp <- addLayer(intertidal_land_high, intertidal_land_high$seasonal_high_temp)
low_temp <- addLayer(intertidal_land_low, intertidal_land_low$seasonal_low_temp)

## combine land and sea temperatures, giving priority to sea temperatures 
raster_intertidal_high <- merge(intertidal_sea_high, high_temp)
raster_intertidal_low <- merge(intertidal_sea_low, low_temp)
names(raster_intertidal_high) <- names(raster_marine_high) 
names(raster_intertidal_low) <- names(raster_marine_low) 
##plot(raster_intertidal_high[[4]])
##plot(raster_intertidal_low)
## write out:
writeRaster(raster_intertidal_low, "data-processed/raster_intertidal_low.grd", 
            overwrite = TRUE, format = "raster")
writeRaster(raster_intertidal_high, "data-processed/raster_intertidal_high.grd", 
            overwrite = TRUE, format = "raster")


## write out mask layer for use in restricting realized ranges:
raster_intertidal_mask <- raster_intertidal_low
raster_intertidal_mask[!is.na(raster_intertidal_mask)] = 1
##plot(raster_intertidal_mask, asp = 1)
writeRaster(raster_intertidal_mask[[1]], "data-processed/raster_intertidal_mask.grd", 
            overwrite = TRUE, format = "raster")