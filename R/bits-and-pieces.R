## code for various purposes that I cannot precisely describe in a succint sentence
library(tidyverse)



###### investigating range duplicates
duplicated <- IUCN$species[which(IUCN$species %in% GBIF$species)]

duplicate <- filter(realized_ranges, species == as.character(duplicated[1]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[1]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+2, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)

duplicate <- filter(realized_ranges, species == as.character(duplicated[2]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[2]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin-5, st_bbox(duplicate)$ymin+10, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)


duplicate <- filter(realized_ranges, species == as.character(duplicated[3]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[3]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin+8, st_bbox(duplicate)$ymin+9, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)




## looking at how complete traits are:
traits <- read_csv("data-raw/globtherm_traits_collated_180617.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")

traits <- traits[traits$genus_species %in% thermal_limits$genus_species,]

## fix wonky column names 
colnames(traits) <- str_replace_all(colnames(traits), pattern = " ", replacement = "_") %>%
  str_replace_all(., pattern = "\\(", replacement = "") %>%
  str_replace_all(., pattern = "\\)", replacement = "") %>%
  str_replace_all(., pattern = "\\;", replacement = "") %>%
  str_replace_all(., pattern = "\\/", replacement = "_") %>%
  str_replace_all(., pattern = "\\,", replacement = "") %>%
  str_replace_all(., pattern = "\\.", replacement = "") %>%
  str_replace_all(., pattern = "\\?", replacement = "")
colnames(traits)


## count how many are NA in each column 
na_count <- sapply(traits, function(y) sum(length(which(is.na(y)))))
na <- as.data.frame(na_count) %>%
  mutate(column = row.names(.)) %>%
  filter(!str_detect(column, "source")) %>%
  filter(!str_detect(column, "Source")) %>%
  filter(!str_detect(column, "notes"))

summary <- ggplot(na, aes(x = reorder(column, -na_count), y = na_count, fill = "orange")) + geom_col() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  coord_flip() + 
  labs(y = "Trait", x = "Number of values missing")

ggsave(summary, filename = "figures/traits_completeness-summary.png", device = "png")



## plotting expectations of model results:
df <- data.frame()
ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "|Latitudinal midpoint|", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60,90), labels = c("0°","30°", "60°", "90°"), limits = c(0,91))

ggsave(filename = "figures/expectations_latitudinal-midpoint.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Habitat type", y = "Range filling") +
scale_x_continuous(breaks = c(0,30,60,90, 120), labels = c("Marine", "Coastal", "Intertidal", "Terrestrial", "Freshwater"), limits = c(0,121))

ggsave(filename = "figures/expectations_habitat-type.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Range size (km2)", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,1000000))

ggsave(filename = "figures/expectations_range-size.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Trophic position", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60,90, 120), labels = c("Primary \n producer", "Herbivore", "Insectivore", "Omnivore", "Carnivore"), limits = c(0,121))

ggsave(filename = "figures/expectations_trophic-position.png", device = "png", height = 2, width = 3.5, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Body size (cm)", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,580))

ggsave(filename = "figures/expectations_body-size.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Dispersal distance category", y = "Range filling") +
  scale_x_continuous(breaks = c(0,30,60), labels = c("0-1 km", "1-10 km", "10+ km"), limits = c(0,61))

ggsave(filename = "figures/expectations_dispersal-distance.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Dispersal type category", y = "Range filling") +
  scale_x_continuous(breaks = c(0,40,80,120,160,200), labels = c("walking", "flying", "non-pelagic \n development \n and sessile adults", "non-pelagic \n development \n and crawling adults", "non-pelagic \n development \n and swimming adults", "pelagic \n  development"), limits = c(0,210))

ggsave(filename = "figures/expectations_dispersal-ability.png", device = "png", height = 3, width = 7, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Warm season dormancy", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_warm-season-dormancy.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Cold season dormancy", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_cold-season-dormancy.png", device = "png", height = 2, width = 3, units = "in")

ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Migratory", y = "Range filling") +
  scale_x_continuous(breaks = c(20, 55), labels = c("Y", "N"), limits = c(0,75))

ggsave(filename = "figures/expectations_migratory.png", device = "png", height = 2, width = 3, units = "in")


ggplot(df) + geom_point() + ylim(0, 2) +
  labs(x = "Acclimation response ratio", y = "Range filling") +
  scale_x_continuous(breaks = waiver(), limits = c(0,1.5))

ggsave(filename = "figures/expectations_range-size.png", device = "png", height = 2, width = 3, units = "in")





### code to deal with CRU portal data
absolute <- nc_open("data-raw/absolute.nc")

lat = ncvar_get(absolute, "lat")
lon = ncvar_get(absolute, "lon")
time = ncvar_get(absolute, "time")
temps = ncvar_get(absolute, "tem")

nc_close(absolute)

## save rasters of each month
raster <- raster(xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
                 crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
x = 1 
while (x  < 13) {
  slice <- temps[ , , x] 
  
  r <- raster(t(slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  raster <- addLayer(raster, r) 
  x = x+1
}

names(raster) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")

## collapse into data frame where each row is a raster square 
tmp <- data.frame(rasterToPoints(raster)) 


## mean high temperature from the warmest month
high_tmp <- tmp %>%
  mutate(seasonal_high = pmax(.$January, .$February, .$March, .$April, .$May, 
                              .$June, .$July, .$August, .$September, .$October, .$November, .$December)) %>%
  select(x, y, seasonal_high) 


## mean low temperature from the coldest month 
low_tmp <- tmp %>%
  mutate(seasonal_low = pmin(.$January, .$February, .$March, .$April, .$May, 
                             .$June, .$July, .$August, .$September, .$October, .$November, .$December)) %>%
  select(x, y, seasonal_low)

## create raster layers of high and low seasonal temperature
r = raster(nrow = nrow(lat), ncol = nrow(lon),
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
raster_high <- rasterize(high_tmp[, 1:2], r, high_tmp[,3], fun=mean)
names(raster_high) <- "seasonal_high"
plot(raster_high, asp = 1)

raster_low <- rasterize(low_tmp[, 1:2], r, low_tmp[,3], fun=mean)
names(raster_low) <- "seasonal_low"
plot(raster_low, asp = 1)


## separate ocean data from land data
countries <- ne_countries(returnclass = "sp") 
countries <- spTransform(countries, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
land_high_tmp <- crop(raster_high, extent(countries)) %>%
  mask(., countries)
ocean_high_tmp <- crop(raster_high, extent(countries)) %>%
  mask(., countries, inverse = TRUE)

land_low_tmp <- crop(raster_low, extent(countries)) %>%
  mask(., countries)
ocean_low_tmp <- crop(raster_low, extent(countries)) %>%
  mask(., countries, inverse = TRUE)

plot(ocean_high_tmp)
plot(ocean_low_tmp)
plot(land_high_tmp)
plot(land_low_tmp)





## original splitting by realm :

##### TERRESTRIAL ##### 
## repeat for marine species after 
## read in thermal limit data for each species that has both thermal tolerance metrics 

thermal_limits <- read_csv("./data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm")

upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Terrestrial")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Terrestrial")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

names_high <- c("seasonal_high", paste(both_upper$Genus, both_upper$Species, sep = "_"))
names_low <- c("seasonal_low", paste(both_lower$Genus, both_lower$Species, sep = "_"))

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
species = 1
while (species < nrow(both_upper) + 1) {
  land_high_tmp <- addLayer(land_high_tmp, land_high_tmp[[1]] - both_upper$thermal_limit[species]) 
  
  species = species + 1
}
names(land_high_tmp) <- names_high
plot(land_high_tmp)

species = 1
while (species < nrow(both_lower) + 1) {
  land_low_tmp <- addLayer(land_low_tmp, land_low_tmp[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(land_low_tmp) <- names_low
plot(land_low_tmp)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
land_high_tmp[land_high_tmp > 0] <- NA
land_low_tmp[land_low_tmp < 0] <- NA

plot(land_high_tmp)
plot(land_low_tmp)


## combine: 
combined <- land_high_tmp[[1]]
i = 2  
while (i < nrow(both_upper) + 1) {
  combined <- addLayer(combined, mask(land_high_tmp[[i]], land_low_tmp[[i]]))
  
  i = i + 1
}

combined <- combined[[-1]]
plot(combined)



##### MARINE ##### 
upper_limits <- thermal_limits %>%
  filter(type == "max") %>%
  filter(realm == "Marine")

lower_limits <- thermal_limits %>%
  filter(type == "min") %>%
  filter(realm == "Marine")

both_upper <- upper_limits[upper_limits$genus_species %in% lower_limits$genus_species,]
both_lower <- lower_limits[lower_limits$genus_species %in% upper_limits$genus_species,]

names_high <- c("seasonal_high", paste(both_upper$Genus, both_upper$Species, sep = "_"))
names_low <- c("seasonal_low", paste(both_lower$Genus, both_lower$Species, sep = "_"))

## create an individual raster layer of difference between thermal limit and seasonal temperature for each species 
species = 1
while (species < nrow(both_upper) + 1) {
  ocean_high_tmp <- addLayer(ocean_high_tmp, ocean_high_tmp[[1]] - both_upper$thermal_limit[species]) 
  
  species = species + 1
}
names(ocean_high_tmp) <- names_high
plot(ocean_high_tmp)

species = 1
while (species < nrow(both_lower) + 1) {
  ocean_low_tmp <- addLayer(ocean_low_tmp, ocean_low_tmp[[1]] - both_lower$thermal_limit[species]) 
  
  species = species + 1
}
names(ocean_low_tmp) <- names_low
plot(ocean_low_tmp)


## exclude raster cells outside of the thermal tolerance (where seasonal_high - Tmax < 0 and where seasonal_low - Tmin < 0)
ocean_high_tmp[ocean_high_tmp > 0] <- NA
ocean_low_tmp[ocean_low_tmp < 0] <- NA

plot(ocean_high_tmp)
plot(ocean_low_tmp)


## combine: 
combined <- ocean_high_tmp[[1]]
i = 2  
while (i < nrow(both_upper) + 1) {
  combined <- addLayer(combined, mask(ocean_high_tmp[[i]], ocean_low_tmp[[i]]))
  
  i = i + 1
}

combined <- combined[[-1]]
plot(combined)


plot(combined$Oncorhynchus_tshawytscha)






## garbage from 03_creating-potential-ranges.R:
###########
land_only <- nc_open("data-raw/HadCRUT.4.6.0.0.median.nc")

lat = ncvar_get(land_only, "latitude")
lon = ncvar_get(land_only, "longitude")
temps_land = ncvar_get(land_only, "temperature_anomaly")

nc_close(land_only)

## figure out which coordinates are on land vs the ocean:
map_df <- crossing(lat, lon) %>%
  arrange(-lat) %>%
  mutate(temp = as.vector(t(temps_land[ , , 1])))

x = 2
while(x < 2048) {
  layer <- as.vector(t(temps_land[ , , x]))
  map_df <- map_df %>%
    mutate(new_temp = layer) %>%
    mutate(temp = ifelse(!is.na(new_temp), new_temp, temp))
  x = x+1
}

## any cells that remain NA have no temperature data and therefore are ocean 
map_df <- map_df %>%
  mutate(is_land = !is.na(map_df$temp)) %>%
  select(lat, lon, is_land) 


## create new ocean raster (NA in land_only) and land raster(!NA in land_only)
land <- subset(map_df, is_land == TRUE) 
ocean <-  subset(map_df, is_land == FALSE) 




## raster to polygon code:
combinedtest <- combined[[-1]]
plot(combinedtest)

copy <- combinedtest
copy[copy < 0] = 1 

plot(copy[[11]], axes = TRUE) 

polygon <- copy[[11]] %>%
  rasterToPolygons(., dissolve = TRUE) %>% 
  st_as_sf()
smooth <- smooth(polygon, method = "ksmooth", smoothness = 2) 

st_crs(smooth) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

plot(smooth, axes = TRUE) 

upper <- st_crop(smooth, xmin = -180, xmax = 180, ymin = 10, ymax = 85)

st_crop(smooth, xmin = -180, xmax = 180, ymax = 10, ymin = -90)

plot(upper, axes = TRUE) 
plot(lower, axes = TRUE) 

upper_raster <- land_high_tmp[[12]]
upper_raster[upper_raster < 0] = 1 
upper_raster <- crop(upper_raster, y = c(-180, 180, -30, 15))

lower_raster <- copy[[11]]
lower_raster <- crop(lower_raster, y = c(-180, 180, -90, -30))


plot(upper_raster)


## write code to restrict clumps to those crossing a set of lines (will be midpoints)
## and then clumps within xx km of those clumps 

potential <- combined[[2]]
plot(potential)

## create multilinestring representing latitudinal midpoints of range polygons 
multiline <- st_multilinestring(list(rbind(c(-170,0),c(-115,0)), 
                                     rbind(c(75,15),c(89,15)), 
                                     rbind(c(10,-85),c(-10,-85)),
                                     rbind(c(-170,-10),c(110,-10))))

polyraster <- clump(potential, directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()

plot(st_geometry(polyraster))
plot(multiline, add = TRUE) 

intersects <- st_intersects(polyraster, range, sparse = FALSE)[,]
sub <- filter(polyraster, intersects == TRUE)

plot(sub)
plot(st_geometry(multiline), add = TRUE, col = "RED")


## checking out how CTmax vs CTmin affects potential range:
i = 1
while (i < nlayers(terr_low)) {
  sp <- str_replace(names(terr_high)[i], "_", " ")
  rr <- realized_ranges[which(as.character(realized_ranges$species) %in% sp),] %>%
    st_transform_proj(crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
  
  par(mfrow=c(1,2))
  plot(terr_high[[i]], main = paste(names(terr_high)[i], "\n Seasonal high", sep = ""))
  plot(st_geometry(rr), add = TRUE, col = "red")
  plot(terr_low[[i]], main = "Seasonal low")
  plot(st_geometry(rr), add = TRUE, col = "red")
  
  dev.copy(png, filename = paste("figures/seasonal-temp-filtering/", names(terr_high)[i], sep = ""), width = 1000, height =600, res = 150);
  dev.off() 
  
  i = i + 1
}


## investigating clumping:
clumped_temps <- raster_intertidal_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps, col = colorRampPalette(brewer.pal(8, "Accent"))(48))

clumped_temps <- raster_marine_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf()
plot(clumped_temps, col =(brewer.pal(4, "Accent")))

clumped_temps <- raster_terr_high %>%
  clump(., directions = 8) %>%
  rasterToPolygons(., dissolve = TRUE) %>%
  st_as_sf(bbox = c(xmin = -180, xmax = 180, ymax = 90, ymin = -90)) 
plot(clumped_temps, col =  colorRampPalette(brewer.pal(8, "Accent"))(61))