## comparing GBIF occurrence data to IUCN ranges 
library(tidyverse)
library(sf)

## read in IUCN ranges
iucn <- st_read("/Volumes/SundayLab/IUCN/FILTERED/IUCN-ectotherms.shp") %>%
  select(binomial, geometry) %>%
  rename(species = binomial) 

## read in gbif ranges
gbif <- st_read("/Volumes/SundayLab/polygons/Filtered occurences ectotherm animals_020817.shp")
spp <- as.character(gbif$species)

## update species synonyms to match
spp[which(spp == "Norops auratus")] <- "Anolis auratus" 
spp[which(spp == "Norops dollfusianus")] <- "Anolis cupreus"
spp[which(spp == "Norops humilis")] <- "Anolis humilis"
spp[which(spp == "Norops lemurinus")] <- "Anolis lemurinus"
spp[which(spp == "Norops limifrons")] <- "Anolis limifrons"
spp[which(spp == "Norops pentaprion")] <- "Anolis pentaprion"
spp[which(spp == "Norops sagrei")] <- "Anolis sagrei"
spp[which(spp == "Norops tropidogaster")] <- "Anolis tropidogaster"
spp[which(spp == "Norops tropidolepis")] <- "Anolis tropidolepis"
spp[which(spp == "Clubiona triviali")] <- "Clubiona trivialis"
spp[which(spp == "Cyprinodon nevadensis")] <- "Cyprinodon nevadensis_amargosae"
spp[which(spp == "Liolaemus bellii")] <- "Liolaemus altissimus"
spp[which(spp == "Liolaemus andinus")] <- "Liolaemus poecilochromus"
spp[which(spp == "Prochilodus lineatus")] <- "Prochilodus scrofa"
spp[which(spp == "Saproscincus mustelinus")] <- "Saproscincus mustelina"
spp[which(spp == "Saproscincus tetradactylus")] <- "Saproscincus tetradactyla"
spp[which(spp == "Pagothenia borchgrevinki")] <- "Trematomus borchgrevinki"
spp[which(spp == "Gnypetoscincus queenslandiae")] <- "Tropidophorus queenslandiae"
spp[which(spp == "Lophozozymus incisus")] <- "Xantho incisus"

## read in filtered occurrence records used to create GBIF ranges 
file = 1

while (file < length(spp) + 1) {
  sp <- read.csv(paste("/Volumes/SundayLab/dropbox filtered occurrences for sWEEP/Filtered occurences ectotherm animals 16 12 31 copy/",
                       spp[file], ".csv", sep = ""))
  
  ## bind occurrence rows to new df:
  if (file == 1) {
    occ <- sp
  }
  else {
    occ <- rbind(occ, sp)
  }
  
  file = file + 1
}


## subset to species for which we have IUCN ranges:
occ_ol <- occ %>%
  filter(species %in% iucn$species)

length(unique(occ_ol$species)) ## 85 spp


## for each range:
range = 1
while (range < length(unique(occ_ol$species))+1) {
  cur <- filter(occ_ol, occ_ol$species == unique(occ_ol$species)[range])
  
  ## retrieve IUCN range:
  poly_i <- iucn %>%
    filter(species == as.character(unique(occ_ol$species)[range]))
  
  ## retrieve GBIF range:
  poly_g <- gbif %>%
    filter(species == as.character(unique(occ_ol$species)[range]))
  
  ## merge extents to plot together:
  e = merge(extent(poly_i),extent(poly_g))
  plot(e, type="n", xlab = "longitude", ylab = "latitude")
  plot(poly_g, col = "darkgoldenrod1", add=TRUE)
  plot(poly_i, col = "darkorange3", add=TRUE)
  
  ## transform occurence data into spatial data:
  pts <- data.frame(lon = cur$lon, lat = cur$lat)
  pts = matrix(ncol = 2, c(cur$lon, cur$lat))
  spat_occ <- st_multipoint(x = pts)
  spat_occ <- st_cast(st_sfc(spat_occ, crs = as.character(crs(gbif))), "POINT")
  plot(spat_occ, add = TRUE)
  
  ## find how many points fall inside the IUCN and the GBIF range 
  in_iucn <- spat_occ[st_within(spat_occ, poly_i) %>% lengths > 0,]
  plot(in_iucn, add = TRUE, col = "yellow")
  
  in_gbif <- spat_occ[st_within(spat_occ, poly_g) %>% lengths > 0,]
  
  ## save plot:
  dev.copy(png, filename = paste("./figures/investigating-iucn-gbif-overlap/", 
                                 str_replace(unique(cur$species), " ", "_"), ".png", sep = ""),
           width = 1000, height = 600)
  dev.off()
  dev.off()
  
  ## add information to data frame:
  if (range == 1) {
    error <- data.frame(species = unique(cur$species), num_occ = nrow(pts), 
                        num_occ_in_iucn = length(in_iucn), num_occ_in_gbif = length(in_gbif))
  }
  else {
    error <- rbind(error, data.frame(species = unique(cur$species), num_occ = nrow(pts), 
                                     num_occ_in_iucn = length(in_iucn), 
                                     num_occ_in_gbif = length(in_gbif)))
  }
  
  ## move to next species:
  range = range + 1
}


## save data frame
write.csv(error, "data-processed/occurence-overlap.csv", row.names= FALSE)

## calc. 



## garbage:
## read in thermal limit data
thermal <- read.csv("data-raw/globtherm_full_dataset_2019.csv")


