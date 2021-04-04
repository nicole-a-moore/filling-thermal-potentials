## collating IUCN range maps
library(tidyverse)
library(rgdal)
library(sf)
library(rnaturalearth)
library(latticeExtra)

thermal_limits <- read_csv("data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm")

## do a bit of cleaning
thermal_species <- unique(paste(thermal_limits$Genus, thermal_limits$Species, sep = " "))
thermal_species[which(str_detect(thermal_species, "<ca>") == TRUE)] <- c("Apeltes quadracus", 
                                                                         "Aulacomya atra",
                                                                         "Cyprinella spiloptera",
                                                                         "Trachemys scripta")

thermal_limits$Genus[which(str_detect(thermal_limits$Genus, "<ca>") == TRUE)] <- "Trachemys"
thermal_limits$Species[which(str_detect(thermal_limits$Species, "<ca>") == TRUE)] <- c("quadracus", "atra",
                                                                                       "spiloptera", "spiloptera")
thermal_limits$genus_species <- paste(thermal_limits$Genus, thermal_limits$Species, sep = "_") 



####################################################
##              IUCN RANGES                       ##
####################################################

## figure out which fish polygon groups to download
families <- c("Pomacanthidae", "Blenniidae"	, "Albulidae", "Elopidae", "Megalopidae", "Chaetodontidae", "Epinephelidae", "Tetraodontidae", "Sparidae", "Centracanthidae", "Acanthuridae", "Syngnathidae", "Aulostomidae", "Centriscidae", "Fistulariidae", "Solenostomidae", "Istiophoridae", "Scombridae" , "Xiphiidae", "Labridae", "Scaridae")
classes <- c("Chondrichthyes", "Myxini")
orders <- c("Clupeiformes")

families_overlap <- families[which(families %in% thermal_limits$Family)]
classes_overlap <- classes[which(classes %in% thermal_limits$Class)]
orders_overlap <- orders[which(orders %in% thermal_limits$Order)]

## download sets that include: "Blenniidae", "Tetraodontidae", "Sparidae", "Labridae","Chondrichthyes", "Clupeiformes"

## sort through downloaded range maps, filter out species we do not have in thermal tolerance data
amphibs_all <- st_read("/Volumes/SundayLab/IUCN/AMPHIBIANS/AMPHIBIANS.shp")
amphibs <- thermal_species[which(thermal_species %in% amphibs_all$binomial)]

amphibs_overlap <- amphibs_all %>%
  filter(binomial %in% amphibs)
rm(amphibs_all,amphibs)

blennies_all <- st_read("/Volumes/SundayLab/IUCN/BLENNIES/BLENNIES.shp")
blennies <- thermal_species[which(thermal_species %in% blennies_all$binomial)]

blennies_overlap <- blennies_all %>%
  filter(binomial %in% blennies)
rm(blennies_all,blennies)

clups_all <- st_read("/Volumes/SundayLab/IUCN/CLUPEIFORMES/CLUPEIFORMES.shp")
clups <- thermal_species[which(thermal_species %in% clups_all$binomial)]

clups_overlap <- clups_all %>%
  filter(binomial %in% clups)
rm(clups_all,clups)

puff_all <- st_read("/Volumes/SundayLab/IUCN/PUFFERFISH/PUFFERFISH.shp")
puff <- thermal_species[which(thermal_species %in% puff_all$binomial)]

puff_overlap <- puff_all %>%
  filter(binomial %in% puff)
rm(puff_all,puff)

reptiles_all <- st_read("/Volumes/SundayLab/IUCN/REPTILES/REPTILES.shp")
reptiles <- thermal_species[which(thermal_species %in% reptiles_all$binomial)]

reptiles_overlap <- reptiles_all %>%
  filter(binomial %in% reptiles)
rm(reptiles_all,reptiles)

seabream_all <- st_read("/Volumes/SundayLab/IUCN/SEABREAMS_PORGIES_PICARELS/SEABREAMS_PORGIES_PICARELS.shp")
seabream <- thermal_species[which(thermal_species %in% seabream_all$binomial)]

seabream_overlap <- seabream_all %>%
  filter(binomial %in% seabream)
rm(seabream_all,seabream)

sharks_all <- st_read("/Volumes/SundayLab/IUCN/SHARKS_RAYS_CHIMAERAS/SHARKS_RAYS_CHIMAERAS.shp")
sharks <- thermal_species[which(thermal_species %in% sharks_all$binomial)]

sharks_overlap <- sharks_all %>%
  filter(binomial %in% sharks)
rm(sharks_all,sharks)

wrasses_all <- st_read("/Volumes/SundayLab/IUCN/WRASSES_PARROTFISHES/WRASSES_PARROTFISHES.shp")
wrasses <- thermal_species[which(thermal_species %in% wrasses_all$binomial)]

wrasses_overlap <- wrasses_all %>%
  filter(binomial %in% wrasses)
rm(wrasses_all,wrasses)


## combine and filter out non-resident ranges:
combined <- rbind(amphibs_overlap, blennies_overlap, clups_overlap, puff_overlap, 
                  reptiles_overlap, seabream_overlap, sharks_overlap, wrasses_overlap) 

## collect same ID number (species) into one MULTIPOLYGON:
combined <- aggregate(combined, list(combined$id_no), function(x) x[1])

## write out to file:
st_write(combined, "/Volumes/SundayLab/IUCN/FILTERED/IUCN-ectotherms.shp", driver = "ESRI Shapefile", 
         append = FALSE)

# ## plot one:
# one_amphib <- filter(combined, binomial == "Eurycea bislineata")
# plot(st_geometry(one_amphib))
# 
# crs(one_amphib)
# 
# ## transform it to the correct projection:
# one_amphib <- st_transform(one_amphib, crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# plot(st_geometry(one_amphib))
# crs(one_amphib)
# 
# countries <- ne_countries(returnclass = "sf") 
# countries <- st_transform(countries, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
# 
# plot(st_geometry(countries))
# plot(st_geometry(one_amphib), add = TRUE, col = "red")


##########################################################
##                      GARD RANGES                     ##
##########################################################
## checking the Meiri data for more reptile distributions 
## data from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.83s7k

## read in ranges:
GARD <- st_read("data-raw/GARD1.1_dissolved_ranges/modeled_reptiles.shp")

## check for species we have thermal limits for:
meiri_ol <- GARD$Binomial[which(GARD$Binomial %in% thermal_species)] ## 287

## taxize and check again:
library(taxize)

## meiri species:
taxa <- GARD$Binomial
taxa1 <- data.frame(binomial = taxa[1:5000])
tsn_search1 <- get_tsn(as.character(taxa1$binomial), accepted = FALSE) ## find tsn for each unique taxa
#saveRDS(tsn_search1, "data-processed/tsn_search1_meiri.rds")
tsn_search1 <- readRDS("data-processed/tsn_search1_meiri.rds")

taxa2 <- data.frame(binomial = taxa[5001:7828])
tsn_search2 <- get_tsn(as.character(taxa2$binomial), accepted = FALSE) 
#saveRDS(tsn_search2, "data-processed/tsn_search2_meiri.rds")
tsn_search2 <- readRDS("data-processed/tsn_search2_meiri.rds")

## check out weirdo Python genus
## throws error from taxize

taxa3 <- data.frame(binomial = taxa[7839:10064])
tsn_search3 <- get_tsn(as.character(taxa3$binomial), accepted = FALSE) 
#saveRDS(tsn_search3, "data-processed/tsn_search3_meiri.rds")
tsn_search3 <- readRDS("data-processed/tsn_search3_meiri.rds")

## combine:
tsns1 <- data.frame(tsn_search1)
tsns2 <- data.frame(tsn_search2)
tsns3 <- data.frame(tsn_search3)
empty <- tsns3[1:10,] 
empty[1:10,] <- NA

tsns <- rbind(tsns1, tsns2, empty, tsns3)
tsns$binomial <- taxa
#saveRDS(tsns, "data-processed/tsn_search_meiri.rds")

found <- tsns %>%
  subset(match == "found")  ## found 7400/10064 spp

report <- lapply(found$ids, itis_acceptname)
report_df <- data.frame(matrix(unlist(report), nrow=7400, byrow=T), stringsAsFactors=FALSE)
#saveRDS(report_df, "data-processed/report_df.rds")
report_df <- readRDS("data-processed/report_df.rds")

found <- found %>%
  mutate(genus_species_corrected = report_df$X2)

merged <- left_join(tsns, found) 

merged <- merged %>%
  mutate(og_name = binomial) %>%
  mutate(binomial = ifelse(!is.na(merged$genus_species_corrected), as.character(.$genus_species_corrected),
                           as.character(.$binomial)))

## thermal limit species:
taxa_tt <- data.frame(binomial = as.character(thermal_species), stringsAsFactors = FALSE)
tsn_search_tt <- get_tsn(as.character(taxa_tt$binomial), accepted = FALSE)
#saveRDS(tsn_search_tt, "data-processed/tsn_search_tt.rds")
tsn_search_tt <- readRDS("data-processed/tsn_search_tt.rds")

tsns_tt <- data.frame(tsn_search_tt)
tsns_tt$binomial <- taxa_tt$binomial

found <- tsns_tt %>%
  subset(match == "found")  ## found  spp

report <- lapply(found$ids, itis_acceptname)
report_df_tt <- data.frame(matrix(unlist(report), nrow=681, byrow=T), stringsAsFactors=FALSE)
#saveRDS(report_df_tt, "data-processed/report_df_tt.rds")
report_df_tt <- readRDS("data-processed/report_df_tt.rds")

found <- found %>%
  mutate(genus_species_corrected = report_df_tt$X2)

merged_tt <- left_join(tsns_tt, found) 

merged_tt <- merged_tt %>%
  mutate(rangetherm_names = binomial) %>%
  mutate(binomial = ifelse(!is.na(merged_tt$genus_species_corrected), as.character(.$genus_species_corrected),
                           as.character(.$binomial)))

## check again:
meiri_ol_tt <- GARD$Binomial[which(merged$binomial %in% merged_tt$binomial)] ## 304
length(which(meiri_ol %in% meiri_ol_tt)) ## okay - all previous species found did not get an updated name 

## Liolaemus kriegi and Liolaemus ceii are same species - have two ranges in GARD
## same with two other species
## for now, leave it

## yay! 304 extra species 
## for now, go with names in rangetherm for simplicity sake
inner <- inner_join(merged, merged_tt, by = "binomial")
GARD_ol <- GARD[which(merged$binomial %in% inner$binomial),]
GARD_ol$Binomial <- inner$rangetherm_names

## write: 
st_write(GARD_ol, "/Volumes/SundayLab/GARD/GARD_ol.shp", driver = "ESRI Shapefile", 
         append = FALSE)


##########################################################
##              ERROR CHECKING GBIF RANGES              ##
##########################################################
GBIF <- st_read("/Volumes/SundayLab/polygons/Filtered occurences ectotherm animals_020817.shp")

## fix species name mistakes:
GBIF$species <- as.character(GBIF$species)
GBIF$species[which(GBIF$species == "Clubiona triviali")] <- "Clubiona trivialis"

GBIF <- GBIF[GBIF$species %in% thermal_species, ] ## get rid of species not in thermal ectotherm data

## save new version:
st_write(GBIF, "/Volumes/SundayLab/polygons/gbif_error-checked.shp", append = FALSE)



##########################################################
##          combining GBIF, GARD and IUCN ranges        ##
##########################################################
IUCN <- st_read("/Volumes/SundayLab/IUCN/FILTERED/IUCN-ectotherms.shp") %>%
  select(binomial, geometry) %>%
  rename(species = binomial) 
GBIF <- st_read("/Volumes/SundayLab/polygons/gbif_error-checked.shp")
GARD <- st_read("/Volumes/SundayLab/GARD/GARD_ol.shp") %>%
  select(Binomial, geometry) %>%
  rename(species = Binomial) 

crs(IUCN)
crs(GBIF)
crs(GARD)

## add source column to identify which ranges are from GBIF vs which are from IUCN:
GBIF$source <- "GBIF"
IUCN$source <- "IUCN"
GARD$source <- "GARD"

## combine:
realized_ranges <- rbind(IUCN, GBIF, GARD) ## okay, we have 830 ectotherm ranges (but really 827)



## taking inventory of unique and overlapping species:
######################################################
length(unique(realized_ranges$species)) ## 478 unique species 
length(which(duplicated(realized_ranges$species))) ## 352 are duplicated
## 140 have only one range source 
dup <- realized_ranges$species[which(duplicated(realized_ranges$species))]

df <- data.frame(species = realized_ranges$species, source = realized_ranges$source) %>%
  count(species) 

tally_of_ol <- df %>%
  count(n)

length(unique(IUCN$species)) ## 310
length(unique(GBIF$species)) ## 216
length(unique(GARD$species)) ## 301 - 3 not unique bc taxonomy (GARD$species[which(duplicated(GARD$species))] )


## write new thermal limits database with only species that we have ranges for 
thermal_limits_new <- thermal_limits[paste(thermal_limits$Genus, thermal_limits$Species, sep = " ") %in% realized_ranges$species,]
length(unique(thermal_limits_new$genus_species)) ## 478

write.csv(thermal_limits_new, "data-processed/thermal-limits_ectotherms-with-ranges_GARD.csv", row.names = FALSE)




##########################################################
##          adding information to realized ranges       ##
##########################################################
## add extra information to realized ranges about realm, which hemisphere species is in, whether realized range crosses the equator

## split into equator-crossers, northern hemisphere and southern hemisphere ranges:
equator <- st_linestring(rbind(c(-180, 0), c(180, 0)))
n_hemi <- st_polygon(list(matrix(c(-180,0,-180,90,180,90,180,0,-180,0),ncol=2, byrow=TRUE)))
s_hemi <- st_polygon(list(matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),ncol=2, byrow=TRUE)))

equator_check <- st_intersects(realized_ranges, equator, sparse = FALSE)[,]
does_not <- filter(realized_ranges, equator_check == FALSE)

in_north <- st_intersects(does_not, n_hemi, sparse = FALSE)[,]
in_south <- st_intersects(does_not, s_hemi, sparse = FALSE)[,]

in_both <- filter(does_not, in_north == TRUE & in_south == TRUE)
## only one in the north and south, mostly in the north so consider to be in_north for now

## create sf object for ranges with information about hemisphere:
northern_ranges <- filter(does_not, in_north == TRUE) %>%
  mutate(hemisphere = "N")

southern_ranges <- filter(does_not, in_north == FALSE & in_south == TRUE) %>%
  mutate(hemisphere = "S")

crosses_equator <- filter(realized_ranges, equator_check == TRUE) %>%
  mutate(hemisphere = "EQUATOR")

realms <- thermal_limits %>%
  select(genus_species, realm) %>%
  filter(!duplicated(.))

realized_ranges <- rbind(northern_ranges, southern_ranges, crosses_equator) %>%
  left_join(realized_ranges, realms, by = c("species" = "genus_species"))

st_write(realized_ranges, "data-processed/realized-ranges_unsplit.shp", append = FALSE)





##########################################################
##                 investigating overlap                ##
##########################################################

## look at species with ranges in IUCN and made by Greta
## all 85 are squamata 
thermal_limits_new <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")
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

duplicate <- filter(realized_ranges, species == as.character(duplicated[4]))
plot(st_geometry(duplicate)[2], col = "blue", main = as.character(duplicated[4]))
plot(st_geometry(duplicate)[1], add = TRUE, col = "red") 
plot(st_geometry(countries), add = TRUE)
legend(st_bbox(duplicate)$xmin+28, st_bbox(duplicate)$ymin+22, legend=c("GBIF", "IUCN"),
       col=c("blue", "red"),  lty=1:2, cex=0.8)


## see how different IUCN ranges are from Greta's polygons:
i = 1 
props <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(props) <- c("species", "IUCN_filled", "overlap", "GBIF_filled")

while(i < length(duplicated) + 1) {
  duplicate <- filter(realized_ranges, species == as.character(duplicated[i]))
  IUCNrange <- as_Spatial(st_geometry(duplicate)[1])
  GBIFrange <- as_Spatial(st_geometry(duplicate)[2])
  
  overlap <- intersect(GBIFrange, IUCNrange)
  overlap <- area(overlap) / 10000000
  
  leftover_IUCN <- (area(IUCNrange)/10000000) - overlap 
  leftover_GBIF <-  (area(GBIFrange)/10000000) - overlap 
  row <- data.frame(species = as.character(duplicated[i]), 
                    IUCN_filled = leftover_IUCN, overlap = overlap, GBIF_filled = leftover_GBIF)
  props <- rbind(props, row)
  
  i = i + 1
}

props <- gather(props, "leftover_type", "area", -species)
props$leftover_type <- factor(props$leftover_type, levels = c("GBIF_filled", "overlap", "IUCN_filled"))

ggplot(props, aes(fill = leftover_type, y = area, x=reorder(species, -area))) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  labs(y = "Area (km^2)", x = "Species", fill = "") +
  scale_fill_manual(values =  c("darkgoldenrod1", "azure4", "darkorange3"), 
                    labels=c("GBIF range not in IUCN range",
                             "Overlap between ranges",
                             "IUCN range not in GBIF range")) +
  theme(axis.text.y = element_text(size = 6))

ggsave(device = "png", filename = "figures/investigating-iucn-gbif-overlap/IUCN-GBIF-range-overlap.png", height = 6, width = 10)

## make plot to explain what overlaps mean:
plot(IUCNrange, col = "darkorange3")
plot(GBIFrange, add = TRUE, col = "darkgoldenrod1")
plot(overlap, add = TRUE, col = "azure4")
plot(st_geometry(countries), add = TRUE)

dev.copy(png, filename = "figures/investigating-iucn-gbif-overlap/IUCN-GBIF-range-overlap_example.png", width = 1000, height = 600);
dev.off()




## garbage:

## create sf that represents the equator (a line)
equator <- st_linestring(rbind(c(-180, 0), c(180, 0)))
plot(st_geometry(equator))

## test equator crossing check:
one_species <- filter(IUCN, species == "Hylarana erythraea")
equator_check <- st_intersects(one_species, equator, sparse = FALSE)[,] ## should return TRUE
another_species <- filter(IUCN, species == "Eurycea bislineata")
equator_check <- st_intersects(another_species, equator, sparse = FALSE)[,] ## should return FALSE

both_species <- rbind(one_species, another_species) 
equator_check <- st_intersects(both_species, equator, sparse = FALSE)[,] ## should return TRUE, FALSE

## yay! it works



## see which ranges intersect the equator:
equator_check <- st_intersects(realized_ranges, equator, sparse = FALSE)[,]
crosses_equator <- filter(realized_ranges, equator_check == TRUE)

## see which remaining ranges have Northern and Southern parts
does_not <- filter(realized_ranges, equator_check == FALSE)
pts <- matrix(c(-180,-90,-180,90,180,90,180,-90,-180,-90),ncol=2, byrow=TRUE)
pts_n <- matrix(c(-180,0,-180,90,180,90,180,0,-180,0),ncol=2, byrow=TRUE)
pts_s <- matrix(c(-180,-90,-180,0,180,0,180,-90,-180,-90),ncol=2, byrow=TRUE)
both <- st_polygon(list(pts))
n_hemi <- st_polygon(list(pts_n))
s_hemi <- st_polygon(list(pts_s))
plot(both)
plot(n_hemi, add = TRUE)
plot(s_hemi, add = TRUE, col = "red")

in_north <- st_intersects(does_not, n_hemi, sparse = FALSE)[,]
in_south <- st_intersects(does_not, s_hemi, sparse = FALSE)[,]

in_both <- filter(does_not, in_north == TRUE & in_south == TRUE)
## only one in the north and south 

plot(st_geometry(countries))
plot(st_geometry(equator), add = TRUE)
plot(st_geometry(in_both)[1], add = TRUE, col = "red") 


## write out pngs of the ranges that cross the equator
i = 1
while (i < length(crosses_equator$species) + 1) {
  
  graphics.off()
  plot(st_geometry(countries), main = crosses_equator$species[i])
  plot(st_geometry(crosses_equator)[i], add = TRUE, col = "red")
  plot(st_geometry(equator), add = TRUE, col = "blue")
  
  dev.copy(png, filename=paste("data-processed/equator-crossers/", crosses_equator$species[i], ".png", sep = ""), width = 1000, height = 600);
  dev.off ();
  
  i = i + 1
}



