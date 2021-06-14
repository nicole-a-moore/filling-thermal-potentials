## for cleaning and manipulating the traits database 
library(tidyverse)
library(sf)

## read in all traits
traits_all <- read_csv("data-raw/globtherm_traits_collated_180617.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv")

## check that no species we have ranges for are missing from the trait database:
length(which(!thermal_limits$genus_species %in% traits_all$genus_species))
## 7 rows species missing! Which ones? 
new <- unique(thermal_limits$genus_species[which(!thermal_limits$genus_species %in% traits_all$genus_species)]) ## 4 species
new_infos <- thermal_limits[which(thermal_limits$genus_species %in% new),] 

types <- new_infos %>%
  select(genus_species, type) %>%
  unique() 

types <- aggregate(types$type, list(types$genus_species), paste, collapse = ", ") %>%
  dplyr::rename("genus_species" = Group.1, "limit_type" = x)

new_infos <- new_infos %>%
  select(-type) %>%
  filter(!duplicated(genus_species)) %>%
  left_join(., types, by = "genus_species") 

## subset to only ectotherm species for which we have thermal limits and range
traits_sub <- traits_all[traits_all$genus_species %in% thermal_limits$genus_species,]

## add missing species:
new_species <- traits_sub[1:length(new),1:40] 
new_species[1:length(new),] <- NA 
new_species <- new_species %>%
  mutate(genus_species = new_infos$genus_species) %>%
  mutate(Genus = new_infos$Genus, Species = new_infos$Species, Family = new_infos$Family,
         Phylum = new_infos$Phylum, Class = new_infos$Class, 
         Order = new_infos$Order, Realm = new_infos$realm) %>%
  mutate('data gatherer' = "Nikki") 

traits_sub <- rbind(traits_sub, new_species)

## make new column saying whether species has one or both thermal limits 
lim_types <- thermal_limits %>%
  select(genus_species, type) %>%
  unique() 

lim_types <- aggregate(lim_types$type, list(lim_types$genus_species), paste, collapse = ", ") %>%
  dplyr::rename("genus_species" = Group.1, "limit_type" = x)

traits_sub <- left_join(traits_sub, lim_types, by = "genus_species")


## write out and start filling in the missing ones!!
write.csv(traits_sub, "data-processed/globtherm_traits_collated_180617_ectotherms-with-limits.csv", row.names = FALSE)


## fix wonky column names 
oldnames <- colnames(traits_sub)
colnames(traits_sub) <- str_replace_all(colnames(traits_sub), pattern = " ", replacement = "_") %>%
  str_replace_all(., pattern = "\\(", replacement = "") %>%
  str_replace_all(., pattern = "\\)", replacement = "") %>%
  str_replace_all(., pattern = "\\;", replacement = "") %>%
  str_replace_all(., pattern = "\\/", replacement = "_") %>%
  str_replace_all(., pattern = "\\,", replacement = "") %>%
  str_replace_all(., pattern = "\\.", replacement = "") %>%
  str_replace_all(., pattern = "\\?", replacement = "")
colnames(traits_sub)

dormancy_sub <- traits_sub %>%
  filter(limit_type == "max, min") %>%
  filter(!is.na(cold_season_dormancy) & !is.na(hot_season_dormancy)) %>%
  mutate(cold_season_dormancy = ifelse(str_detect(cold_season_dormancy, "N") | 
                                         str_detect(cold_season_dormancy, "no"), 
                                       "No", ifelse(str_detect(cold_season_dormancy, "Y"), 
                                                    "Yes", cold_season_dormancy)
                                       
  )) %>%
  mutate(hot_season_dormancy = ifelse(str_detect(hot_season_dormancy, "N") | 
                                        str_detect(hot_season_dormancy, "no"), 
                                      "No", hot_season_dormancy))


unique(dormancy_sub$cold_season_dormancy)
unique(dormancy_sub$hot_season_dormancy)

# ## check how many are yes and no
# length(which(dormancy_sub$cold_season_dormancy == "Yes")) ## 39/141 
# length(which(dormancy_sub$hot_season_dormancy == "Yes")) ## 0/141
# 
# ## check how many of those are terrestrial
# terr_sub <- dormancy_sub %>%
#   filter(Realm == "Terrestrial") 
# 
# length(which(terr_sub$cold_season_dormancy == "Yes")) ## 39/112, all of them!

# ## plot dormancy Y/N against realized range latitudinal midpoint (should end up with 141 species)
# realized_ranges <- st_read("data-processed/realized-ranges_split.shp")
# 
# lat_mps <- realized_ranges %>%
#   filter(!duplicated(species)) %>%
#   mutate(species = str_replace_all(species, " ", "_")) %>%
#   filter(species %in% dormancy_sub$genus_species) %>%  ## have 140 since 1 species is duplicated in traits data... oops
#   as.data.frame(.) %>%
#   select(species, lttdnl_) %>%
#   rename("lat_mp" = lttdnl_)
#   

# dormancy_sub <- left_join(dormancy_sub, lat_mps, by = c("genus_species" = "species"))
# 
# dormancy_sub %>%
#   filter(Realm == "Terrestrial") %>%
#   ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy, col = Order)) + geom_point() ## very interesting.... seems like they are all squamata
# 
# unique(dormancy_sub$Order[which(dormancy_sub$cold_season_dormancy == "Yes")]) ## it is true!
# 
# ## now look at only squamata: 
# dormancy_sub %>%
#   filter(Order == "Squamata") %>% ## 41/111 squamata species are dormant in the cold seasons
#   ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy, col = Order)) + geom_point()
# 
# ## print out names 
# dormancy_sub %>%
#   filter(Order == "Squamata")  %>%
#   filter(cold_season_dormancy == "Yes") %>%
#   .$genus_species




## updated traits:
new_traits <- read.csv("./data-processed/globtherm_traits_collated_180617_ectotherms-with-limits_filled.csv") %>%
  mutate(genus_species = str_replace_all(.$genus_species, "_", " "))
colnames(new_traits) <- str_replace_all(colnames(new_traits), pattern = "\\.", "_")
colnames(new_traits) <- str_replace_all(colnames(new_traits), pattern = "\\__", "_")


dormancy_sub <- new_traits %>%
  mutate(cold_season_dormancy_ = ifelse(str_detect(cold_season_dormancy_, "N") | 
                                         str_detect(cold_season_dormancy_, "no"), 
                                       "No", ifelse(str_detect(cold_season_dormancy_, "Y"), 
                                                    "Yes", cold_season_dormancy_))) %>%
  mutate(hot_season_dormancy_ = ifelse(str_detect(hot_season_dormancy_, "N") | 
                                        str_detect(hot_season_dormancy_, "no"), 
                                      "No", ifelse(str_detect(hot_season_dormancy_, "Y"), 
                                                   "Yes", hot_season_dormancy_))) %>%
  mutate(cold_season_dormancy_ = ifelse(is.na(cold_season_dormancy_), "No", cold_season_dormancy_)) %>%
  mutate(hot_season_dormancy_ = ifelse(is.na(hot_season_dormancy_), "No", hot_season_dormancy_))

## clean trophic position:
dormancy_sub <- dormancy_sub %>%
  mutate(Trophic_position = ifelse(Trophic_position == "insectivore, omnivore", "omnivore", 
                                   ifelse(Trophic_position == "insectivore, carnivore", "carnivore", 
                                          ifelse(Trophic_position == "herbiovore", "herbivore", 
                                                 as.character(Trophic_position)))))

## clean dispersal distance category:
dormancy_sub <- dormancy_sub %>%
  mutate(dispersal_distance_category = ifelse(dispersal_distance_category == "1-10km", "1-10", 
                                              ifelse(dispersal_distance_category == "1--10", "1-10",
                                                     ifelse(dispersal_distance_category == "0-1km", "0-1",
                                                 as.character(dispersal_distance_category)))))

## fix duplicated species by getting rid of ones where realm is not the realm listed in globtherm
dup_rows <- dormancy_sub[which(dormancy_sub$genus_species %in% 
                                 dormancy_sub$genus_species
                               [which(duplicated(dormancy_sub$genus_species))]), ]
dormancy_sub[28,] <- dormancy_sub[79,]
dormancy_sub$Realm[28] <- "Intertidal"
dormancy_sub <- dormancy_sub[-c(79, 29), ]

## make sure all realms are the same as in thermal limit data:
r_lims <- select(thermal_limits, genus_species, realm) %>%
  mutate(genus_species = str_replace_all(genus_species, "_", " ")) %>%
  filter(!duplicated(.))
r_traits <- select(dormancy_sub, genus_species, Realm)

realm_merge <- left_join(r_traits, r_lims) %>%
  mutate(Realm = realm) %>%
  select(-realm)%>%
  left_join(., select(dormancy_sub, -Realm))

## write to file:
write.csv(realm_merge, "./data-processed/wrangled-traits.csv", row.names = FALSE)





# dormancy_sub <- dormancy_sub %>%
#   filter(limit_type == "max, min") %>%
#   filter(!is.na(cold_season_dormancy_) & !is.na(hot_season_dormancy_)) 
# 
# unique(dormancy_sub$cold_season_dormancy_)
# unique(dormancy_sub$hot_season_dormancy_)
# 
# ## check how many are yes and no
# length(which(dormancy_sub$cold_season_dormancy_ == "Yes")) ## 119/186 
# length(which(dormancy_sub$hot_season_dormancy_ == "Yes")) ## 5/186
# 
# lat_mps <- realized_ranges %>%
#   filter(!duplicated(species)) %>%
#   filter(species %in% dormancy_sub$genus_species) %>%  ## have 185 since 1 species is duplicated in traits data... oops
#   as.data.frame(.) %>%
#   select(species, lttdnl_) %>%
#   rename("lat_mp" = lttdnl_)
# 
# dormancy_sub <- left_join(dormancy_sub, lat_mps, by = c("genus_species" = "species"))
# 
# dormancy_sub %>%
#   filter(Realm == "Terrestrial") %>%
#   ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy_, col = Order)) + geom_point() 
# 
# unique(dormancy_sub$Order[which(dormancy_sub$cold_season_dormancy == "Yes")]) ## no longer all squamata!
# 
# ## now look at only squamata: 
# dormancy_sub %>%
#   filter(Order == "Squamata") %>%
#   ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy_, col = Order)) + geom_point()
# ## more higher latitude species have "Yes" after updating 
# 
# 
# ## look at logic for 'No' ... did we guess these based on climate?
# dormancy_sub %>%
#   filter(Realm == "Terrestrial") %>%
#   filter(cold_season_dormancy_ == "No") %>%
#   select(source_logic_for_cold_season_dormancy) %>%
#   View(.)
# ## about half have a cited source... leave for now
# 
# 
# ## look at logic for 'Yes'... did we guess these based on climate?
# dormancy_sub %>%
#   filter(Realm == "Terrestrial") %>%
#   filter(cold_season_dormancy_ == "Yes") %>%
#   select(source_logic_for_cold_season_dormancy) %>%
#   View(.)
# ## mostly cited
# 
# ## come back with a more systematic way of classifying maybe...
# 
# 

# #### creating pie charts for presentation #####
# library(RColorBrewer)
# 
# mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
# new_traits %>%
#   count(Realm) %>%
#   arrange(desc(n)) %>%
#   mutate(prop = n/sum(n)*100) %>%
#   mutate(prop.y = cumsum(prop) - 0.5*prop) %>%
#   ggplot(., aes(x = "", y = prop, fill = Realm)) +
#   labs(fill = "Realm") +
#   geom_bar(width = 1, stat = "identity", color = "white") +
#   coord_polar("y", start = 0 ) +
#   geom_text(aes(y = prop.y, x = 0.05, label = n), color = "white") +
#   scale_fill_manual(values = mycols, labels = c("Freshwater (n=13)", "Intertidal (n=25)", 
#                                                 "Marine (n=48)", 
#                                                 "Terrestrial (n=355)")) +
#   theme_void() +
#   theme(plot.margin = unit(c(1.5,1.5,1,5,1,5), "cm"))
# 
# realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
#   mutate(range_id = paste(species, source, sep = "_"))
# 
# both <- new_traits %>%
#   filter(limit_type == 'max, min')
# 
# source_cols <- c("darkorange3", "darkgoldenrod1", "azure4")
# 
# realized_ranges %>%
#   st_drop_geometry() %>%
#   mutate(source = ifelse(.$species %in% .$species[which(duplicated(.$species))], 
#                          "Both", 
#                          as.character(source))) %>%
#   mutate(source = factor(source, levels = c("IUCN", "GBIF", "Both"))) %>%
#   count(source) %>%
#   mutate(n = c(225, 129, 85))%>%
#   arrange(-n) %>%
#   mutate(prop = n/sum(n)*100) %>%
#   mutate(prop.y = cumsum(prop) - 0.5*prop) %>%
#   ggplot(., aes(x = "", y = prop, fill = source)) +
#   labs(fill = "Range source") +
#   geom_bar(width = 1, stat = "identity", color = "white") +
#   coord_polar("y", start = 0 ) +
#   geom_text(aes(y = prop.y, x = 0.05, label = n), color = "white") +
#   scale_fill_manual(values = source_cols, labels = c("IUCN (n=225)", "GBIF (n=129)", "Both (n=85)")) +
#   theme_void()
# 
#  
# cols <- c("grey", "darkorange3")
# 
# new_traits %>%
#   count(filled = !is.na(maximum_body_size_SVL_HBL_cm_)) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# cols <- c("goldenrod", "darkorange3", "grey")
# 
# new_traits %>%
#   mutate(dispersal_distance_category = ifelse(is.na(dispersal_distance_category), "NA",
#                                                     ifelse(str_detect(logic_source_for_dispersal_distance, "no info") 
#                                               & !is.na(dispersal_distance_category), "guess",
#                                               "known"))) %>%
#   count(filled = dispersal_distance_category)%>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
#          
# 
# cols <- c("darkorange3", "grey")
# 
# new_traits %>%
#   mutate(dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_ = ifelse(is.na(dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_), 'NA', "filled")) %>%
#   count(filled = dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_)%>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# cols <- c("darkorange3", "grey")
# 
# new_traits %>%
#   mutate(Trophic_position = ifelse(is.na(Trophic_position), 'NA', "filled")) %>%
#   count(filled = Trophic_position) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# cols <- c("darkorange3", "grey")
# 
# new_traits %>%
#   mutate(highland_specialist = ifelse(is.na(highland_specialist), 'NA', "filled")) %>%
#   count(filled = highland_specialist) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# cols <- c("darkorange3", "grey")
# 
# new_traits %>%
#   mutate(Pelagic_habitat_for_marine_species_ = ifelse(is.na(Pelagic_habitat_for_marine_species_), 'NA', "filled")) %>%
#   count(filled = Pelagic_habitat_for_marine_species_) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# cols <- c("darkorange3", "grey")
# 
# new_traits %>%
#   mutate(migratory_ = ifelse(is.na(migratory_), 'NA', "filled")) %>%
#   count(filled = migratory_) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# 
# cols <- c("goldenrod", "darkorange3", "grey")
# 
# new_traits %>%
#   mutate(hot_season_dormancy_ = ifelse(is.na(hot_season_dormancy_), NA, 
#                                        ifelse(str_detect(source_logic_for_hot_season_dormancy, "no info") 
#                                               & !is.na(hot_season_dormancy_), "guess",
#                                               "known"))) %>%
#   count(filled = hot_season_dormancy_) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 
# 
# new_traits %>%
#   mutate(cold_season_dormancy_ = ifelse(is.na(cold_season_dormancy_), NA, 
#                                        ifelse(str_detect(source_logic_for_cold_season_dormancy, "no info") 
#                                               & !is.na(cold_season_dormancy_), "guess",
#                                               "known"))) %>%
#   count(filled = cold_season_dormancy_) %>%
#   ggplot(., aes(y = n, x = c(1), fill = filled)) + 
#   geom_bar(stat = "identity", fill = cols, width = 0.01) +
#   scale_x_continuous(lim = c(0.9,1.1)) +
#   theme_classic() + coord_flip() 




### cleaning and merging depth and elevation data with rest of trait 
## read in all traits
traits_all <- read_csv("./data-processed/wrangled-traits.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

d_and_e <- read_csv("./data-raw/globtherm_traits_collated_180617_depth_and_elev_Sarah.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

## fix column names 
cols <- str_replace_all(colnames(d_and_e), " ", '_') 
cols <- str_replace_all(cols, "\\/", '_') 
colnames(d_and_e)  <- cols

## fix realm column:
d_and_e <- select(d_and_e, -Realm) %>%
  left_join(., select(traits_all, genus_species, Realm))

## make sure no names got changed:
traits_all$genus_species[which(!traits_all$genus_species %in% d_and_e$genus_species)]

## ELEV ##
e <- d_and_e %>%
  filter(Realm == 'Terrestrial')

## DEPTH ##
## make sure depths are negative 
d <- d_and_e %>%
  filter(Realm == 'Marine') %>%
  mutate(upper_depth_limit = ifelse(upper_depth_limit > 0, 
                                    (upper_depth_limit)*-1,
                                    upper_depth_limit)) %>%
  mutate(lower_depth_limit = ifelse(lower_depth_limit > 0, 
                                    (lower_depth_limit)*-1,
                                    lower_depth_limit)) 

## get rid of duplicated Mya arenaria:
d <- d[-which(d$genus_species == "Mya_arenaria" & d$Order == "Veneroida"),]

## join with rest of data:
d_and_e <- d_and_e %>%
  filter(Realm != "Marine") %>%
  rbind(., d)

## looks good! merge with rest of data:
traits <- select(traits_all, -source_logic_for_pelagic_habitat, -Pelagic_habitat_for_marine_species_) %>%
  left_join(., d_and_e, by = c("genus_species", "Realm", "data_gatherer", 
                                                "Genus", "Species", "Family", "Kingdom", 
                                                "Phylum", "Class", "Order")) %>%
  mutate(notes = ifelse(!is.na(notes.x) & !is.na(notes.y), 
                        paste(notes.x, notes.y, sep = ", "),
                        ifelse(!is.na(notes.x), notes.x,
                               ifelse(!is.na(notes.y), notes.y, NA)))) %>% # combine notes columns
  select(-notes.x, -notes.y) %>%
  dplyr::rename("pelagic_habitat_for_marine_species" = `Pelagic_habitat_(for_marine_species)`)
  
## get rid of duplicated Macoma balthica:
traits <- traits[-which(traits$genus_species == "Macoma_balthica" & traits$habitat_data_gatherer == "Sarah"),]

## write: 
write.csv(traits, "./data-processed/wrangled-traits.csv", row.names = FALSE)



##############################################################################################
##    write out new traits database with new GARD species and missed GBIF species           ##
##############################################################################################
## this version of the thermal tolerance database has all terrestrial, marine and species we have a realized range for: 
thermal_limits_allspp <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges_taxized.csv")

## read in all old traits
traits_all <- read_csv("data-raw/globtherm_traits_collated_180617.csv") %>%
  mutate(genus_species = paste(.$Genus, .$Species, sep = "_"))

## look for traits for GBIF species that were missed
missed <- unique(thermal_limits_allspp$genus_species[which(!thermal_limits_allspp$genus_species %in% traits$genus_species)])
## add ones still missing:
missed <- append(as.character(missed), c("Agabus_ramblae", "Apteropanorpa_tasmanica", "Armases_ricordi", "Brachyrhaphis_episcopi", "Dallia_pectoralis", "Deroceras_reticulatum", "Eciton_burchellii", "Etheostoma_flabellare", "Gambusia_nicaraguensis",  "Gammarus_locusta", "Lygisaurus_foliorum", "Margarites_refulgens", 'Onoba_gelida', 'Ophionereis_fasciata', 'Orchomenella_pinguides', 'Paramoera_walkeri'))
missed_traits <- traits_all[which(traits_all$genus_species %in% missed),] 
missed <- missed[which(!missed %in% missed_traits$genus_species)]
missed_traits <- filter(missed_traits, !str_detect(missed_traits$Realm, "Freshwater")) ## get rid if freshwater spp

missing_cols <- colnames(missed_traits)[which(!colnames(missed_traits) %in% colnames(traits))]
missed_traits <- missed_traits[, -which(colnames(missed_traits) %in% missing_cols[19:29])] 
missing_cols <- which(!colnames(missed_traits) %in% colnames(traits))
colnames(missed_traits)[missing_cols] <- colnames(traits)[c(3,11,12,14,15,16,17,18,19
                                                            ,20,21,22,23,25,26,28,42,29,13)]

missed_traits$maximum_body_size_SVL_HBL_cm_ <- as.numeric(missed_traits$maximum_body_size_SVL_HBL_cm_ )

## add the species in this database that ARE NOT in traits to the bottom of traits:
new_infos <- thermal_limits_allspp[which(thermal_limits_allspp$genus_species %in% missed),] 

types <- new_infos %>%
  select(genus_species, type) %>%
  unique() 

types <- aggregate(types$type, list(types$genus_species), paste, collapse = ", ") %>%
  dplyr::rename("genus_species" = Group.1, "limit_type" = x)

new_infos <- new_infos %>%
  select(-type) %>%
  filter(!duplicated(genus_species)) %>%
  left_join(., types, by = "genus_species") 

missed_therm <- thermal_limits_allspp[which(thermal_limits_allspp$genus_species 
                                            %in% missed_traits$genus_species),] %>%
  select(genus_species, type) %>%
  unique() 

missed_therm <- aggregate(missed_therm$type, list(missed_therm$genus_species),
                          paste, collapse = ", ") %>%
  dplyr::rename("genus_species" = Group.1, "limit_type" = x)


missed_traits <- left_join(missed_traits, missed_therm, by = "genus_species") %>%
  mutate(notes = ifelse(is.na(limit_type), 
                       paste(as.character(notes), ", missing from thermal limit data!", sep = ""),
                       as.character(notes)))

## add missing species:
new_species <- traits[1:nrow(new_infos),1:ncol(traits)] 
new_species[1:nrow(new_infos),] <- NA 
new_species <- new_species %>%
  mutate(genus_species = new_infos$genus_species)%>%
  mutate(Genus = new_infos$Genus, Species = new_infos$Species, Family = new_infos$Family,
         Phylum = new_infos$Phylum, Class = new_infos$Class, 
         Order = new_infos$Order, Realm = new_infos$realm, limit_type = new_infos$limit_type) 

traits_new <- rbind(traits, new_species) %>%
  filter(Realm != "Freshwater")

traits_new <- full_join(traits_new, missed_traits)

traits_new$newly_added <- ifelse(traits_new$genus_species %in% missed_traits$genus_species |
                                   traits_new$genus_species %in% new_species$genus_species, "Y", NA)

## write to fill to fill in:
write.csv(traits_new, "data-processed/rangetherm-traits_all-spp.csv", row.names = FALSE)


## read filled version of traits:
traits <- read.csv("data-raw/rangetherm-traits_all-spp_filled.csv")




## write out database for collecting Tsel data:
tpref <- traits_new %>%
  filter(Realm == "Terrestrial") %>%
  select(Genus, Species, Kingdom, Phylum, Class, Order, Family, limit_type) %>%
  mutate(Tpref = NA, ref_Tpref = NA, n_Tpref = NA, ref_from_ref = NA, notes_Tpref = NA, )

write.csv(tpref, "data-processed/Tpref.csv", row.names = FALSE)


## read in filled tpref:
tpref <- read.csv("data-raw/Tpref_filled.csv") %>%
  mutate(Tpref = as.character(Tpref))

## 1. for all Tprefs reported as a range of temperatures, midpoint of range:
tpref <- tpref %>%
  mutate(has_hyphen = ifelse(str_detect(Tpref, "\\-"), as.character(Tpref), NA)) %>%
  mutate(first = ifelse(!is.na(has_hyphen), 
                        str_split_fixed(has_hyphen, "\\-", n=2)[,1], 
                        NA)) %>%
  mutate(second = ifelse(!is.na(has_hyphen), 
                        str_split_fixed(has_hyphen, "\\-", n=2)[,2], 
                        NA)) %>%
  mutate(mean = (as.numeric(first) + as.numeric(second)) / 2)  %>%
  mutate(Tpref = ifelse(!is.na(has_hyphen), mean, as.character(Tpref))) %>%
  select(-has_hyphen, -first, -second, -mean) %>%
  mutate(genus_species = paste(Genus, Species, sep = "_")) 
  

## 2. for species with multiple measures of Tpref/Tb, take the weighted mean
## (if sample size is unknown, exclude from the weighted mean)

dups_tpref <- tpref %>%
  filter(genus_species %in% .$genus_species[duplicated(.$genus_species)]) %>%
  filter(!is.na(Tpref) & !is.na(n_tpref)) %>% ## remove unknown sample sizes
  group_by(genus_species) %>%
  mutate(ref_Tpref = paste(ref_Tpref, collapse = " "),
         ref_from_ref = paste(ref_from_ref, collapse = " "),
         notes_Tpref = paste(notes_Tpref, collapse = " ")) %>%
  summarise(Tpref = round(weighted.mean(as.numeric(Tpref), n_tpref), 2), ref_Tpref, ref_from_ref,
            notes_Tpref, limit_type, Tb, ref_Tb, notes_Tb) %>%
  ungroup() %>%
  filter(!duplicated(.))

dups_tb <- tpref %>%
  filter(genus_species %in% .$genus_species[duplicated(.$genus_species)]) %>%
  filter(!is.na(n_Tb) & !is.na(Tb)) %>% ## remove unknown sample sizes
  group_by(genus_species) %>%
  mutate(notes_Tb = paste(notes_Tb, collapse = " "))  %>%
  summarise(Tb = round(weighted.mean(as.numeric(Tb), n_Tb), 2), ref_Tb,
            notes_Tb, limit_type,  Tpref, ref_Tpref, ref_from_ref, notes_Tpref) %>%
  ungroup() %>%
  filter(!duplicated(.))

## rejoin:
tprefs_clean <- tpref %>%
  filter(!genus_species %in% .$genus_species[duplicated(.$genus_species)])%>%
  select(genus_species, Tpref, ref_Tpref, ref_from_ref, notes_Tpref, limit_type, Tb, ref_Tb, notes_Tb) %>% 
  rbind(., dups_tb, dups_tpref) %>%
  filter(!is.na(Tpref) | !is.na(Tb))

## write out:
write.csv(tprefs_clean, "data-processed/Tpref_clean.csv", row.names = FALSE)



### looking at field Tb versus lab Tpref
meiri <- read.csv("data-raw/lizard_body_temperatures_natural_history_and_life-history_traits.csv") %>%
  mutate(genus_species = str_replace_all(species, " ", "_"))

ol <- left_join(tprefs_clean, meiri, by = c("genus_species")) %>%
  mutate(Tpref = as.numeric(Tpref))

ol %>%
  ggplot(., aes(x = Tpref, y = mean.body.temperature, col = Family)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Preferred temperature (C)", y = "Mean field body temperature (C)", col = "Family")


tb <- meiri %>%
  filter(genus_species %in% traits_new$genus_species & !genus_species %in% ol$genus_species) 


one_lim <- acc_data %>%
  filter(genus_species %in% acc_data$genus_species[duplicated(acc_data$genus_species)]) %>%
  group_by(genus_species) %>%
  filter(genus_species %in% names(Te_sun_max) | realm %in% c("Marine", "Intertidal")) %>%
  filter(any(is.na(acclimation_temperature)))

both_lim <- acc_data %>%
  filter(genus_species %in% acc_data$genus_species[duplicated(acc_data$genus_species)]) %>%
  group_by(genus_species) %>%
  filter(genus_species %in% names(Te_sun_max) | realm %in% c("Marine", "Intertidal")) %>%
  filter(!any(is.na(acclimation_temperature)))

all <- thermal_limits_allspp %>%
  filter(genus_species %in% thermal_limits_allspp$genus_species[duplicated(thermal_limits_allspp$genus_species)]) %>%
  filter(genus_species %in% names(Te_sun_max) | realm %in% c("Marine", "Intertidal")) 
