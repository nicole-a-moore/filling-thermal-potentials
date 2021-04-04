## for  cleaning and manipulating the traits database 
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
  rename("genus_species" = Group.1, "limit_type" = x)

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
  rename("genus_species" = Group.1, "limit_type" = x)

traits_sub <- left_join(traits_sub, lim_types, by = "genus_species")


## write out and start filling in the missing ones!!
write.csv(traits_sub, "data-processed/globtherm_traits_collated_180617_ectotherms-with-limits.csv", row.names = FALSE)




###### inspecting dormancy trait completeness #####
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

## check how many are yes and no
length(which(dormancy_sub$cold_season_dormancy == "Yes")) ## 39/141 
length(which(dormancy_sub$hot_season_dormancy == "Yes")) ## 0/141

## check how many of those are terrestrial
terr_sub <- dormancy_sub %>%
  filter(Realm == "Terrestrial") 

length(which(terr_sub$cold_season_dormancy == "Yes")) ## 39/112, all of them!





## plot dormancy Y/N against realized range latitudinal midpoint (should end up with 141 species)
realized_ranges <- st_read("data-processed/realized-ranges_split.shp")

lat_mps <- realized_ranges %>%
  filter(!duplicated(species)) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  filter(species %in% dormancy_sub$genus_species) %>%  ## have 140 since 1 species is duplicated in traits data... oops
  as.data.frame(.) %>%
  select(species, lttdnl_) %>%
  rename("lat_mp" = lttdnl_)
  
dormancy_sub <- left_join(dormancy_sub, lat_mps, by = c("genus_species" = "species"))

dormancy_sub %>%
  filter(Realm == "Terrestrial") %>%
  ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy, col = Order)) + geom_point() ## very interesting.... seems like they are all squamata

unique(dormancy_sub$Order[which(dormancy_sub$cold_season_dormancy == "Yes")]) ## it is true!

## now look at only squamata: 
dormancy_sub %>%
  filter(Order == "Squamata") %>% ## 41/111 squamata species are dormant in the cold seasons
  ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy, col = Order)) + geom_point()

## print out names 
dormancy_sub %>%
  filter(Order == "Squamata")  %>%
  filter(cold_season_dormancy == "Yes") %>%
  .$genus_species




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
                                                   "Yes", hot_season_dormancy_)))
## write to file:
write.csv(dormancy_sub, "./data-processed/wrangled-traits.csv", row.names = FALSE)

dormancy_sub <- dormancy_sub %>%
  filter(limit_type == "max, min") %>%
  filter(!is.na(cold_season_dormancy_) & !is.na(hot_season_dormancy_)) 

unique(dormancy_sub$cold_season_dormancy_)
unique(dormancy_sub$hot_season_dormancy_)

## check how many are yes and no
length(which(dormancy_sub$cold_season_dormancy_ == "Yes")) ## 119/186 
length(which(dormancy_sub$hot_season_dormancy_ == "Yes")) ## 5/186

lat_mps <- realized_ranges %>%
  filter(!duplicated(species)) %>%
  filter(species %in% dormancy_sub$genus_species) %>%  ## have 185 since 1 species is duplicated in traits data... oops
  as.data.frame(.) %>%
  select(species, lttdnl_) %>%
  rename("lat_mp" = lttdnl_)

dormancy_sub <- left_join(dormancy_sub, lat_mps, by = c("genus_species" = "species"))

dormancy_sub %>%
  filter(Realm == "Terrestrial") %>%
  ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy_, col = Order)) + geom_point() 

unique(dormancy_sub$Order[which(dormancy_sub$cold_season_dormancy == "Yes")]) ## no longer all squamata!

## now look at only squamata: 
dormancy_sub %>%
  filter(Order == "Squamata") %>%
  ggplot(., aes(x = abs(lat_mp), y = cold_season_dormancy_, col = Order)) + geom_point()
## more higher latitude species have "Yes" after updating 


## look at logic for 'No' ... did we guess these based on climate?
dormancy_sub %>%
  filter(Realm == "Terrestrial") %>%
  filter(cold_season_dormancy_ == "No") %>%
  select(source_logic_for_cold_season_dormancy) %>%
  View(.)
## about half have a cited source... leave for now


## look at logic for 'Yes'... did we guess these based on climate?
dormancy_sub %>%
  filter(Realm == "Terrestrial") %>%
  filter(cold_season_dormancy_ == "Yes") %>%
  select(source_logic_for_cold_season_dormancy) %>%
  View(.)
## mostly cited

## come back with a more systematic way of classifying maybe...



#### creating pie charts for presentation #####
library(RColorBrewer)

mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
new_traits %>%
  count(Realm) %>%
  arrange(desc(n)) %>%
  mutate(prop = n/sum(n)*100) %>%
  mutate(prop.y = cumsum(prop) - 0.5*prop) %>%
  ggplot(., aes(x = "", y = prop, fill = Realm)) +
  labs(fill = "Realm") +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0 ) +
  geom_text(aes(y = prop.y, x = 0.05, label = n), color = "white") +
  scale_fill_manual(values = mycols, labels = c("Freshwater (n=13)", "Intertidal (n=25)", 
                                                "Marine (n=48)", 
                                                "Terrestrial (n=355)")) +
  theme_void() +
  theme(plot.margin = unit(c(1.5,1.5,1,5,1,5), "cm"))

realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
  mutate(range_id = paste(species, source, sep = "_"))

both <- new_traits %>%
  filter(limit_type == 'max, min')

source_cols <- c("darkorange3", "darkgoldenrod1", "azure4")

realized_ranges %>%
  st_drop_geometry() %>%
  mutate(source = ifelse(.$species %in% .$species[which(duplicated(.$species))], 
                         "Both", 
                         as.character(source))) %>%
  mutate(source = factor(source, levels = c("IUCN", "GBIF", "Both"))) %>%
  count(source) %>%
  mutate(n = c(225, 129, 85))%>%
  arrange(-n) %>%
  mutate(prop = n/sum(n)*100) %>%
  mutate(prop.y = cumsum(prop) - 0.5*prop) %>%
  ggplot(., aes(x = "", y = prop, fill = source)) +
  labs(fill = "Range source") +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0 ) +
  geom_text(aes(y = prop.y, x = 0.05, label = n), color = "white") +
  scale_fill_manual(values = source_cols, labels = c("IUCN (n=225)", "GBIF (n=129)", "Both (n=85)")) +
  theme_void()

 
cols <- c("grey", "darkorange3")

new_traits %>%
  count(filled = !is.na(maximum_body_size_SVL_HBL_cm_)) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

cols <- c("goldenrod", "darkorange3", "grey")

new_traits %>%
  mutate(dispersal_distance_category = ifelse(is.na(dispersal_distance_category), "NA",
                                                    ifelse(str_detect(logic_source_for_dispersal_distance, "no info") 
                                              & !is.na(dispersal_distance_category), "guess",
                                              "known"))) %>%
  count(filled = dispersal_distance_category)%>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 
         

cols <- c("darkorange3", "grey")

new_traits %>%
  mutate(dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_ = ifelse(is.na(dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_), 'NA', "filled")) %>%
  count(filled = dispersal_type_walking_may_need_to_be_reworded_to_something_that_encaptures_slithering_)%>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

cols <- c("darkorange3", "grey")

new_traits %>%
  mutate(Trophic_position = ifelse(is.na(Trophic_position), 'NA', "filled")) %>%
  count(filled = Trophic_position) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

cols <- c("darkorange3", "grey")

new_traits %>%
  mutate(highland_specialist = ifelse(is.na(highland_specialist), 'NA', "filled")) %>%
  count(filled = highland_specialist) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

cols <- c("darkorange3", "grey")

new_traits %>%
  mutate(Pelagic_habitat_for_marine_species_ = ifelse(is.na(Pelagic_habitat_for_marine_species_), 'NA', "filled")) %>%
  count(filled = Pelagic_habitat_for_marine_species_) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

cols <- c("darkorange3", "grey")

new_traits %>%
  mutate(migratory_ = ifelse(is.na(migratory_), 'NA', "filled")) %>%
  count(filled = migratory_) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 


cols <- c("goldenrod", "darkorange3", "grey")

new_traits %>%
  mutate(hot_season_dormancy_ = ifelse(is.na(hot_season_dormancy_), NA, 
                                       ifelse(str_detect(source_logic_for_hot_season_dormancy, "no info") 
                                              & !is.na(hot_season_dormancy_), "guess",
                                              "known"))) %>%
  count(filled = hot_season_dormancy_) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 

new_traits %>%
  mutate(cold_season_dormancy_ = ifelse(is.na(cold_season_dormancy_), NA, 
                                       ifelse(str_detect(source_logic_for_cold_season_dormancy, "no info") 
                                              & !is.na(cold_season_dormancy_), "guess",
                                              "known"))) %>%
  count(filled = cold_season_dormancy_) %>%
  ggplot(., aes(y = n, x = c(1), fill = filled)) + 
  geom_bar(stat = "identity", fill = cols, width = 0.01) +
  scale_x_continuous(lim = c(0.9,1.1)) +
  theme_classic() + coord_flip() 
