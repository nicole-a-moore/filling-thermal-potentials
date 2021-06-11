# load libraries ----------------------------------------------------------
library(tidyverse)
library(broom)


#### Step 2: read in the data
globtherm_hasrealm<-read.csv(file="data-raw/GlobalTherm_18_6_with_Realm.csv", header=T)
dim(globtherm_hasrealm)


#merge published globtherm with this one that has realm
#eventually want these columns: thermy,realm, body_size, mobility
globtherm_published<-read.csv(file="data-raw/Published_GlobalTherm_upload_10_11_17.csv", header=T)
merged<-merge(globtherm_published, globtherm_hasrealm, by = c("Genus","Species"), all.x=T)

globtherm<-merged[1:46] %>%
  mutate(Realm=merged$Realm, Realm.1=merged$Realm.1, Realm.2=merged$Realm.2)


#ditch .x from merging
names(globtherm)
names(globtherm)<-gsub(".x","",names(globtherm),fixed = TRUE)

#for pretreatment, "F" = field fresh - but no info on what time of year samples collected.
#there are 1060 empty cells and 160 F's
#make all empty cells and "F" cells into NAs
#make all others numeric
#also, recategorize max_metric to only 2 categories, lethal and critical limits
globtherm<-globtherm %>%
  mutate(max_pretreatmentnum=ifelse(max_pretreatment %in% 
                                      c("", "F", "NA"), NA, as.numeric(max_pretreatment))) %>%
  mutate(max_pretreatmentnum=ifelse(min_pretreatment %in% 
                                      c("", "F", "NA"), NA, as.numeric(min_pretreatment))) 



#  mutate(max_metric_2 = as.factor(ifelse(max_metric %in% c("ctmax", "ct50", "rmax", "UTNZ"), "ct", "lt"))) %>%
#  mutate(min_metric_2 = as.factor(ifelse(globtherm$min_metric %in% c("ctmin", "LTNZ"), "ct", "lt")))
names(globtherm)

#recategorize so there is a realm  and a "thermy" for each and put birds into terrestrial endo group

globtherm<-globtherm %>%
  mutate(thermy=case_when(globtherm$Phylum %in% c("Ascomycota", 
                                        "Basidiomycota", 
                                        "Rhodophyta", 
                                        "Chlorophyta", 
                                        "Streptophyta", 
                                        "Jungermanniopsida", 
                                        "Polytrichopsida", "Phaeophyceae") ~ "plant",
                          globtherm$Class %in% c("Actinopteri", "Amphibia", "Arachnida", "Asteroidea","Bivalvia",
                         "Branchiopoda","Chondrichthyes","Collembola","Echinoidea","Gastropoda",
                         "Insecta", "Lepidosauria", "Malacostraca", "Maxillopoda", "Oligochaeta", "Ophiuroidea",
                         "Ostracoda", "Polychaeta", "Polyplacophora", "Rhynchonellata", "Pycnogonida", "Archelosauria", "Ascidiacea") ~ "ectotherm",
                         globtherm$Class %in% c("Mammalia","Aves") ~ "endotherm",
            TRUE ~ "flag"))


#make all birds terrestrial
globtherm<-globtherm %>% 
  mutate(realm = replace(Realm, Class == "Aves", "Terrestrial"))


#make all marine and intertidal elevations = 0
temp_marine<-globtherm %>%
  filter(realm %in% c("Marine", "Intertidal")) %>%
  mutate(elevation_min_0=0)%>%
  mutate(elevation_max_0=0)


globtherm<-rbind(globtherm %>%
                   filter(!realm %in% c("Marine", "Intertidal")) %>%
                   mutate(elevation_min_0=elevation_min) %>%
                   mutate(elevation_max_0=elevation_max),
                 temp_marine)
#check
globtherm %>%
  filter(globtherm$realm %in% c("Marine", "Intertidal")) %>%
  select(elevation_min_0, elevation_max_0)


write_csv(globtherm, "data-processed/globtherm_processed.csv")

