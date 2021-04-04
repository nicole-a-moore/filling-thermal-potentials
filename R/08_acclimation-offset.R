## calculating acclimation offset for each species 
library(tidyverse)
library(broom)

## check how many species we have are in intratherm:
intratherm <- read.csv("/Users/nikkimoore/Documents/intra-therm/data-processed/intratherm-with-elev.csv")
traits <- read.csv("data-processed/wrangled-traits.csv")

intra_sp <- unique(intratherm$genus_species)
tp_sp <- unique(traits$genus_species)

overlap <- tp_sp[tp_sp %in% intra_sp]
length(overlap) ## 36 of our species are in intratherm, most only have tmax estimates

## plot those ARRs
intratherm %>%
  filter(genus_species %in% overlap, !is.na(acclim_temp), !is.na(parameter_value), 
         genus_species %in% genus_species[which(duplicated(genus_species))]) %>%
  filter(parameter_tmax_or_tmin == 'tmax') %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")

intratherm %>%
  filter(genus_species %in% overlap, !is.na(acclim_temp), !is.na(parameter_value), 
         genus_species %in% genus_species[which(duplicated(genus_species))]) %>%
  filter(parameter_tmax_or_tmin == 'tmin') %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")

## fit models:
ol_arr <- intratherm %>%
  mutate(acclim_temp = as.numeric(as.character(acclim_temp))) %>%
  filter(genus_species %in% overlap, !is.na(acclim_temp), !is.na(parameter_value), 
         genus_species %in% genus_species[which(duplicated(genus_species))]) %>%
  group_by(genus_species) %>%
  do(tidy(lm(parameter_value ~ acclim_temp, data = .), conf.int = TRUE))


## investigate how many tmax vs tmin estimates in intratherm
intratherm %>%
  count(parameter_tmax_or_tmin)
## 2632 tmax, 158 tmin

intratherm %>%
  filter(parameter_tmax_or_tmin == "tmin") %>%
  select(acclim_temp) %>%
  unique(.) %>%
  nrow(.) ## 27 acclim temps for tmin

intratherm %>%
  filter(parameter_tmax_or_tmin == "tmax") %>%
  select(acclim_temp) %>%
  unique(.) %>%
  nrow(.) ## 185 acclim temps for tmax


## fit model to find global ARR?
