## plotting range filling metrics
library(tidyverse)
library(sf)
library(gridExtra)


metrics <- read.csv("data-processed/range-filling-quantifications.csv") 

thermal_limits <- read_csv("data-raw/globtherm_full_dataset_2019.csv") %>%
  filter(thermy == "ectotherm") %>%
  mutate(genus_species = str_replace_all(.$genus_species, "_", " ")) 

realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") 

rr <- realized_ranges %>%
  mutate(realized_range_area = units::set_units(st_area(.), km^2)) %>%
  as.data.frame(.) %>%
  mutate(range = paste(str_replace_all(species, " ", "."), source, sep = "_"))

data <- left_join(metrics, thermal_limits, by = c("species" = "genus_species")) %>%
  left_join(., rr, by = c("range","species", "source"))


no_duplicates <- data %>%
  select(-extreme_enviro_temp, -acclimation_offset, -thermal_limit, -metric, -type, 
         -collection_latitude, - absolute_latitude, -elevation, -collection_longitude, 
         -start_temp, -acclimation_temperature, -ramping_rate, -exposure_duration) %>%
  filter(!duplicated(.))





## plot across latitude, across realized and potential range size?, across realm
h1 <- no_duplicates %>%
  ggplot(., aes(x = log(equ_of + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(equatorward overfilling + 1)")

h2 <- no_duplicates %>%
  ggplot(., aes(x = log(pol_of + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(poleward overfilling + 1)")

h3 <- no_duplicates %>%
  ggplot(., aes(x = log(equ_op + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(equatorward overpredicting + 1)")

h4 <- no_duplicates %>%
  ggplot(., aes(x = log(pol_op + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(poleward overpredicting + 1)")

g1 <- grid.arrange(h1, h2, h3, h4, nrow = 2)

ggsave(g1, path = "./figures/", filename = "histogram_unweighted-metrics.png", device = "png", 
       height = 6, width = 8)


hw1 <- no_duplicates %>%
  ggplot(., aes(x = log(weighted_equ_of + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(weighted equatorward overfilling + 1)")

hw2 <- no_duplicates %>%
  ggplot(., aes(x = log(weighted_pol_of + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(weighted poleward overfilling + 1)")

hw3 <- no_duplicates %>%
  ggplot(., aes(x = log(weighted_equ_op + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(weighted equatorward overpredicting + 1)")

hw4 <- no_duplicates %>%
  ggplot(., aes(x = log(weighted_pol_op + 1))) + geom_histogram() +
  labs(y = "Frequency", x = "ln(weighted poleward overpredicting + 1)")

gw1 <- grid.arrange(hw1, hw2, hw3, hw4, nrow = 2)

ggsave(gw1, path = "./figures/", filename = "histogram_weighted-metrics.png", device = "png", 
       height = 6, width = 8)


## proportion now:
hp1 <- no_duplicates %>%
  ggplot(., aes(x = equ_of_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Proportion equatorward overfilling")

hp2 <- no_duplicates %>%
  ggplot(., aes(x = pol_of_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Proportion poleward overfilling")

hp3 <- no_duplicates %>%
  ggplot(., aes(x = equ_op_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Proportion equatorward overpredicting")

hp4 <- no_duplicates %>%
  ggplot(., aes(x = pol_op_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Proportion poleward overpredicting")

gp1 <- grid.arrange(hp1, hp2, hp3, hp4, nrow = 2)

ggsave(gp1, path = "./figures/", filename = "histogram_proportion-metrics.png", device = "png", 
       height = 6, width = 8)



hpw1 <- no_duplicates %>%
  ggplot(., aes(x = weighted_equ_of_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Weighted proportion equatorward overfilling")

hpw2 <- no_duplicates %>%
  ggplot(., aes(x = weighted_pol_of_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Weighted proportion poleward overfilling")

hpw3 <- no_duplicates %>%
  ggplot(., aes(x = weighted_equ_op_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Weighted proportion equatorward overpredicting")

hpw4 <- no_duplicates %>%
  ggplot(., aes(x = weighted_pol_op_prop)) + geom_histogram() +
  labs(y = "Frequency", x = "Weighted proportion poleward overpredicting")

gpw1 <- grid.arrange(hpw1, hpw2, hpw3, hpw4, nrow = 2)

ggsave(gpw1, path = "./figures/", filename = "histogram_weighted-proportion-metrics.png", device = "png", 
       height = 6, width = 8)




## plot across realized range latitudinal midpoint:
across_lat <- no_duplicates %>%
  gather(key = "of_op_metric", value = "value", c(equ_of_prop, equ_op_prop, pol_op_prop, pol_of_prop)) %>%
  ggplot(., aes(x = abs(lat_mp), y = value, col = of_op_metric)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "|latitduinal midpoint of realized range|", y = "proportion range filling metric", colour = "metric:")

ggsave(across_lat, path = "./figures/", filename = "across-latitude-unweighted-prop.png", device = "png", 
       height = 6, width = 8)

across_lat_w <- no_duplicates %>%
  gather(key = "of_op_metric", value = "value", c(weighted_equ_of_prop, weighted_equ_op_prop, weighted_pol_op_prop, weighted_pol_of_prop)) %>%
  ggplot(., aes(x = abs(lat_mp), y = value, col = of_op_metric)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "|latitduinal midpoint of realized range|", y = "weighted proportion range filling metric", colour = "metric:")

ggsave(across_lat_w, path = "./figures/", filename = "across-latitude-weighted-prop.png", device = "png", 
       height = 6, width = 8)




## across range size:
no_duplicates %>%
  filter(as.numeric(realized_range_area) < max(as.numeric(realized_range_area))) %>%
  gather(key = "of_op_metric", value = "value",  c(equ_of, equ_op, pol_op, pol_of)) %>%
  ggplot(., aes(x = as.numeric(realized_range_area), y = log(value), col = of_op_metric)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "area of realized range (km^2)", y = "ln(range filling metric)", colour = "metric:")

no_duplicates %>%
  filter(as.numeric(realized_range_area) < max(as.numeric(realized_range_area))) %>%
  gather(key = "of_op_metric", value = "value",  c(equ_of_prop, equ_op_prop, pol_op_prop, pol_of_prop)) %>%
  ggplot(., aes(x = as.numeric(realized_range_area), y = value, col = of_op_metric)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "area of realized range (km^2)", y = "proportion range filling metric", colour = "metric:")


no_duplicates %>%
  filter(as.numeric(realized_range_area) < max(as.numeric(realized_range_area))) %>%
  gather(key = "of_op_metric", value = "value", c(weighted_equ_of_prop, 
                                                  weighted_pol_of_prop, 
                                                  weighted_equ_op_prop,
                                                  weighted_pol_op_prop)) %>%
  ggplot(., aes(x = as.numeric(realized_range_area), y = value, col = of_op_metric)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "area of realized range (km^2)", y = "weighted proportion range filling metric", colour = "metric:")



## Q-Q plots: