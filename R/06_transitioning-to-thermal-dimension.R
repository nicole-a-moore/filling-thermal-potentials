## transitioning into the thermal dimension
library(tidyverse)
library(raster)
library(sf)
library(gridExtra)
select <- dplyr::select

## goal: to measure thermal range overfilling and underfilling in the thermal dimension

## this requires:
## 1. geting max and minimum temperatutes in the raster squares under each realized range (maybe need average values later?)
##    a. combine marine and terrestrial temperatute data into one layer 

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## read in seasonal high and low temp data:
terr_seasonal_high <- read.csv("data-processed/terrestrial_seasonal-max-temps_6mo-dormancy.csv")
terr_seasonal_low <- read.csv("data-processed/terrestrial_seasonal-min-temps_6mo-dormancy.csv")
marine_seasonal_high <- read.csv("data-processed/marine_seasonal-max-temps_6mo.csv") 
marine_seasonal_low <- read.csv("data-processed/marine_seasonal-min-temps_6mo.csv")

## rasterize:
raster_terr_high <- rasterize(terr_seasonal_high[, 1:2], r, terr_seasonal_high[,3:4], fun=mean) 
raster_terr_high[is.infinite(raster_terr_high)] <- NA
names(raster_terr_high) <- c("seasonal_high_temp", "hot_dormancy_6mo")
##plot(raster_terr_high, asp = 1)
raster_marine_high <- rasterize(marine_seasonal_high[, 1:2], r, marine_seasonal_high[,3:4], fun=mean)
raster_marine_high[is.infinite(raster_marine_high)] <- NA
names(raster_marine_high) <- c("seasonal_high_temp", "hot_dormancy_6mo")
##plot(raster_marine_high, asp = 1)

raster_terr_low <- rasterize(terr_seasonal_low[, 1:2], r, terr_seasonal_low[,3:4], fun=mean)
raster_terr_low[is.infinite(raster_terr_low)] <- NA
names(raster_terr_low) <- c("seasonal_low_temp", "cold_dormancy_6mo")
##plot(raster_terr_low, asp = 1)
raster_marine_low <- rasterize(marine_seasonal_low[, 1:2], r, marine_seasonal_low[,3:4], fun=mean)
raster_marine_low[is.infinite(raster_marine_low)] <- NA
names(raster_marine_low) <-  c("seasonal_low_temp", "cold_dormancy_6mo")
##plot(raster_marine_low, asp = 1)

## merge marine high and low temps into a single raster layer:
## if cells overlap at the land-ocean boundary, keep the sea temperatures 
high_temps <- merge(raster_marine_high, raster_terr_high)
names(high_temps) <- c("seasonal_high_temp", "hot_dormancy_6mo")
##plot(high_temps[[1]], asp = 1)
low_temps <- merge(raster_marine_low, raster_terr_low)
names(low_temps) <- c("seasonal_low_temp", "cold_dormancy_6mo")
##plot(low_temps, asp = 1)


## plot anomalies between extreme and extreme w dormancy temps:
low <- low_temps[[2]] - low_temps[[1]]

high <- high_temps[[1]] - high_temps[[2]]

low_df <- data.frame(rasterToPoints(low)) 
names(low_df) <- c("x", "y", "temp_dif")

l_plot <- low_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = temp_dif)) + 
  coord_fixed() +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", 
       fill = "|Temperature \ndifference|") +
  scale_fill_gradient(low = "white", high = "cornflowerblue",
                      na.value = "black") +
  scale_x_continuous(limits = c(-180, 180),
                     expand = c(0,0),
                     breaks = c(-150, -120, -90, 
                                -60, -30, 0, 30,
                                60, 90, 120, 150), 
                     labels = c('150°W', '120°W', '90°W', 
                                '60°W','30°W', '0°', 
                                '30°E', '60°E', '90°E', 
                                '120°E', '150°E')) +
  scale_y_continuous(limits = c(-90, 90),
                     expand = c(0,0),
                     breaks = c(-90, -60, -30, 0, 30, 60, 90),
                     labels = c('90°S', '60°S', '30°S', '0°'
                                , '30°N', '60°N', '90°N')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank()
  )

high_df <- data.frame(rasterToPoints(high)) 
names(high_df) <- c("x", "y", "temp_dif")

h_plot <- high_df %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = temp_dif)) + 
  coord_fixed() +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", 
       fill = "|Temperature \ndifference|") +
  scale_fill_gradient(low = "white", high = "brown1",
                      na.value = "black") +
  scale_x_continuous(limits = c(-180, 180),
                     expand = c(0,0),
                     breaks = c(-150, -120, -90, 
                                -60, -30, 0, 30,
                                60, 90, 120, 150), 
                     labels = c('150°W', '120°W', '90°W', 
                                '60°W','30°W', '0°', 
                                '30°E', '60°E', '90°E', 
                                '120°E', '150°E')) +
  scale_y_continuous(limits = c(-90, 90),
                     expand = c(0,0),
                     breaks = c(-90, -60, -30, 0, 30, 60, 90),
                     labels = c('90°S', '60°S', '30°S', '0°'
                                , '30°N', '60°N', '90°N')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank()
  )
  
ggsave(l_plot, path = "figures/",
       file = "experienced-temp-anomalies_cold.png", 
       height = 4, width = 9, units = "in", device = "png")

ggsave(h_plot, path = "figures/",
       file = "experienced-temp-anomalies_hot.png", 
       height = 4, width = 9, units = "in", device = "png")





## b. extract temperatures under each realized range for each species 
## read in rasterized realized ranges
rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")
traits <- read.csv("./data-processed/wrangled-traits.csv")

## loop through ranges and extract temperatures under each 
## if species is dormant, extract 'experienced temperatures' (ie. with temperatures when presumably dormant removed)
i = 1
realized_temps <- c()
while (i < nlayers(rasterized_rrs) + 1) {
  range <- rasterized_rrs[[i]]
  
  ## use realized range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps[[1]], range)
  highs_dormancy <- mask(high_temps[[2]], range)
  
  lows <- mask(low_temps[[1]], range)
  lows_dormancy <- mask(low_temps[[2]], range)
  
  ## extract temperatures assuming no dormancy:
  if (length(which(is.na(values(highs)) == FALSE)) > 0) {
    high_vals <- data.frame(high_or_low = "high", temps = values(highs)[which(is.na(values(highs)) == FALSE)],
                            type = "realized", 
                            dormancy = "none")
    ## plot(highs)
  }
  else {
    high_vals <- data.frame(high_or_low = "high", temps = NA, type = "realized", dormancy = "none")
  }
  if (length(which(is.na(values(lows)) == FALSE)) > 0) {
    low_vals <- data.frame(high_or_low = "low", temps = values(lows)[which(is.na(values(lows)) == FALSE)],
                           type = "realized",
                           dormancy = "none")
    ## plot(lows)
  }
  else {
    low_vals <- data.frame(high_or_low = "low", temps = NA,type = "realized", dormancy = "none")
  }
  
  ## see if species goes dormant:
  species <- str_split_fixed(names(range), pattern = "\\_", n = 2)[1,1] %>%
    str_replace_all("\\.", " ")

  species_traits <- traits[which(traits$genus_species == species),]
  cold_dormancy <- species_traits$cold_season_dormancy_ == "Yes"
  hot_dormancy <- species_traits$hot_season_dormancy_ == "Yes"
  
  ## if the species goes dormant, must extract temperatures considering and not considering dormancy 
  ## if hot and cold dormant, use both dormancy temperatures:
  if(!is.na(hot_dormancy) && hot_dormancy == TRUE &&
     !is.na(cold_dormancy) && cold_dormancy == TRUE &&
     length(which(is.na(values(highs_dormancy)) == FALSE)) > 0 &
     length(which(is.na(values(lows_dormancy)) == FALSE)) > 0) {
    
    high_vals$dormancy <- "both"
    low_vals$dormancy <- "both"
    
    high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                             temps = values(highs_dormancy)[which(is.na(values(highs_dormancy)) == FALSE)],
                                             type = "realized_dormancy",
                                             dormancy = "both"))
    
    low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                                           temps = values(lows_dormancy)[which(is.na(values(lows_dormancy)) == FALSE)],
                                           type = "realized_dormancy",
                                           dormancy = "both"))
  }
  else if (!is.na(hot_dormancy) && hot_dormancy == TRUE &&
           length(which(is.na(values(highs_dormancy)) == FALSE)) > 0) {
    
    high_vals$dormancy <- "hot"
    low_vals$dormancy <- "hot"
    
    high_vals <- rbind(high_vals, data.frame(high_or_low = "high", 
                                             temps = values(highs_dormancy)[which(is.na(values(highs_dormancy)) == FALSE)],
                                             type = "realized_dormancy",
                                             dormancy = "hot"))
  }
  else if(!is.na(cold_dormancy) && cold_dormancy == TRUE &&
          length(which(is.na(values(lows_dormancy)) == FALSE)) > 0) {
    
    high_vals$dormancy <- "cold"
    low_vals$dormancy <- "cold"
    
    low_vals <- rbind(low_vals, data.frame(high_or_low = "low", 
                                           temps = values(lows_dormancy)[which(is.na(values(lows_dormancy)) == FALSE)],
                                           type = "realized_dormancy",
                                           dormancy = "cold"))
  }
  else if (!is.na(cold_dormancy) && cold_dormancy == TRUE &&
           length(which(is.na(values(lows_dormancy)) == FALSE)) <= 0) {
    high_vals$dormancy <- "cold"
    low_vals$dormancy <- "cold"
  }
  else if (!is.na(hot_dormancy) && hot_dormancy == TRUE &&
           length(which(is.na(values(highs_dormancy)) == FALSE)) <= 0) {
    high_vals$dormancy <- "hot"
    low_vals$dormancy <- "hot"
  }
    
  if(i == 1) {
    realized_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(rasterized_rrs)[i]) %>%
      select(range, everything())
  }
  else {
    realized_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = names(rasterized_rrs)[i]) %>%
      select(range, everything()) %>%
      rbind(realized_temps, .)
  }
  
  i = i + 1
}

## save:
write.csv(realized_temps, "./data-processed/thermal-dimension_realized-temps_dormancy.csv", row.names = FALSE)
realized_temps <- read.csv("./data-processed/thermal-dimension_realized-temps_dormancy.csv") 

## 2. comparing these temperatures to thermal tolerance limits in various ways 
##    a. for species with both thermal tolerance limits, get the temperatures 
##        in their potential range using the same method

## read in potential ranges:
potential_ranges <- readRDS("data-processed/potential_ranges_notcutatequator_dormancy.rds")

## loop through ranges and extract temperatures under each 
## if species is dormant, extract 'experienced temperatures' (ie. with temperatures when presumably dormant removed)
i = 1
potential_temps <- c()
while (i < nlayers(potential_ranges) + 1) {
  range <- potential_ranges[[i]]
  
  ## use potential range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps[[1]], range)
  highs_dormancy <- mask(high_temps[[2]], range)
  
  lows <- mask(low_temps[[1]], range)
  lows_dormancy <- mask(low_temps[[2]], range)
  
  ## extract temperatures assuming no dormancy:
  if (length(which(is.na(values(highs)) == FALSE)) > 0) {
    high_vals <- data.frame(high_or_low = "high", 
                            temps = values(highs)[which(is.na(values(highs)) == FALSE)],
                            type = "potential",
                            dormancy = "none")
    ## plot(highs)
  }
  else {
    high_vals <- data.frame(high_or_low = "high", temps = NA, type = "potential", 
                            dormancy = "none")
  }
  if (length(which(is.na(values(lows)) == FALSE)) > 0) {
    low_vals <- data.frame(high_or_low = "low", temps = values(lows)[which(is.na(values(lows)) == FALSE)],
                           type = "potential", 
                           dormancy = "none")
    ## plot(lows)
  }
  else {
    low_vals <- data.frame(high_or_low = "low", temps = NA, type = "potential", 
                           dormancy = "none")
  }
  
  ## see if species goes dormant: 
  species <- names(range)
  dormant = str_detect(species, "dormancy") && !str_detect(species, "no_dormancy")
  
  ## if the species goes dormant, must extract temperatures considering and not considering dormancy 
  if(dormant) {
    split <- str_split_fixed(species, "\\_", n = 7)
    source <- split[1,7]
    species <- paste(split[1,1], split[1,2], sep = ".")
    dormancy_type <- split[,5]
    which_pr <- split[,3]
    range_id <- paste(species, source, sep = "_")
    
    if (which_pr == "0") {
      high_vals$dormancy = dormancy_type
      low_vals$dormancy = dormancy_type
    }
    
    ## if hot dormant, use hot dormancy temperatures:
    if (length(which(is.na(values(highs_dormancy)) == FALSE)) > 0 &
        length(which(is.na(values(lows_dormancy)) == FALSE)) > 0 &
        which_pr == "6" &&
        dormancy_type == "both") {
      high_vals <- data.frame(high_or_low = "high", 
                                               temps = values(highs_dormancy)[which(is.na(values(highs_dormancy)) == FALSE)],
                                               type = "potential_dormancy",
                                               dormancy = "both")
      
      low_vals <- data.frame(high_or_low = "low", 
                                             temps = values(lows_dormancy)[which(is.na(values(lows_dormancy)) == FALSE)],
                                             type = "potential_dormancy",
                                             dormancy = "both")
      
      ## plot(highs_dormancy)
      ## plot(lows_dormancy)
    }
    else if (length(which(is.na(values(highs_dormancy)) == FALSE)) > 0 & 
             which_pr == "6" &&
             dormancy_type == "hot"){
      high_vals <- data.frame(high_or_low = "high", 
                                               temps = values(highs_dormancy)[which(is.na(values(highs_dormancy)) == FALSE)],
                                               type = "potential_dormancy",
                                               dormancy = "hot")
      
      if (length(which(is.na(values(lows_dormancy)) == FALSE)) > 0) {
        low_vals <- data.frame(high_or_low = "low", 
                               temps = values(lows)[which(is.na(values(lows)) == FALSE)],
                               type = "potential_dormancy",
                               dormancy = "hot")
      }
      else {
        low_vals$dormancy = "hot"
      }
      
    }
    else if (length(which(is.na(values(lows_dormancy)) == FALSE)) > 0 & 
             which_pr == "6" &&
             dormancy_type == "cold") {
      low_vals <- data.frame(high_or_low = "low", 
                                             temps = values(lows_dormancy)[which(is.na(values(lows_dormancy)) == FALSE)],
                                             type = "potential_dormancy",
                                             dormancy = "cold")
      if (length(which(is.na(values(highs_dormancy)) == FALSE)) > 0) {
        high_vals <- data.frame(high_or_low = "high", 
                                temps = values(highs)[which(is.na(values(highs)) == FALSE)],
                                type = "potential_dormancy",
                                dormancy = "cold")
        
      }
      else {
        high_vals$dormancy = dormancy_type
      }
      
    }
  }
  else {
    if (str_detect(species, "no_dormancy")) {
      split <- str_split_fixed(species, "_", n = 5)
      range_id <- paste(paste(split[1,1], split[1,2], sep = "."), split[1,5], sep = "_")
    }
    else {
      split <- str_split_fixed(species, "_", n = 3)
      range_id <- paste(paste(split[1,1], split[1,2], sep = "."), split[1,3], sep = "_")
    }
    
  }
  
  if(i == 1) {
    potential_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = range_id) %>%
      select(range, everything())
  }
  else {
    potential_temps <- rbind(high_vals, low_vals) %>%
      mutate(range = range_id) %>%
      select(range, everything()) %>%
      rbind(potential_temps, .)
  }
  
  i = i + 1
}

## save:
write.csv(potential_temps, "./data-processed/thermal-dimension_potential-temps_dormancy.csv", row.names = FALSE)
#potential_temps <- read.csv("./data-processed/thermal-dimension_potential-temps_dormancy.csv")


## plot them all: 
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = "."))

all_temps <- rbind(realized_temps, potential_temps)
num = 1 
while (num < length(unique(potential_temps$range)) + 1) {
  
  h <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "high")
  
  l <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "low") 
  
  ## get thermal limits:
  lims <- thermal_limits[which(thermal_limits$genus_species == str_split_fixed(h$range[1], "_", n = 2)[1,1]),]
  ctmax <- lims$thermal_limit[which(lims$type == "max")]
  ctmin <- lims$thermal_limit[which(lims$type == "min")]
  
  hplot <- h %>%
    ggplot(., aes(x = temps, fill = type)) +
    labs(fill = "Range:", y = "Frequency", x = "Extreme high temperature", 
         title = paste("(", unique(l$dormancy), " season dormancy)", sep = "")) + 
    geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) +
    scale_fill_manual(values = c("yellow", "darkgrey"))
  
  lplot <- l %>%
    ggplot(., aes(x = temps, fill = type)) +
    labs(fill = "Range:", y = "Frequency", x = "Extreme low temperature", 
         title = h$range[1]) + 
    geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) +
    scale_fill_manual(values = c("yellow", "darkgrey"))
  
  hlplot <- grid.arrange(lplot, hplot, ncol = 2)
  
  ggsave(hlplot, path = "figures/thermal-dimension/", 
         filename = paste(l$range[1], ".png",sep = ""), 
         height = 5, width = 11, units = "in", device = "png")
  
  num = num + 1
}



## plot potential and realized thermal niches with and without dormancy to see how they change:
num = 1 
while (num < length(unique(potential_temps$range)) + 1) {
  
  h <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "high")
  
  l <- all_temps %>%
    filter(range == unique(potential_temps$range)[num], high_or_low == "low") 
  
  ## get thermal limits:
  lims <- thermal_limits[which(thermal_limits$genus_species == str_split_fixed(h$range[1], "_", n = 2)[1,1]),]
  ctmax <- lims$thermal_limit[which(lims$type == "max")]
  ctmin <- lims$thermal_limit[which(lims$type == "min")]
  
  if (unique(h$dormancy) != "none") {
    c_dormancy_plot <- l %>%
      filter(type %in% c("potential_dormancy", "realized_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal low temperature") + 
      geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    c_no_dormancy_plot <- l %>%
      filter(type %in% c("potential", "realized")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal low temperature") + 
      geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    h_dormancy_plot <- h %>%
      filter(type %in% c("potential_dormancy", "realized_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal high temperature", 
           title = paste("(", unique(h$dormancy), " season dormancy)", sep = "")) + 
      geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    h_no_dormancy_plot <- h %>%
      filter(type %in% c("potential", "realized")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal high temperature", 
           title = unique(h$range)) + 
      geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    prplot <- grid.arrange(h_no_dormancy_plot, h_dormancy_plot, 
                           c_no_dormancy_plot, c_dormancy_plot, 
                           ncol = 2)
    
    ggsave(prplot, path = "figures/thermal-dimension/", 
           filename = paste(l$range[1], "_dormancy-p-vs-r.png",sep = ""), 
           height = 5, width = 11, units = "in", device = "png")
    
    c_potential_plot <- l %>%
      filter(type %in% c("potential", "potential_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal low temperature") + 
      geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    c_realized_plot <- l %>%
      filter(type %in% c("realized", "realized_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal low temperature") + 
      geom_vline(xintercept = ctmin, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    h_potential_plot <- h %>%
      filter(type %in% c("potential", "potential_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal high temperature", 
           title = paste("(", unique(h$dormancy), " season dormancy)", sep = "")) + 
      geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    h_realized_plot <- h %>%
      filter(type %in% c("realized", "realized_dormancy")) %>%
      ggplot(., aes(x = temps, fill = type)) +
      labs(fill = "Range:", y = "Frequency", x = "Seasonal high temperature", 
           title = unique(h$range)) + 
      geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) + geom_histogram(bins = 100, position ="identity", alpha = .6) 
    
    dndplot <- grid.arrange(h_realized_plot, h_potential_plot,
                            c_realized_plot, c_potential_plot, 
                           ncol = 2)
    
    ggsave(dndplot, path = "figures/thermal-dimension/", 
           filename = paste(l$range[1], "_dormancy-d-vs-nd.png",sep = ""), 
           height = 5, width = 11, units = "in", device = "png")
    
    num = num + 1
  } 
  else {
    num = num + 1
  }
}





## try plotting in 2D thermal niche space!!!

## first, reorganize temperature data
## we need to know which high and low temperatures belong to each raster cell

rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")

## loop through ranges and extract temperatures under each 
i = 1
realized_temps_2D <- c()
while (i < nlayers(rasterized_rrs) + 1) {
  range <- rasterized_rrs[[i]]
  
  ## use realized range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps, range)
  lows <- mask(low_temps, range)
  vals <- data.frame(range = names(rasterized_rrs)[i], type = "realized", high_temp = values(highs), 
                     low_temp = values(lows)) %>%
    filter(!is.na(high_temp) | !is.na(low_temp)) 
  
  if (nrow(vals) == 0) { ## if no temperatures 
    vals <- data.frame(range = names(rasterized_rrs)[i], type = "realized", high_temp = NA, 
                       low_temp = NA) 
  }
  
  
  if(i == 1) {
    realized_temps_2D <- vals
  }
  else {
    realized_temps_2D <- rbind(realized_temps_2D, vals)
  }
  i = i + 1
}

i = 1
potential_temps_2D <- c()
while (i < nlayers(potential_ranges) + 1) {
  range <- potential_ranges[[i]]
  
  ## use realized range raster as a mask to extract temperatures underneath
  highs <- mask(high_temps, range)
  lows <- mask(low_temps, range)
  vals <- data.frame(range = names(potential_ranges)[i], type = "potential", high_temp = values(highs), 
                     low_temp = values(lows)) %>%
    filter(!is.na(high_temp) | !is.na(low_temp)) 
  
  if (nrow(vals) == 0) { ## if no temperatures 
    vals <- data.frame(range = names(potential_ranges)[i], type = "potential", high_temp = NA, 
                       low_temp = NA) 
  }
  
  
  if(i == 1) {
    potential_temps_2D <- vals
  }
  else {
    potential_temps_2D <- rbind(potential_temps_2D, vals)
  }
  i = i + 1
}



all_temps_2D <- rbind(realized_temps_2D, potential_temps_2D)

all_temps_2D %>%
  filter(range == "Chamaeleo.dilepis_IUCN") %>%
  ggplot(., aes(x = high_temp, y = low_temp)) +
  geom_bin2d(aes(color = type))

lims <- thermal_limits[which(thermal_limits$genus_species == "Chamaeleo.dilepis"),]
ctmax <- lims$thermal_limit[which(lims$type == "max")]
ctmin <- lims$thermal_limit[which(lims$type == "min")]


all_temps_2D %>%
  filter(range == "Chamaeleo.dilepis_IUCN") %>%
  ggplot(., aes(x = high_temp, y = low_temp, colour = type)) +
  stat_density2d(geom = "density2d", aes(color = type, alpha=..level..),
                 size=2,
                 contour=TRUE) +
  labs(x = "Seasonal high temperature", y = "Seasonal low temperature", colour = "Range", alpha = "Density", title = "Chamaeleo.dilepis_IUCN") + geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) +
  geom_hline(yintercept = ctmin, linetype="dotted", size = 0.5) 


## try for all and save to a folder:

## get thermal limits:
lims <- thermal_limits[which(thermal_limits$genus_species == str_split_fixed(last(names(potential_ranges)), "_", n = 2)[1,1]),]
ctmax <- lims$thermal_limit[which(lims$type == "max")]
ctmin <- lims$thermal_limit[which(lims$type == "min")]








num = 1 
while (num < length(unique(potential_temps$range)) + 1) {
  p <- all_temps_2D %>%
    filter(range == unique(potential_temps$range)[num]) %>%
    filter(type == "potential")
  
  r <- all_temps_2D %>%
    filter(range == unique(potential_temps$range)[num]) %>%
    filter(type == "realized")
  
  ## get thermal limits:
  lims <- thermal_limits[which(thermal_limits$genus_species == str_split_fixed(p$range[1], "_", n = 2)[1,1]),]
  ctmax <- lims$thermal_limit[which(lims$type == "max")]
  ctmin <- lims$thermal_limit[which(lims$type == "min")]
  
  pplot <- p %>%
    ggplot(., aes(x = high_temp, y = low_temp)) +
    geom_bin2d(binwidth = 0.25) +
    scale_fill_gradient(low="lightblue", high="darkblue") +
    labs(fill = "Potential range frequency") 
    
  rplot <- r %>%
    ggplot(., aes(x = high_temp, y = low_temp)) +
    geom_bin2d(binwidth = 0.25) +
    scale_fill_gradient(low= "pink", high="purple") +
    labs(fill = "Realized range frequency")
  
  
    
  prplot <- ggplot(data = layer_data(pplot),
                                 aes(x = x, y = y, fill = fill)) +
    geom_tile() +
    geom_tile(data = layer_data(rplot)) +
    scale_fill_identity() +
    theme_bw() +
    geom_vline(xintercept = ctmax, linetype="dotted", size = 0.5) +
    geom_hline(yintercept = ctmin, linetype="dotted", size = 0.5) + 
    labs(x = "Seasonal high temperature (C)", y = "Seasonal low temperature (C)", title = p$range[1])
  
  library(cowplot)
  prplot <- plot_grid(prplot,
            plot_grid(get_legend(pplot),
                      get_legend(rplot),
                      ncol = 1),
            nrow = 1,
            rel_widths = c(1, 0.7)) 
  
  ggsave(prplot, path = "figures/thermal-dimension_2D/", 
         filename = paste(p$range[1], ".png",sep = "_"), 
         height = 6, width = 11, units = "in", device = "png")
  
  num = num + 1
}
