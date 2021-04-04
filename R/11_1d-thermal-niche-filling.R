### calculating 1D potential thermal niche filling for term paper
## with sensitivity to elevation and dormancy
library(tidyverse)

###################################################
##                 read in data                  ##
###################################################
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges.csv") %>%
  mutate(genus_species = paste(Genus, Species, sep = " "))


## calculate thermal niche filling for all sensitivity sets separately 
###################################################
##         elevation + dormancy adjustment       ##
###################################################
## split potnetial and realized
niche_elev_x_dormancy <- read.csv("./data-processed/thermal-niche/niche_elev_x_dormancy.csv")

r_niches <- niche_elev_x_dormancy %>%
  filter(type == 'realized_elev_x_dormancy')
p_niches <- niche_elev_x_dormancy %>%
  filter(type == 'potential_elev_x_dormancy')

## collapse all temperatures in niches to only temperatures at the limits of the niches
niche_lims_elev_x_dormancy <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = 'elev_x_dormancy')

## calculate niche filling metrics:
filling_elev_x_dormancy <- calculate_filling(niche_lims_elev_x_dormancy)

###################################################
##             dormancy adjustment               ##
###################################################
niche_dormancy <- read.csv("./data-processed/thermal-niche/niche_dormancy.csv")

## split potnetial and realized
r_niches <- niche_dormancy %>%
  filter(type == 'realized_dormancy')
p_niches <- niche_dormancy %>%
  filter(type == 'potential_dormancy')

niche_lims_dormancy <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = 'dormancy')

filling_dormancy <- calculate_filling(niche_lims_dormancy)

###################################################
##                 elev adjustment               ##
###################################################
niche_elev <- read.csv("./data-processed/thermal-niche/niche_elev.csv")

## split potnetial and realized
r_niches <- niche_elev %>%
  filter(type == 'realized_elevation')
p_niches <- niche_elev %>%
  filter(type == 'potential_elevation')

niche_lims_elev <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = 'elevation')

filling_elev <- calculate_filling(niche_lims_elev)

###################################################
##                 no adjustments                ##
###################################################
niche_reg <- read.csv("./data-processed/thermal-niche/niche_reg.csv")

## split potnetial and realized
r_niches <- niche_reg %>%
  filter(type == 'realized_reg')
p_niches <- niche_reg %>%
  filter(type == 'potential_reg')

niche_lims_reg <- extract_niche_limits(r_niches, p_niches, thermal_limits, type = "reg")

filling_reg <- calculate_filling(niche_lims_reg)



###################################################
##         exporting data for analysis           ##
###################################################
filling_all <- rbind(filling_reg, filling_elev, filling_dormancy, filling_elev_x_dormancy)

write.csv(filling_all, "data-processed/thermal-niche/thermal-niche-filling-metrics.csv", row.names = FALSE)


###################################################
##          metric calculating functions         ##
###################################################

## function that takes all temperatutes in a species realized and potential thermal niche, ctmax and ctmin,
## and extracts niche limits
## collapses potential_niche to a single minimum and maximum available temperature 
## collapses realized_niche to a single minimum and maximum occupied temperature 
extract_niche_limits <- function(r_niches, p_niches, thermal_limits, type) {

  ## loop through each range:
  range = 1
  while (range < length(unique(r_niches$range)) + 1) {
    
    range_id <- unique(r_niches$range)[range]
   
    ## extract the p and r niche temps for range:
    rniche <- r_niches[which(r_niches$range == range_id), ]
    pniche <- p_niches[which(p_niches$range == range_id), ]
    
    ## get species thermal limits
    type = str_split_fixed(unique(rniche$type), '_', n=2)[1,2]
    if(type == 'elevation') {
      split <- str_split_fixed(range_id, '_', n = 4) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
    }
    else if (type == 'elev_x_dormancy') {
      split <- str_split_fixed(range_id, '_', n = 5) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
    }
    else if(type == 'dormancy') {
      split <- str_split_fixed(range_id, "_", n = 4)
      species <- paste(split[1,1], split[1,2], sep = " ")
    }
    else {
      split <- str_split_fixed(range_id, '_', n = 3) 
      species <- paste(split[1,1], split[1,2], sep = ' ')
    }
    lims <- thermal_limits[which(thermal_limits$genus_species == species),]
    ctmax <- filter(lims, type == 'max') %>%
      .$thermal_limit
    ctmin <- filter(lims, type == 'min') %>%
      .$thermal_limit
    
    ## extract niche limits:
    r_upper <- max(rniche$temps, na.rm = TRUE)
    r_lower <- min(rniche$temps, na.rm = TRUE)
    p_upper <- max(pniche$temps, na.rm = TRUE)
    p_lower <- min(pniche$temps, na.rm = TRUE)
    
    ## add to df:
    if (range == 1) {
      niche_lims <- data.frame(range = range_id, 
                               type = type,
                               dormancy = unique(rniche$dormancy),
                               realm = unique(rniche$realm),
                               r_niche_upper = r_upper,
                               p_niche_upper = p_upper,
                               ctmax = ctmax,
                               r_niche_lower = r_lower,
                               p_niche_lower = p_lower, 
                               ctmin = ctmin)
    }
    else {
      niche_lims <- rbind(niche_lims, data.frame(range = range_id, 
                               type = type,
                               dormancy = unique(rniche$dormancy),
                               realm = unique(rniche$realm),
                               r_niche_upper = r_upper,
                               p_niche_upper = p_upper,
                               ctmax = ctmax,
                               r_niche_lower = r_lower,
                               p_niche_lower = p_lower, 
                               ctmin = ctmin))
    }
    
    range = range + 1
    
  }
  
  return(niche_lims)
}


## function that takes a species realized + potential thermal niche limits and calculates:
##    1. warm niche underfilling: °C of available thermal niche above the maximum realized temperature but 
##       below CTmax, multiplied by -1
##    2. warm niche overfilling: °C of occupied thermal niche above the maximum available temperature 
##    3. cold niche underfilling: °C of available thermal niche below the minimum realized temperature but 
##       above CTmin, multiplied by -1 
##    4. cold niche overfilling: °C of occupied thermal niche below the minimum available temperature 
##    5. warm fundamental niche overfilling: °C of occupied thermal niche above the ctmax
##    6. cold fundamental niche overfilling: °C of occupied thermal niche below the ctmin
calculate_filling <- function(niche_lims) {
  ## warm under: (p niche upper - r niche upper) - if negative, assign 0 and if positive, multiply by -1
  ## warm over: (r niche upper - p niche upper) - if negative, assign 0 and if positive, leave
  ## warm fundamental under: (r niche upper - ctmax) - if positive, assign 0 and if negative, multiply by -1
  ## warm fundamental over: (r niche upper - ctmax) - if negative, assign 0 and if positive, leave
  ## cold under: (p niche lower - r niche lower) - if positive, assign 0 and if negative, leave
  ## cold over: (r niche lower - p niche lower) - if positive, assign 0 and if negative, multiply by -1
  ## cold fundamental over: (r niche lower - ctmin) - if positive, assign 0 and if positive, leave
  ## cold fundamental under: (r niche lower - ctmin) - if negative, assign 0 and if positive, multiply by -1
  niche_filling <- niche_lims %>%
    mutate(warm_under = ifelse((p_niche_upper - r_niche_upper) < 0, 0, -(p_niche_upper - r_niche_upper))) %>%
    mutate(warm_over = ifelse((r_niche_upper - p_niche_upper) < 0, 0, (r_niche_upper - p_niche_upper))) %>%
    mutate(f_warm_over = ifelse((r_niche_upper - ctmax) < 0, 0, (r_niche_upper - ctmax))) %>%
    mutate(f_warm_under = ifelse((r_niche_upper - ctmax) > 0, 0, -(r_niche_upper - ctmax))) %>%
    mutate(cold_under = ifelse((p_niche_lower - r_niche_lower) > 0, 0, (p_niche_lower - r_niche_lower))) %>%
    mutate(cold_over = ifelse((r_niche_lower - p_niche_lower) > 0, 0, -(r_niche_lower - p_niche_lower))) %>%
    mutate(f_cold_over = ifelse((r_niche_lower - ctmin) > 0, 0, (r_niche_upper - ctmin))) %>%
    mutate(f_cold_under = ifelse((r_niche_lower - ctmin) < 0, 0, -(r_niche_lower - ctmin))) 
  
  return(niche_filling)
}
