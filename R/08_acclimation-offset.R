## estimating 'acclimatized' thermal limits for each species 
library(tidyverse)
library(broom)

#########################################################
##       Collating ARRs from existing databases        ##
#########################################################
## check how many species we have are in intratherm:
intratherm <- read.csv("data-raw/ARR-data/intratherm-with-elev.csv") %>%
  mutate(genus_species = str_replace_all(genus_species, ' ', "_"))
taxize <- read.csv("data-processed/globtherm_taxize-key.csv") %>%
  mutate(binomial = str_replace_all(binomial, ' ', "_")) %>%
  mutate(acc_name = str_replace_all(acc_name, ' ', "_"))
traits <- read.csv("data-processed/rangetherm-traits_all-spp.csv")

intra_sp <- unique(intratherm$genus_species) 

taxize <- filter(taxize, binomial %in% traits$genus_species)
  
## check if species in our thermal data are in intratherm:
tp_correct <- unique(taxize$acc_name) 
tp_incorrect <- unique(taxize$binomial[which(!taxize$binomial %in% taxize$acc_name)])
overlap_intra <- tp_correct[tp_correct %in% intra_sp] ## 36 are in intratherm under their most correct name
overlap_incorrect <- tp_incorrect[tp_incorrect %in% intra_sp] ## 0 are in intratherm under their incorrect name

## change name in intratherm to be same as in rangetherm:
tax <- taxize %>%
  filter(acc_name %in% overlap_intra) 

intra_ol <- left_join(intratherm, tax, by = c("genus_species" = "acc_name")) %>%
  mutate(rangetherm_name =
           ifelse(is.na(binomial), as.character(genus_species), as.character(binomial))) %>%
  filter(rangetherm_name %in% tax$binomial)

## if data is from comte/rohr, remove from intratherm overlap:
intra_ol <- filter(intra_ol, !original_compilation %in% c("Rohr", "Comte"))

## how many of our species are in other data sets that were filtered down to create intratherm?
comte <- read.csv("data-raw/ARR-data/comte-all.csv") %>%
  mutate(Species = str_replace_all(Species, ' ', "_"))
rohr <- read.csv("data-raw/ARR-data/amphib-rohr.csv") %>%
  mutate(genus_species1 = paste(Genus1, Species1, sep = "_")) %>%
  mutate(genus_species2 = paste(Genus2, Species2, sep = "_")) %>%
  mutate(genus_species3 = paste(Genus3, Species3, sep = "_")) 

comte_sp <- unique(comte$Species)
rohr_sp1 <- unique(rohr$genus_species1)
rohr_sp2 <- unique(rohr$genus_species2)
rohr_sp2 <- unique(rohr$genus_species3)
rohr_sp2 <- rohr_sp2[!rohr_sp2 %in% rohr_sp1]
rohr_sp3 <- unique(rohr$genus_species3)
rohr_sp3 <- rohr_sp3[!rohr_sp3 %in% c(rohr_sp1, rohr_sp2)]

overlap_comte <- tp_correct[tp_correct %in% comte_sp] # 20 species are under acc name
overlap_incorrect_c <- tp_incorrect[tp_incorrect %in% comte_sp] # 0 under incorrect name
overlap_rohr_sp1 <- tp_correct[tp_correct %in% rohr_sp1] # 33 species are under acc name
overlap_incorrect_r1 <- tp_incorrect[tp_incorrect %in% rohr_sp1] # 9 species are under rangetherm name
overlap_rohr_sp2 <- tp_correct[tp_correct %in% rohr_sp2] # 4 species are under acc name
overlap_incorrect_r2 <- tp_incorrect[tp_incorrect %in% rohr_sp2] # 0 under incorrect name

## make a rangetherm names column and subset data 
comte_ol <- filter(comte, Species %in% overlap_comte) %>%
  mutate(rangetherm_name = Species)

rohr_ol <- filter(rohr, (genus_species1 %in% overlap_rohr_sp1) | 
                    (genus_species2 %in% overlap_rohr_sp2)) %>%
  mutate(rangetherm_name = ifelse(genus_species1 %in% overlap_rohr_sp1, 
                                  as.character(genus_species1), 
                                  as.character(genus_species2)))
  

## select only necessary columns: genus_species, acclimation temperature, thermal limit estimate
intra_ol <- intra_ol %>%
  select(rangetherm_name, genus_species, acclim_temp, ref, parameter_value, 
         parameter_tmax_or_tmin, metric_type, life_stage.x) %>%
  dplyr::rename("lifestage" = life_stage.x) %>%
  dplyr::rename("original_genus_species" = genus_species, "genus_species" = rangetherm_name)

comte_ol <- comte_ol %>%
  select(rangetherm_name, Species, Temperature.of.acclimation...C., Thermal.limit...C., Source,
         Methodology, Life.stage) %>%
  mutate(parameter_tmax_or_tmin = "tmax") %>%
  mutate(metric_type = ifelse(Methodology == "dynamic", "critical", "lethal")) %>%  # dynamic = CTmax, static = LT
  dplyr::rename("genus_species" = rangetherm_name, "original_genus_species" = Species, 
         "acclim_temp" = Temperature.of.acclimation...C., "parameter_value" = Thermal.limit...C., 
         "ref" = Source, "lifestage" = Life.stage) %>%
  select(-Methodology)

rohr_ol <- rohr_ol %>%
  select(rangetherm_name, genus_species1, Raw.CTM1, Raw.CTM2, Stage, Acclim.temp, Citation) %>%
  mutate(parameter_tmax_or_tmin = "tmax") %>%
  mutate(metric_type = "critical") %>%
  gather(key = "CT", value = "parameter_value", c(Raw.CTM1, Raw.CTM2)) %>%
  dplyr::rename("genus_species" = rangetherm_name, "original_genus_species" = genus_species1, 
                "acclim_temp" = Acclim.temp, 
                "ref" = Citation, "lifestage" = Stage) %>%
  select(-CT)

## combine datasets:
arr_data <- rbind(intra_ol, comte_ol, rohr_ol) 

## filter to species with known accliamtion temperatures, parameter values and at least two estimates
arr_data <- filter(arr_data, !is.na(acclim_temp), 
                   !is.na(parameter_value), 
                   genus_species %in% genus_species[which(duplicated(genus_species))])

tally <- arr_data %>%
  group_by(genus_species, parameter_tmax_or_tmin) %>%
  tally(length(unique(acclim_temp))) %>%
  filter(n != 1) %>%
  mutate(temp = paste(genus_species, parameter_tmax_or_tmin))

arr_data <- arr_data  %>%
  mutate(temp = paste(genus_species, parameter_tmax_or_tmin)) %>%
  filter(temp %in% tally$temp) %>%
  select(-temp)


##### how many species with only lower thermal limit in rangetherm now have a published upper limit?
## or vice versa?
tlims <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges_taxized.csv")

tmax <- filter(arr_data, parameter_tmax_or_tmin == "tmax")
min <- filter(tlims, type == "min", !duplicated(genus_species))

min$genus_species[which(min$genus_species %in% tmax$genus_species)] ## Acris crepitans

tmin <- filter(arr_data, parameter_tmax_or_tmin == "tmin")
max <- filter(tlims, type == "max", !duplicated(genus_species))

max$genus_species[which(max$genus_species %in% tmin$genus_species)] ## Anolis cybotes, Takydromus hsuehshanensis

##### how many species are in these other data bases that aren't in our globtherm?
gtherm <- read.csv("data-processed/globtherm_taxized.csv")

new_from_intra <- filter(intratherm, (!genus_species %in% gtherm$genus_species) & 
                           (!genus_species %in% gtherm$acc_name), 
                         !original_compilation %in% c("Rohr", "Comte"))
length(unique(new_from_intra$genus_species)) #15 species

new_from_rohr <- filter(rohr, (!(genus_species1 %in% gtherm$genus_species) & 
                          !genus_species1 %in% gtherm$acc_name), 
                        (!(genus_species2 %in% gtherm$genus_species) & 
                           !genus_species2 %in% gtherm$acc_name), 
                        (!(genus_species3 %in% gtherm$genus_species) & 
                           !genus_species3 %in% gtherm$acc_name))
length(unique(new_from_rohr$genus_species1)) #166 species

new_from_comte <- filter(comte, (!Species %in% gtherm$genus_species) & 
                                 (!Species %in% gtherm$acc_name))

length(unique(new_from_comte$Species)) # 341 species




#########################################################
##                Explore the ARR data                 ##
#########################################################
length(unique(arr_data$genus_species)) ## 35 species 

## plot ARRs
arr_data %>%
  filter(parameter_tmax_or_tmin == 'tmax') %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")
## one has a negative slope - probably need to remove it 

arr_data %>%
  filter(parameter_tmax_or_tmin == 'tmin') %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none")

## plan: explore how potential thermal range changes if we only use these species and only change thermal max

##### thinking it through:

## plot: where do the thermal limits we have fall along these ARRs?
tlims <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges_taxized.csv")

tlims_arr <- tlims %>%
  filter(genus_species %in% arr_data$genus_species) 

arr_upper <- arr_data %>%
  filter(parameter_tmax_or_tmin == 'tmax') 
arr_lower <- arr_data %>%
  filter(parameter_tmax_or_tmin == 'tmin')

upper <- filter(tlims_arr, type == "max", genus_species %in% arr_upper$genus_species)
lower <- filter(tlims_arr, type == "min", genus_species %in% arr_lower$genus_species)

library(RColorBrewer)
cols = length(unique(upper$genus_species))
pal_upper <- as.vector(colorRampPalette((brewer.pal(9, "Reds")))(cols))

arr_upper %>%
  filter(genus_species %in% upper$genus_species) %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none") + 
  labs(x = "Acclimation temperature (C)", y = "Upper thermal limit (C)") +
  geom_point(data = upper, aes(x = acclimation_temperature, y = thermal_limit)) +
  scale_color_manual(values = pal_upper)

cols = length(unique(lower$genus_species))
pal_lower <- as.vector(colorRampPalette((brewer.pal(9, "Blues")))(cols))

arr_lower %>%
  filter(genus_species %in% lower$genus_species) %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value, group = genus_species, col = genus_species)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none") +
  geom_point(data = lower, aes(x = acclimation_temperature, y = thermal_limit)) +
  labs(x = "Acclimation temperature (C)", y = "Lower thermal limit (C)") +
  scale_color_manual(values = pal_lower)

## calculate: distance between max and min thermal limit within sampled acclimation temps
# ie. how much could each thermal limit change by in each location due to acclimation? would it be worth it?

split <- split(arr_upper, f = arr_upper$genus_species) 

for (i in 1:length(split)) {
  data = split[[i]]
  
  ## get max and min simpled temps:
  new <- data.frame(acclim_temp = c(min(data$acclim_temp), max(data$acclim_temp)))
  
  lm <- lm(data = data, parameter_value ~ acclim_temp)
  
  new$pred_parameter_val <- predict(lm, new)
  
  if (i == 1) {
    minmax <- data.frame(genus_species = rep(unique(data$genus_species), 2), cbind(new))
  }
  else {
    minmax_new <- data.frame(genus_species = rep(unique(data$genus_species), 2), cbind(new))
    minmax = rbind(minmax, minmax_new)
  }
}

dif <- minmax %>%
  group_by(genus_species) %>%
  mutate(diff = pred_parameter_val - lag(pred_parameter_val, 
                                         default = first(pred_parameter_val)))

dif <- dif[seq(from = 2, to = nrow(dif), by = 2), ]

ggplot(dif, aes(x = diff)) + geom_histogram()



###############################################################
##                    fit global linear model                ##
###############################################################
## collate all ARR data together, making sure none are duplicated:
# gunderson 
# 
# morley <- read_csv("data-raw/ARR-data/morley.csv")[8:380,] 
# colnames(morley) <- morley[1,]
# morley <- morley[-1,]
# 
# rohr <- read.csv("data-raw/ARR-data/amphib-rohr.csv") %>%
#   mutate(genus_species1 = paste(Genus1, Species1, sep = "_")) %>%
#   mutate(genus_species2 = paste(Genus2, Species2, sep = "_")) %>%
#   mutate(genus_species3 = paste(Genus3, Species3, sep = "_")) 


## low hanging fruit:
intratherm <- read.csv("data-raw/ARR-data/intratherm-with-elev.csv") %>%
  mutate(genus_species = str_replace_all(genus_species, ' ', "_")) %>%
  select(genus_species, acclim_temp, ref, parameter_value, 
         parameter_tmax_or_tmin, metric_type, life_stage.x) %>%
  dplyr::rename("lifestage" = life_stage.x)

comte <- read.csv("data-raw/ARR-data/comte-all.csv") %>%
  mutate(Species = str_replace_all(Species, ' ', "_")) %>%
  select(Species, Temperature.of.acclimation...C., Thermal.limit...C., Source,
         Methodology, Life.stage) %>%
  mutate(parameter_tmax_or_tmin = "tmax") %>%
  mutate(metric_type = ifelse(Methodology == "dynamic", "critical", "lethal")) %>%  # dynamic = CTmax, static = LT
  dplyr::rename("genus_species" = Species, 
                "acclim_temp" = Temperature.of.acclimation...C., "parameter_value" = Thermal.limit...C., 
                "ref" = Source, "lifestage" = Life.stage) %>%
  select(-Methodology)

## combine datasets:
arr <- rbind(intratherm, comte) 

## filter to species with known accliamtion temperatures, parameter values and at least two estimates
arr <- filter(arr, !is.na(acclim_temp), 
                   !is.na(parameter_value), 
                   genus_species %in% genus_species[which(duplicated(genus_species))])

tally <- arr %>%
  group_by(genus_species, parameter_tmax_or_tmin) %>%
  tally(length(unique(acclim_temp))) %>%
  filter(n != 1) %>%
  mutate(temp = paste(genus_species, parameter_tmax_or_tmin))

arr <- arr  %>%
  mutate(temp = paste(genus_species, parameter_tmax_or_tmin)) %>%
  filter(temp %in% tally$temp) %>%
  select(-temp)

upper <- arr %>%
  filter(parameter_tmax_or_tmin == "tmax") 
lower <- arr %>%
  filter(parameter_tmax_or_tmin == "tmin") %>%
  filter(metric_type == "critical")

library(nlme)
upper_model <- upper %>%
  do(tidy(lme(data = ., parameter_value ~ acclim_temp, random = ~1|genus_species))) 

upper_slope <- upper_model %>%
  filter(term=='acclim_temp') %>%
  .$estimate %>%
  unique(.)
lower_slope <- tidy(lme(data = lower, parameter_value ~ acclim_temp, 
                             random = ~1|genus_species)) %>%
  filter(term == "acclim_temp") %>%
  .$estimate %>%
  unique(.)

# upper_quantile <- quantile(upper$acclim_temp, na.rm=TRUE, c(0.05, 0.95))
# lower_quantile <- quantile(lower$acclim_temp, na.rm=TRUE, c(0.05, 0.95))

## predict values so they can be plotted:
library(MuMIn)
upper_new_data <- data.frame(expand_grid(acclim_temp = seq(min(upper$acclim_temp),
                                                    max(upper$acclim_temp),
                                                    length.out = 100)))
lower_new_data <- data.frame(expand_grid(acclim_temp = seq(min(lower$acclim_temp),
                                                           max(lower$acclim_temp),
                                                           length.out = 100)))

pred_upper <- upper %>%
  lme(data = ., parameter_value ~ acclim_temp, random = ~1|genus_species) %>%
  predict(., upper_new_data, level = 0)
pred_lower <- lower %>%
  lme(data = ., parameter_value ~ acclim_temp, random = ~1|genus_species) %>%
  predict(., lower_new_data, level = 0)

fitted_pred_upper <- upper_new_data %>%
  mutate(thermal_limit = pred_upper, type = "max") 
fitted_pred_lower <- lower_new_data %>%
  mutate(thermal_limit = pred_lower, type = "min")


upper %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point(col = "darkred") +
  geom_smooth(data = fitted_pred_upper, aes(x = acclim_temp, y = thermal_limit), inherit.aes = F,
              col = "black") + 
  #geom_vline(xintercept = upper_quantile) +
  labs(x = "Acclimation temperature (C)", y = "Upper thermal limit (C)", col = "Metric:")

lower %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point(col = "darkblue") +
  geom_smooth(data = fitted_pred_lower, aes(x = acclim_temp, y = thermal_limit), inherit.aes = F,
              col = "black") +
  # geom_vline(xintercept = lower_quantile) +
  labs(x = "Acclimation temperature (C)", y = "Lower thermal limit (C)", col = "Metric:")



## solve for species-specific intercepts:
acc_data <- tlims %>%
  mutate(ARR_equ_slope = ifelse(!is.na(acclimation_temperature) & type == "max", 
                                upper_slope, 
                                ifelse(!is.na(acclimation_temperature) & type == "min",
                                       lower_slope,
                                       NA))) %>%
  mutate(ARR_equ_int = thermal_limit - ARR_equ_slope*acclimation_temperature) %>%
  mutate(max_upper_acc_temp = max(upper$acclim_temp), min_upper_acc_temp = min(upper$acclim_temp),
         max_lower_acc_temp = max(lower$acclim_temp), min_lower_acc_temp = min(lower$acclim_temp)) %>%
  group_by(genus_species) %>%
  filter(any(!is.na(ARR_equ_slope)))

# acc_data$ARR_equ_lowlim <- ifelse(acc_data$type == "max", upper_quantile[1], 
#                                lower_quantile[1])
# acc_data$ARR_equ_highlim <- ifelse(acc_data$type == "max", upper_quantile[2], 
#                                 lower_quantile[2])


## how many?
acc_data %>%
  group_by(genus_species)  %>%
  filter(n() == 2) 

## plot global arr and all of our data:
acc_data %>%
  filter(type ==  "max") %>%
  ggplot(., aes(x = acclimation_temperature, y = thermal_limit)) + geom_point(col = "darkred") +
  geom_smooth(data = fitted_pred_upper, aes(x = acclim_temp, y = thermal_limit), inherit.aes = F,
              col = "black") + 
  labs(x = "Acclimation temperature (C)", y = "Upper thermal limit (C)")

acc_data %>%
  filter(type ==  "min") %>%
  ggplot(., aes(x = acclimation_temperature, y = thermal_limit)) + geom_point(col = "darkblue") +
  geom_smooth(data = fitted_pred_lower, aes(x = acclim_temp, y = thermal_limit), inherit.aes = F,
              col = "black") + 
  labs(x = "Acclimation temperature (C)", y = "Lower thermal limit (C)", col = "Metric:")


write.csv(acc_data, "data-processed/acclimation-data.csv", row.names = FALSE)

###############################################################
##            fit global non-linear mixed effect model       ##
###############################################################
## parameterize nonlinear model:
input <- filter(upper, !is.na(acclim_temp) & !is.na(parameter_value)) 

initial <- getInitial(parameter_value ~ SSasympOff(acclim_temp, A, lrc, c0), 
                      data = input)

acc_hot <- nls(parameter_value ~ SSasympOff(acclim_temp, A, lrc, c0), 
               data = input,
               control =  nls.control(maxiter = 50000, minFactor = 0),
               trace = TRUE)

coef(acc_hot)

model <- function(tlim, acclim_temp){
  a  <- tlim[1]
  b <- tlim[2]
  c  <- tlim[3]
  f  <-  a + (b-a)*exp(1)^(-c*acclim_temp)
  return(f)
}

hot_data <- upper %>%
  group_by(acclim_temp) %>%
  dplyr::summarize(parameter_value = mean(parameter_value, na.rm = TRUE))

hot_data %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = upper_quantile)

acc_hot <- nls(parameter_value ~ model(tlim, acclim_temp), 
               start = list(tlim=c(a=53.51939380, b=20.75928673, c=0.02690044)),
               data = hot_data,
               control =  nls.control(maxiter = 10000, minFactor = 0))
coef(acc_hot)
# a = 53.51955
# b = 20.75928673 
# c = 0.02690044 



new.df <- data.frame(acclim_temp=seq(0,45,by=0.1))
new.df$pred1 <- predict(acc_hot, newdata=new.df)

hot_data %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = upper_quantile) +
  geom_line(data=new.df, aes(x=acclim_temp,y=pred1), colour="#339900", size=1)

upper %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = upper_quantile) +
  geom_line(data=new.df, aes(x=acclim_temp,y=pred1), colour="#339900", size=1)



cold_data <- lower %>% 
  filter(metric_type != "lethal") %>%
  group_by(acclim_temp) %>%
  dplyr::summarize(parameter_value = mean(parameter_value, na.rm = TRUE))

lower %>%
  filter(metric_type != "lethal") %>%
  ggplot(., aes(x = acclim_temp, y = parameter_value)) + geom_point() +
  geom_smooth(method = "lm") + geom_vline(xintercept = upper_quantile)

acc_cold <- nls(parameter_value ~ model(tlim, acclim_temp), 
                start = list(tlim=c(a= 21.58533181, b=-3.90798748 , c=0.01942988 )),
                data = lower, 
                control =  nls.control(maxiter = 10000, minFactor = 0))
coef(acc_cold)
# a = 21.58533181
# b = -3.90798748 
# c = 0.01942988  





## garbage:


te_data <- readRDS("data-processed/Te_allspp.rds")

split <- split(arr_upper, f = arr_upper$genus_species) 

acc_te_data <- list()
for (i in 1:length(split)) {
  data = split[[i]]
  
  ## get sp:
  sp <- data$genus_species[1]
  
  ## get og tlims:
  og_lims <- tlims_arr[which(tlims_arr$genus_species == sp),]
  ## if does not have min and max, skip for now:
  if (nrow(og_lims) != 2 & og_lims$type == "min") {
    acc_te_data[[i]] <- NA
  }
  else if(og_lims$realm != "Terrestrial") {
    acc_te_data[[i]] <- NA
  }
  else {
    ## get max and min sampled temps:
    min <- min(data$acclim_temp)
    max <- max(data$acclim_temp)
    min_limit <- minmax$pred_parameter_val[which(minmax$acclim_temp == min & 
                                                   minmax$genus_species == sp)]
    max_limit <- minmax$pred_parameter_val[which(minmax$acclim_temp == max &
                                                   minmax$genus_species == sp)]
    
    ## get acclimation temperatures
    te <- te_data[[which(names(te_data) == sp)]]
    acc_temps <- te_data[[which(names(te_data) == sp)]] %>%
      select(hot_acc_temp) %>%
      dplyr::rename("acclim_temp" = hot_acc_temp)
    
    lm <- lm(data = data, parameter_value ~ acclim_temp)
    
    acc_temps$pred_upper_lim <- predict(lm, acc_temps)
    
    ## if acclimation temp is beyond those sampled, assign max/min limit:
    te$hot_acc_limit <- ifelse(acc_temps$acclim_temp < min, min_limit, 
                               ifelse(acc_temps$acclim_temp > max, max_limit, 
                                      acc_temps$pred_upper_lim))
    
    ## compare to global acclimation response ratio:
    tlims_sp <- filter(tlims, genus_species == sp & type == "max")
    if (nrow(tlims_sp) != 0) {
      te$hot_acc_limit_global <- tlims_sp$ARR_equ_slope*te$hot_acc_temp +
        tlims_sp$ARR_equ_int
      min_limit <- tlims_sp$ARR_equ_slope*tlims_sp$ARR_equ_lowlim +
        tlims_sp$ARR_equ_int
      max_limit <- tlims_sp$ARR_equ_slope*tlims_sp$ARR_equ_highlim +
        tlims_sp$ARR_equ_int
      te$hot_acc_limit_global <- ifelse(te$hot_acc_temp < tlims_sp$ARR_equ_lowlim, min_limit,
                                               ifelse(te$hot_acc_temp > tlims_sp$ARR_equ_highlim, 
                                                      max_limit,
                                                      te$hot_acc_limit_global))
    }
    else {
      te$hot_acc_limit_global = NA
    }
    tlims_sp <- filter(tlims, genus_species == sp & type == "min")
    if (nrow(tlims_sp) != 0) {
      te$cold_acc_limit_global <- tlims_sp$ARR_equ_slope*te$cold_acc_temp +
        tlims_sp$ARR_equ_int
      min_limit <- tlims_sp$ARR_equ_slope*tlims_sp$ARR_equ_lowlim +
        tlims_sp$ARR_equ_int
      max_limit <- tlims_sp$ARR_equ_slope*tlims_sp$ARR_equ_highlim +
        tlims_sp$ARR_equ_int
      te$cold_acc_limit_global <- ifelse(te$cold_acc_temp < tlims_sp$ARR_equ_lowlim, min_limit,
                                                ifelse(te$cold_acc_temp > tlims_sp$ARR_equ_highlim, 
                                                       max_limit,
                                                       te$cold_acc_limit_global))
    } 
    else {
      te$cold_acc_limit_global = NA
    }
    
    acc_te_data[[i]] <- te
  }
}

acc_te_data <- acc_te_data[-which(is.na(acc_te_data))] 


test <- do.call(rbind.data.frame, acc_te_data)

ggplot(test, aes(x = hot_acc_limit, y = hot_acc_limit_global)) + geom_point() + geom_abline()





#### plan for acclimation later: ####
# 1. make acclimation temp grids for marine and intertidal spp. 
#   - calculate means in 30 days before max and min temps
# terr: 
#   terr_acc_cold <- stack("data-processed/terr_acc_cold.grd")
#   terr_acc_hot <- stack("data-processed/terr_acc_hot.grd")

# 
# 2. turn calculated acclimation temps into grids for terrestrial spp
# 
# 3. make function 
#   -  sees what realm spp is from
#   -  calculates a grid of 'acclimated' thermal limits 
#   -  filters by acclimated tolerance 


