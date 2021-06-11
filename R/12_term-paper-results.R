## analysis and figures for term paper
library(tidyverse)
library(cowplot)
library(lme4)
library(nlme)

##########################################################
##              prepping for model fitting              ##
##########################################################
uofill <- read.csv("data-processed/thermal-niche/thermal-niche-filling-metrics.csv") 
traits <- read.csv("data-processed/ectotherm-traits_all-spp.csv")
thermal_limits <- read.csv("data-processed/thermal-limits_ectotherms-with-ranges_taxized.csv")%>%
  select(genus_species, type, metric) %>%
  unique(.)

## make one column for both cold and warm filling, with edge_type specifying whether cold or warm
## make species and source columns 
data <- uofill %>%
  mutate(genus_species = paste(str_split_fixed(.$range, '_', 3)[,1], 
                         str_split_fixed(.$range, '_', 3)[,2], sep = '_')) %>%
  mutate(source = str_split_fixed(.$range, '_', 3)[,3]) %>%
  select(range, genus_species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GBIF"), ordered = TRUE)) %>%
  arrange(genus_species, type, source) 


## give priority to IUCN ranges (but test sensitivity later)
data <- data %>%
  mutate(temp = paste(genus_species, type, sep = '')) %>%
  filter(!duplicated(temp)) %>%
  select(-temp) %>%
  rename("sensitivity_type" = type) 

data <- data %>%
  mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
  mutate(perfect_cold = ifelse(cold_over == 0 & cold_under == 0, 0, NA)) %>%
  gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                          cold_over, perfect_warm, perfect_cold)) %>%
  filter(filling_value != 0 | filling_value == 0 & (filling_type %in% c("perfect_cold","perfect_warm"))) %>%
  select(range, sensitivity_type,
         genus_species, source, realm, dormancy, filling_type, filling_value, lat_mp, range_area_km2) %>%
  mutate(edge_type = ifelse(str_detect(filling_type, "warm"), "warm", "cold")) %>%
  mutate(type = ifelse(str_detect(filling_type, "max"), "max", "min")) %>%
  arrange(genus_species, filling_type, sensitivity_type, source) %>%
  filter(!is.infinite(filling_value)) ## get rid of infinite filling values:


## merge with species traits and geographic predictors:
data <- left_join(data, traits, by = c("genus_species")) %>%
  left_join(., thermal_limits, by = c("genus_species", "type"))

## write out:
write.csv(data, "data-processed/thermal-niche/thermal-niche-filling-metrics_model-ready.csv", row.names = FALSE)

## get stats about taxa: 
data %>%
  filter(realm == "terrestrial") %>%
  filter(!duplicated(range)) %>%
  count(Phylum)


##########################################################
##                     model fitting                    ##
##########################################################
colnames(data)[29] <- "dispersal_distance_category"
colnames(data)[31] <- "dispersal_ability_category"

## split data by sensitivity type
sens_types <- group_split(data, sensitivity_type)

## for now, fit model to version that incorporates elevation permittivity and dormancy:
elev_x_dormancy <- as.data.frame(sens_types[[2]])

## try making dispersal distance continuous:
cont_disp <- elev_x_dormancy %>%
  mutate(dispersal_distance_continuous = ifelse(dispersal_distance_category == "0-1", 1, 
                                                ifelse(dispersal_distance_category == "1-10", 10,
                                                       ifelse(dispersal_distance_category == "10-100", 100,
                                                              ifelse(dispersal_distance_category == "100", 1000, 
                                                                     NA)))))

cont_disp$abs_lat_mp <- abs(cont_disp$lat_mp)
levels(cont_disp$realm) <-  c("intertidal", "marine", "terrestrial")

## fitting to cold and warm filling separately
cold <- cont_disp %>%
  filter(edge_type == "cold")

warm <- cont_disp %>%
  filter(edge_type == "warm")

model_fit_cold <- lme(filling_value ~ abs_lat_mp + 
                        realm + 
                        maximum_body_size_SVL_HBL_cm_ + 
                        dispersal_distance_continuous  + 
                        metric, 
                      
                      random = ~1|Class/Order/Family/Genus, 
                      
                      data = cold, na.action = na.exclude)

summary(model_fit_cold)

model_fit_warm <- lme(filling_value ~ abs_lat_mp + 
                        realm + 
                        maximum_body_size_SVL_HBL_cm_ + 
                        dispersal_distance_continuous  + 
                        metric, 
                      
                      random = ~1|Class/Order/Family/Genus, 
                      
                      data = warm, na.action = na.exclude)

summary(model_fit_warm)

model_fit_warm <- lme(filling_value ~ abs_lat_mp*realm +
                        maximum_body_size_SVL_HBL_cm_ + 
                        dispersal_distance_continuous  + 
                        metric, 
                      
                      random = ~1|Class/Order/Family/Genus, 
                      
                      data = warm, na.action = na.exclude)

summary(model_fit_warm)




## variable correlation values - colinearity, within models do the variables act the same way? vic 
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v 
}

vif.lme(model_fit_warm)
vif.lme(model_fit_cold)

## frequentist approach - coeffs and confidence intervals "certainty"
## check out area of geographic range 

# filling_value ~ lat_mp + realm + range_area_km2 + Trophic_position + 
#   maximum_body_size_SVL_HBL_cm_ + dispersal_distance_category*edge_type +
#   dispersal_ability_category*edge_type + hot_season_dormancy_*edge_type +
#   cold_season_dormancy_*edge_type + metric
# 
# random = ~1|Class/Order/Family/Genus


##########################################################
##                 plotting model results               ##
##########################################################
## making predictions using model results 
library(MuMIn)

new_data <- data.frame(expand_grid(abs_lat_mp = seq(min(cont_disp$abs_lat_mp),
                                                                    max(cont_disp$abs_lat_mp),
                                                                   length.out = 180),
                               realm = c('intertidal', 'marine', "terrestrial"),
                               maximum_body_size_SVL_HBL_cm_ = 16.23225,
                               dispersal_distance_continuous = 1,
                               metric = c("ct", 'LT50')))

pred_cold <- predict(model_fit_cold, new_data, level = 0, se.fit = T)
pred_warm <- predict(model_fit_warm, new_data, level = 0, se.fit = T)

fitted_pred_cold <- new_data %>%
  mutate(filling_value = pred_cold$fit) %>%
  mutate(filling_value_SE = pred_cold$se.fit)

fitted_pred_warm<- new_data %>%
  mutate(filling_value = pred_warm$fit) %>%
  mutate(filling_value_SE = pred_warm$se.fit)

##########################################################
##                   figure prepping                    ##
##########################################################
## plot of model predictions:
cont_disp_int <- cont_disp %>%
  filter(realm == "intertidal") 

realm_cold_int <- fitted_pred_cold %>%
  filter(realm == "intertidal") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_int$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_int$abs_lat_mp))
  
realm_warm_int <- fitted_pred_warm %>%
  filter(realm == "intertidal") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_int$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_int$abs_lat_mp))

int_pred_realm <- cont_disp_int %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "",
       x = "Realized range latitudinal midpoint (°N/S)")  +
  scale_x_continuous(breaks = c(seq(0, 90, 30)), labels = c("0", "30", "60", "90"), 
                     limits = c(0, 66)) +
  geom_line(data = realm_warm_int, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = realm_warm_int, inherit.aes = F, aes(x = abs_lat_mp,
                                                      ymin = (filling_value-1.96*filling_value_SE),
                                    ymax = filling_value+1.96*filling_value_SE), alpha=0.1, 
              fill = "#D55E00") +
  geom_line(data = realm_cold_int, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), col = "#0072B2") +
  geom_ribbon(data = realm_cold_int, inherit.aes = F, aes(x = abs_lat_mp,
                                                      ymin = (filling_value-1.96*filling_value_SE),
                                                      ymax = filling_value+1.96*filling_value_SE), 
                                                      alpha=0.1, fill = "#0072B2") +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) +
  scale_y_continuous(limits = c(-34, 45))

cont_disp_mar <- cont_disp %>%
  filter(realm == "marine")

realm_cold_mar <- fitted_pred_cold %>%
  filter(realm == "marine") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_mar$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_mar$abs_lat_mp))

realm_warm_mar <- fitted_pred_warm %>%
  filter(realm == "marine") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_mar$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_mar$abs_lat_mp))

mar_pred_realm <- cont_disp_mar %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "",
       x = "Realized range latitudinal midpoint (°N/S)")  +
  scale_x_continuous(breaks = c(seq(0, 90, 30)), labels = c("0", "30", "60", "90"),
                     limits = c(0, 66)) +
  geom_line(data = realm_warm_mar, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = realm_warm_mar, inherit.aes = F, aes(x = abs_lat_mp,
                                                          ymin = (filling_value-1.96*filling_value_SE),
                                                          ymax = filling_value+1.96*filling_value_SE),
                                                          fill = "#D55E00", alpha=0.1) +
  geom_line(data = realm_cold_mar, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), col = "#0072B2") +
  geom_ribbon(data = realm_cold_mar, inherit.aes = F, aes(x = abs_lat_mp,
                                                          ymin = (filling_value-1.96*filling_value_SE),
                                                          ymax = filling_value+1.96*filling_value_SE), 
                                                          alpha=0.1, fill = "#0072B2") +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) +
  scale_y_continuous(limits = c(-34, 45))

cont_disp_ter <- cont_disp %>%
  filter(realm == "terrestrial")

realm_cold_ter <- fitted_pred_cold %>%
  filter(realm == "terrestrial") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_ter$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_ter$abs_lat_mp))

realm_warm_ter <- fitted_pred_warm %>%
  filter(realm == "terrestrial") %>%
  filter(metric == "ct") %>%
  filter(abs_lat_mp > min(cont_disp_ter$abs_lat_mp) & 
           abs_lat_mp < max(cont_disp_ter$abs_lat_mp))

ter_pred_realm <- cont_disp_ter %>%
  ggplot(., aes(x = abs_lat_mp, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "Shortfall or excess of temperatures occupied (°C)",
       x = "Realized range latitudinal midpoint (°N/S)")  +
  scale_x_continuous(breaks = c(seq(0, 90, 30)), labels = c("0", "30", "60", "90"), 
                     limits = c(0, 66)) +
  geom_line(data = realm_warm_ter, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = realm_warm_ter, inherit.aes = F, aes(x = abs_lat_mp,
                                                          ymin = (filling_value-1.96*filling_value_SE),
                                                          ymax = filling_value+1.96*filling_value_SE),
              alpha=0.1, fill = "#D55E00") +
  geom_line(data = realm_cold_ter, inherit.aes = F, aes(x = abs_lat_mp, y = filling_value), col = "#0072B2") +
  geom_ribbon(data = realm_cold_ter, inherit.aes = F, aes(x = abs_lat_mp,
                                                          ymin = (filling_value-1.96*filling_value_SE),
                                                          ymax = filling_value+1.96*filling_value_SE), 
              alpha=0.1, fill = "#0072B2") +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) +
  scale_y_continuous(limits = c(-34, 45))

pred_realm <- ggdraw() +
  draw_plot(ter_pred_realm, 0, 0, 1/3, 1) +
  draw_plot(int_pred_realm, 1/3, 0, 1/3, 1) +
  draw_plot(mar_pred_realm, 2/3, 0, 1/3, 1) +
  draw_plot_label(label = c("a)", "b)", "c)"),
                  x = c(0, 1/3, 2/3),
                  y = c(0.99, 0.99, 0.99), size = 10, 
                  color = "grey30")

ggsave(pred_realm, path = "figures/term-paper/", 
       filename = "predictions-across-lat-realm.png", 
       device = "png", width = 11, height = 5)


## split range id into species and source:
uofill <- uofill %>%
  filter(realm == "Terrestrial") %>%
  mutate(species = paste(str_split_fixed(.$range, '_', 3)[,1], 
                         str_split_fixed(.$range, '_', 3)[,2], sep = ' ')) %>%
  mutate(source = str_split_fixed(.$range, '_', 3)[,3]) %>%
  select(range, species, source, everything()) %>%
  mutate(source = factor(.$source, levels = c("IUCN", "GBIF"), ordered = TRUE)) %>%
  arrange(species, type, source)
  
uofill_iucn <- uofill %>%
  mutate(temp = paste(species, type, sep = '')) %>%
  filter(!duplicated(temp)) %>%
  select(-temp)


## split by sensitivity analysis groups
types <- group_split(uofill_iucn, type)

elev_x_dormancy <- types[[1]]
Te_best <- types[[2]]
Te_best_subset <- Te_best[which(Te_best$species %in% Te_tpref$species),]
Te_tpref <- types[[3]]
Te_best_tpref <- filter(Te_best, Te_best$species %in% Te_tpref$species)
#dormancy <- types[[1]]
#elev_x_dormancy <- types[[2]]
#elevation <-types[[3]]
#reg <-types[[4]]

length(which(is.infinite(Te_sun$p_niche_upper))) ## 13 species have no potential thermal range

## prep dataset for figure making:
fig <- Te_best_subset %>%
  filter(realm == "Terrestrial") %>%
  filter(!is.infinite(r_niche_upper),
         !is.infinite(r_niche_lower),
         !is.infinite(p_niche_upper),
         !is.infinite(r_niche_lower)) %>%
  mutate(upper_inner = ifelse(warm_under < 0, "r_niche_upper", "p_niche_upper")) %>%
  mutate(lower_inner = ifelse(cold_under < 0, "r_niche_lower", "p_niche_lower")) %>%
  gather(key = "niche_limit_type", value = "limit_value", c(r_niche_upper, r_niche_lower, p_niche_upper, 
                                                            p_niche_lower, ctmax, ctmin))%>%
  mutate(r_p_c = ifelse(str_detect(niche_limit_type, "r_niche"), "r_niche", ifelse(str_detect(niche_limit_type, 
                                                                                       "p_niche"),"p_niche",
                        "ctlim"))) %>%
  mutate(w_c = ifelse(str_detect(niche_limit_type, "upper") | str_detect(niche_limit_type, "max"), 
                      "warm", "cold")) %>%
  mutate(line_group1 = ifelse((cold_under < 0) & (w_c == "cold") & (r_p_c != 'ctlim'), 'cold_under', 
                             ifelse((cold_over > 0) & (w_c == "cold") & (r_p_c != 'ctlim'), 'cold_over',
                                    ifelse((warm_under < 0) & (w_c == "warm") & (r_p_c != 'ctlim'), 'warm_under',
                                           ifelse((warm_over > 0) & (w_c == "warm") & (r_p_c != 'ctlim'),
                                                  "warm_over", NA))))) %>%
  mutate(line_type = ifelse(str_detect(line_group1, "over") & (r_p_c != "ctlim"), 
                            "Overfilling", ifelse(str_detect(line_group1, "under") & (r_p_c != "ctlim"), 
                                                  "Underfilling", NA))) %>%
  mutate(line_group2 = paste(range, line_group1)) %>%
  mutate(thin_line1 = ifelse((upper_inner == "r_niche_upper") & (niche_limit_type == "r_niche_upper"),"inner",
                            ifelse((upper_inner == "p_niche_upper") & (niche_limit_type == "p_niche_upper"),
                                   "inner", ifelse((lower_inner == "p_niche_lower") & 
                                                     (niche_limit_type == "p_niche_lower"), "inner", 
                                                   ifelse((lower_inner == "r_niche_lower") & 
                                                            (niche_limit_type == "r_niche_lower"), 
                                                          "inner", NA))))) %>%
  mutate(thin_line2 = paste(range, thin_line1)) 

fig_filling <- fig %>%
  filter(!is.na(thin_line1))

fig_overfilling <- fig %>%
  filter(!is.na(line_group1)) %>%
  filter(line_type == "Overfilling") 

fig_underfilling <- fig %>%
  filter(!is.na(line_group1)) %>%
  filter(line_type == "Underfilling") 

## plot filling across realms:
beauty <- ggplot(data = fig, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
  scale_colour_manual(values = c("steelblue", 'red4'), labels = c("Cold limit", "Warm limit")) +
  guides(color = FALSE) +
  geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
  geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
  geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
  scale_shape_manual(values=c(3, 1, 19), labels = c("Fundamental", "Potential",
                                                    "Realized")) + 
  labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
       x = "Realized range latitudinal midpoint (°N)") +
  geom_point() + 
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66))

ggsave(beauty, path = "figures/term-paper/", filename = "niche-filling-across-lat.png", 
       device = "png", width = 11, height = 6)

no_leg <- ggplot(data = fig, aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
  scale_colour_manual(values = c("steelblue", 'red4'), labels = c("Cold limit", "Warm limit")) +
  guides(color = FALSE, shape = FALSE, linetype = FALSE) +
  geom_line(data = fig_underfilling, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5) +
  geom_line(data = fig_overfilling, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1) +
  geom_line(data = fig_filling, aes(group = thin_line2), col = "grey80", size = 0.5) + 
  scale_shape_manual(values=c(3, 1, 19), labels = c("Fundamental", "Potential",
                                                    "Realized")) + 
  labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "Temperature (°C)", 
       x = "Realized range latitudinal midpoint (°N)") +
  geom_point() + 
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
  scale_y_continuous(limits = c(-42, 50))

fig_filling_t <- fig_filling %>%
  filter(realm == "terrestrial") 
  
fig_overfilling_t <- fig_overfilling %>%
  filter(realm == "terrestrial")

fig_underfilling_t <- fig_underfilling %>%
  filter(realm == "terrestrial")

terr <- fig %>%
  filter(realm == "terrestrial") %>%
  ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
  scale_colour_manual(values = c("steelblue", 'red4'), labels = c("Cold limit", "Warm limit")) +
  guides(color = FALSE, shape = FALSE, linetype = FALSE) +
  geom_line(data = fig_underfilling_t, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
  geom_line(data = fig_overfilling_t, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
  geom_line(data = fig_filling_t, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
  scale_shape_manual(values=c(3, 1, 19), labels = c("Fundamental", "Potential",
                                                    "Realized")) + 
  labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
       x = "") +
  geom_point(size = 0.5) + 
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
  scale_y_continuous(limits = c(-42, 50))

fig_filling_m <- fig_filling %>%
  filter(realm == "marine") 

fig_overfilling_m <- fig_overfilling %>%
  filter(realm == "marine")

fig_underfilling_m <- fig_underfilling %>%
  filter(realm == "marine")

marine <- fig %>%
  filter(realm == "marine") %>%
  ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
  scale_colour_manual(values = c("steelblue", 'red4'), labels = c("Cold limit", "Warm limit")) +
  guides(color = FALSE, shape = FALSE, linetype = FALSE) +
  geom_line(data = fig_underfilling_m, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
  geom_line(data = fig_overfilling_m, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
  geom_line(data = fig_filling_m, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
  scale_shape_manual(values=c(3, 1, 19), labels = c("Fundamental", "Potential",
                                                    "Realized")) + 
  labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
       x = "") +
  geom_point(size = 0.5) + 
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
  scale_y_continuous(limits = c(-42, 50))

fig_filling_i <- fig_filling %>%
  filter(realm == "intertidal") 

fig_overfilling_i <- fig_overfilling %>%
  filter(realm == "intertidal")

fig_underfilling_i <- fig_underfilling %>%
  filter(realm == "intertidal")

int <- fig %>%
  filter(realm == "intertidal") %>%
  ggplot(data = ., aes(x = lat_mp, y = limit_value, col = w_c, shape = r_p_c)) + 
  scale_colour_manual(values = c("steelblue", 'red4'), labels = c("Cold limit", "Warm limit")) +
  guides(color = FALSE, shape = FALSE, linetype = FALSE) +
  geom_line(data = fig_underfilling_i, aes(group = line_group2, col = w_c, linetype = "Underfilling"), size = 0.5*0.5) +
  geom_line(data = fig_overfilling_i, aes(group = line_group2, col = w_c, linetype = "Overfilling"), size = 1*0.5) +
  geom_line(data = fig_filling_i, aes(group = thin_line2), col = "grey80", size = 0.5*0.5) + 
  scale_shape_manual(values=c(3, 1, 19), labels = c("Fundamental", "Potential",
                                                    "Realized")) + 
  labs(shape = 'Thermal niche limit', linetype = "", col = '', y = "",
       x = "") +
  geom_point(size = 0.5) + 
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-90, 90, 30)), limits = c(-66, 66)) +
  scale_y_continuous(limits = c(-42, 50))


beauty_comb <- ggdraw() +
  draw_plot(no_leg, 0, 0, 0.75, 1) +
  draw_plot(marine, 0.75, 2/3, 0.25, 1/3) +
  draw_plot(int, 0.75, 1/3, 0.25, 1/3) +
  draw_plot(terr, 0.75, 0, 0.25, 1/3) +
  draw_plot_label(label = c("a)", "b)", "c)", "d)"),
                  x = c(0, 0.75, 0.75, 0.75),
                  y = c(0.99, 0.99, (2/3-0.01), (1/3-0.01)), size = 10, 
                  color = "grey30")

ggsave(beauty_comb, path = "figures/term-paper/", 
       filename = "niche-filling-across-lat-and-realm.png", 
       device = "png", width = 11, height = 5)


## distribution of filling vals:
filling <- Te_best %>%
  mutate(perfect_warm = ifelse(warm_over == 0 & warm_under == 0, 0, NA)) %>%
  mutate(perfect_cold = ifelse(cold_over == 0 & cold_under == 0, 0, NA)) %>%
  gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                          cold_over, perfect_warm, perfect_cold)) %>%
  filter(filling_value != 0 | filling_value == 0 & (filling_type %in% c("perfect_cold","perfect_warm"))) %>%
  select(range, species, source, realm, dormancy, filling_type, filling_value, lat_mp) %>%
  mutate(warm_cold= ifelse(str_detect(filling_type, "warm"), "warm", "cold")) 

perf_fill <- filling %>%
  filter(filling_value == 0)

warm_fig <- filling %>%
  filter(warm_cold == "warm") %>%
  filter(!is.infinite(filling_value)) %>%
  ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_segment(aes(x = abs(lat_mp), 
                                                                                xend = abs(lat_mp), 
                                                                                y = 0,
                                                                                yend = filling_value)) +
  scale_color_manual(values = c('red4')) + theme_bw() +
  guides(colour = FALSE) + 
  labs(y = "", x = "Absolute realized range latitudinal midpoint (°N)") + 
  scale_y_continuous(limits = c(-24, 45)) 
+
  geom_point(data = perf_fill, shape = 5)


filling %>%
  filter(warm_cold == "warm") %>%
  filter(!is.infinite(filling_value)) %>%
  ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_point() +
  scale_color_manual(values = c('red4')) + theme_bw() +
  guides(colour = FALSE) + 
  labs(y = "", x = "Absolute realized range latitudinal midpoint (°N/S)") + 
  scale_y_continuous(limits = c(-24, 10)) 

cold_fig <- filling %>%
  filter(warm_cold == "cold") %>%
  filter(!is.infinite(filling_value)) %>%
  ggplot(., aes(x = abs(lat_mp), y = filling_value, col = warm_cold)) + geom_segment(aes(x = abs(lat_mp), 
                                                                                       xend = abs(lat_mp), 
                                                                                       y = 0,
                                                                                       yend = filling_value)) +
  scale_color_manual(values = c("steelblue")) + theme_bw()  +
  guides(colour = FALSE) +
  labs(y = "Shortfall or excess of temperatures occupied (°C)", 
       x = "Absolute realized range latitudinal midpoint (°N)") + 
  scale_y_continuous(limits = c(-24, 45)) 

fig <- ggdraw() +
  draw_plot(cold_fig, 0, 0, 0.5, 1) +
  draw_plot(warm_fig, 0.5, 0, 0.5, 1) +
  draw_plot_label(label = c("a)", "b)"),
                  x = c(0, 0.5),
                  y = c(0.99, 0.99), size = 10, 
                  color = "grey30")

ggsave(fig, path = "figures/term-paper/", 
       filename = "niche-filling-distribution-cold-warm.png", 
       device = "png", width = 9, height = 5)




##########################################################
##                        garbage                       ##
##########################################################
## dispersal distance:
disp_cold <- fitted_pred_cold %>%
  group_by(dispersal_distance_continuous, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE))

disp_warm <- fitted_pred_warm %>%
  group_by(dispersal_distance_continuous, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE))

pred_disp_plot <- cont_disp %>%
  ggplot(., aes(x = dispersal_distance_continuous, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "Shortfall or excess of temperatures occupied (°C)",
       x = "Dispersal distance (km)")  +
  geom_line(data = disp_warm, inherit.aes = F, aes(x = dispersal_distance_continuous, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = disp_warm, inherit.aes = F, aes(x = dispersal_distance_continuous,
                                                     ymin = (filling_value-1.96*filling_value_SE),
                                                     ymax = filling_value+1.96*filling_value_SE, 
                                                     group = realm), alpha=0.1) +
  geom_line(data = disp_cold, inherit.aes = F, aes(x = dispersal_distance_continuous, y = filling_value), 
            col = "#0072B2") +
  geom_ribbon(data = disp_cold, inherit.aes = F, aes(x = dispersal_distance_continuous,
                                                     ymin = (filling_value-1.96*filling_value_SE),
                                                     ymax = filling_value+1.96*filling_value_SE, 
                                                     group = realm), alpha=0.1) +
  facet_wrap(~ realm) +
  theme_bw()


## body size:
bs_cold_mar <- fitted_pred_cold %>%
  filter(realm == "marine") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ > min(cont_disp_mar$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ < max(cont_disp_mar$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

bs_warm_mar <- fitted_pred_warm %>%
  filter(realm == "marine") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ > min(cont_disp_mar$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ < max(cont_disp_mar$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

mar_pred_bs <- cont_disp_mar %>%
  ggplot(., aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "Shortfall or excess of temperatures occupied (°C)",
       x = "Body length (cm)")  +
  geom_line(data = bs_warm_mar, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = bs_warm_mar, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  geom_line(data = bs_cold_mar, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#0072B2") +
  geom_ribbon(data = bs_cold_mar, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) 
#scale_y_continuous(limits = c(-28, 45))

bs_cold_int <- fitted_pred_cold %>%
  filter(realm == "intertidal") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ >= min(cont_disp_int$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ <= max(cont_disp_int$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

bs_warm_int <- fitted_pred_warm %>%
  filter(realm == "intertidal") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ >= min(cont_disp_int$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ <= max(cont_disp_int$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

int_pred_bs <- cont_disp_int %>%
  ggplot(., aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "Shortfall or excess of temperatures occupied (°C)",
       x = "Body length (cm)")  +
  geom_line(data = bs_warm_int, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = bs_warm_int, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  geom_line(data = bs_cold_int, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#0072B2") +
  geom_ribbon(data = bs_cold_int, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) 
#scale_y_continuous(limits = c(-28, 45))

bs_cold_ter <- fitted_pred_cold %>%
  filter(realm == "terrestrial") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ >= min(cont_disp_ter$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ <= max(cont_disp_ter$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

bs_warm_ter <- fitted_pred_warm %>%
  filter(realm == "terrestrial") %>%
  group_by(maximum_body_size_SVL_HBL_cm_, realm) %>%
  mutate(filling_value = mean(filling_value)) %>%
  mutate(filling_value_SE = mean(filling_value_SE)) %>%
  filter(maximum_body_size_SVL_HBL_cm_ >= min(cont_disp_ter$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE) & 
           maximum_body_size_SVL_HBL_cm_ <= max(cont_disp_ter$maximum_body_size_SVL_HBL_cm_, na.rm = TRUE))

ter_pred_bs <- cont_disp_ter %>%
  ggplot(., aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value, col = edge_type, shape = metric)) +
  geom_point() +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), labels = c("Cold niche edge", "Warm niche edge")) +
  labs(col = "", shape = "", y = "Shortfall or excess of temperatures occupied (°C)",
       x = "Body length (cm)")  +
  geom_line(data = bs_warm_ter, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#D55E00") +
  geom_ribbon(data = bs_warm_ter, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  geom_line(data = bs_cold_ter, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value), 
            col = "#0072B2") +
  geom_ribbon(data = bs_cold_ter, inherit.aes = F, aes(x = maximum_body_size_SVL_HBL_cm_,
                                                       ymin = (filling_value-1.96*filling_value_SE),
                                                       ymax = filling_value+1.96*filling_value_SE, 
                                                       group = realm), alpha=0.1) +
  theme_bw() +
  guides(colour = FALSE, shape = FALSE) 




cont_disp %>%
  ggplot(., aes(x = abs(lat_mp), y = filling_value, col = edge_type)) + geom_point() + 
  geom_smooth(method = "lm") 

cont_disp %>%
  filter(range_area_km2 < 30000000) %>%
  ggplot(., aes(x = range_area_km2, y = filling_value, col = edge_type)) + geom_point() + 
  geom_smooth(method = "lm") 

cont_disp %>%
  ggplot(., aes(x = dispersal_distance_continuous, y = filling_value, col = edge_type)) + geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = maximum_body_size_SVL_HBL_cm_, y = filling_value, col = edge_type)) + geom_point() + 
  geom_smooth(method = "lm") + 
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = hot_season_dormancy_, y = filling_value, col = edge_type)) + geom_boxplot() +
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = cold_season_dormancy_, y = filling_value, col = edge_type)) + geom_boxplot() + 
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = Trophic_position, y = filling_value, col = edge_type)) + geom_boxplot() + 
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = metric, y = filling_value, col = edge_type)) + geom_boxplot() + 
  facet_wrap(~ realm)

cont_disp %>%
  ggplot(., aes(x = dispersal_distance_continuous, y = filling_value, col = edge_type)) + geom_point() + 
  facet_wrap(~ dispersal_ability_category)


lims <- elev_x_dormancy %>%
  gather(key = "niche_limit_type", value = "limit_value", c(r_niche_upper, r_niche_lower, p_niche_upper, 
                                                            p_niche_lower)) %>%
  select(range, species, source, realm, dormancy, niche_limit_type, limit_value, lat_mp) %>%
  mutate(r_p = ifelse(str_detect(niche_limit_type, "r_niche"), "r_niche", "p_niche"))


## plot:
filling %>%
  ggplot(., aes(x = filling_value, fill = warm_cold)) + 
  geom_histogram(position = position_dodge())  

filling %>%
  ggplot(., aes(x = filling_value, fill = warm_cold)) + 
  geom_histogram(position = position_dodge())  +
  facet_wrap(~ realm)

lims %>%
  ggplot(., aes(x = limit_value, fill = r_p)) + 
  geom_histogram(position = position_dodge())  

lims %>%
  ggplot(., aes(x = limit_value, fill = r_p)) + 
  geom_histogram(position = position_dodge())  +
  facet_wrap(~ realm)

## make plot showing each spp. warm and cold over/underfilling
filling %>%
  filter(warm_cold == "cold") %>%
  ggplot(., aes(x = lat_mp, y = filling_value, fill = warm_cold)) + 
  geom_point()  

filling %>%
  filter(warm_cold == "cold") %>%
  ggplot(., aes(x = lat_mp, y = filling_value, fill = warm_cold)) + 
  geom_point() +
  facet_wrap(~ realm)

filling %>%
  filter(warm_cold == "warm") %>%
  ggplot(., aes(x = lat_mp, y = filling_value, fill = warm_cold)) + 
  geom_point()  

filling %>%
  filter(warm_cold == "warm") %>%
  ggplot(., aes(x = lat_mp, y = filling_value, fill = warm_cold)) + 
  geom_point() +
  facet_wrap(~ realm)


## okay, so species warm underfill and cold overfill 

## rearrange 

### breadth plot
breadth <- elev_x_dormancy %>%
  gather(key = "niche_limit_type", value = "limit_value", c(r_niche_upper, r_niche_lower, p_niche_upper, 
                                                            p_niche_lower)) %>%
  select(range, species, source, realm, dormancy, niche_limit_type, limit_value, range_area_km2, lat_mp) %>%
  mutate(r_p = ifelse(str_detect(niche_limit_type, "r_niche"), "r_niche", "p_niche")) %>%
  mutate(line_group = paste(range, r_p))
breadth %>%
  ggplot(., aes(x = lat_mp, y = limit_value, col = r_p)) + geom_point() +
  geom_line(aes(group = line_group))

clim_avail <- elev_x_dormancy %>%
  gather(key = "niche_limit_type", value = "limit_value", c(p_niche_upper, p_niche_lower, ctmin, ctmax)) %>%
  select(range, species, source, realm, dormancy, niche_limit_type, limit_value, range_area_km2, lat_mp)
clim_avail %>%
  ggplot(., aes(x = lat_mp, y = limit_value, col = niche_limit_type)) + geom_point()

filling <- elev_x_dormancy %>%
  gather(key = "filling_type", value = "filling_value", c(warm_under, warm_over, cold_under, 
                                                          cold_over)) %>%
  select(range, species, source, realm, dormancy, filling_type, filling_value, range_area_km2, lat_mp) %>%
  mutate(warm_cold= ifelse(str_detect(filling_type, "warm"), "warm", "cold")) %>%
  mutate(line_group = paste(range, warm_cold))
filling %>%
  ggplot(., aes(x = lat_mp, y = filling_value, col = warm_cold)) + geom_point() +
  geom_line(aes(group = line_group), position = position_dodge())


## see how potential niche and thermal limits compare (potential vs fundamental niche)
clim  <-  elev_x_dormancy %>%
  filter(!is.infinite(r_niche_upper),
         !is.infinite(r_niche_lower),
         !is.infinite(p_niche_upper),
         !is.infinite(r_niche_lower)) %>%
  mutate(avail_space_warm = p_niche_upper - ctmax) %>%
  mutate(avail_space_cold = ctmin - p_niche_lower) %>%
  gather(key = "space_type", value = "value", c(avail_space_warm, avail_space_cold)) 

clim_space <- clim %>%
  ggplot(data = ., aes(x = abs(lat_mp), y = abs(value), col = space_type)) + 
  scale_colour_manual(values = c("skyblue", 'red4')) +
  labs(shape = 'Niche limit type', col = '') +
  geom_point() +
  labs(col = '', y = "Degrees of unavailable niche space (°C)", x = "Realized range latitudinal midpoint") +
  theme_bw() +
  scale_x_continuous(breaks = c(seq(0, 90, 30)), 
                     labels = c("0°", "30°", "60°", "90°")) + 
  guides(color = FALSE)

## does ctmax increase with lat? no
fig %>%
  filter(r_p_c == "ctlim") %>%
  filter(w_c == "warm") %>%
  ggplot(., aes(x = abs(lat_mp), y = limit_value)) + geom_point()


## does ctmin increase with lat? yes
fig %>%
  filter(r_p_c == "ctlim") %>%
  filter(w_c == "cold") %>%
  ggplot(., aes(x = abs(lat_mp), y = limit_value)) + geom_point()




## explore cold and warm dornancy 
nodor <- as.data.frame(sens_types[[4]]) %>%
  mutate(dispersal_distance_continuous = ifelse(dispersal_distance_category == "0-1", 1, 
                                                ifelse(dispersal_distance_category == "1-10", 10,
                                                       ifelse(dispersal_distance_category == "10-100", 100,
                                                              ifelse(dispersal_distance_category == "100", 1000, 
                                                                     NA)))))


cont_disp %>%
  mutate(realm == "terrestrial") %>%
  filter(filling_type == "cold_under" | filling_type == "cold_over") %>%
  ggplot(., aes(x = cold_season_dormancy_, y = filling_value)) + 
  geom_boxplot()

nodor %>%
  mutate(realm == "terrestrial") %>%
  filter(filling_type == "cold_under" | filling_type == "cold_over") %>%
  ggplot(., aes(x = cold_season_dormancy_, y = filling_value)) + 
  geom_boxplot()

cont_disp %>%
  mutate(realm == "terrestrial") %>%
  filter(filling_type == "warm_under" | filling_type == "warm_over") %>%
  ggplot(., aes(x = hot_season_dormancy_, y = filling_value)) + 
  geom_boxplot()

nodor %>%
  mutate(realm == "terrestrial") %>%
  filter(filling_type == "warm_under" | filling_type == "warm_over") %>%
  ggplot(., aes(x = hot_season_dormancy_, y = filling_value)) + 
  geom_boxplot()
