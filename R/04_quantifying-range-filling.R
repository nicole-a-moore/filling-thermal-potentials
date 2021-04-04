## calculating potential range overfilling and potential range overpredicting for species with both thermal limits 
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(rnaturalearth)
library(fasterize)
library(rmapshaper)
select <- dplyr::select




## read in rasterized realized ranges and realized range shapefiles:
realized_ranges <- st_read("data-processed/realized-ranges_unsplit.shp") %>%
  mutate(range_id = paste(species, source, sep = "_"))

rasterized_rrs <- readRDS("data-processed/rasterized_rrs.rds")

## read in potential ranges:
potential_ranges <- readRDS("data-processed/potential_ranges_notcutatequator.rds")

range_filling <- data.frame(range = names(potential_ranges), equ_of = NA, pol_of = NA,
                            equ_op = NA, pol_op = NA, weighted_equ_of = NA, weighted_pol_of = NA,
                            weighted_equ_op = NA, weighted_pol_op = NA,
                            equ_of_prop = NA, pol_of_prop = NA,
                            equ_op_prop = NA, pol_op_prop = NA, 
                            weighted_equ_of_prop = NA, 
                            weighted_pol_of_prop = NA,
                            weighted_equ_op_prop = NA, 
                            weighted_pol_op_prop= NA )

i = 1
while (i < nlayers(potential_ranges) + 1) {
  ## for each pair of realized and potential range rasters
  pr <- potential_ranges[[i,]]
  rr <- rasterized_rrs[[which(names(rasterized_rrs) == names(pr))]] 
  
  ## get latitudinal midpoint of realized range and potential range:
  pr_pts <- rasterToPoints(pr, spatial = TRUE)
  pr_lmp <- (ymax(pr_pts) + ymin(pr_pts)) / 2
  rr_pts <- rasterToPoints(rr, spatial = TRUE)
  rr_lmp <- (ymax(rr_pts) + ymin(rr_pts)) / 2
  
  ## get minimum and maximum latitude of realized range and potential range:
  if (pr_lmp > 0) {
    pr_max <- ymax(pr_pts)
    pr_min <- ymin(pr_pts)
  }
  else {
    pr_max <- ymin(pr_pts)
    pr_min <- ymax(pr_pts)
  }
  
  if (rr_lmp > 0) {
    rr_max <- ymax(rr_pts)
    rr_min <- ymin(rr_pts)
  }
  else {
    rr_max <- ymin(rr_pts)
    rr_min <- ymax(rr_pts)
  }
  
  ## remove cells in potential range that overlap realized range and vice vera:
  masked_pr <- mask(pr, rr, updatevalue = NA, maskvalue = 1) ## set cells in pr that are overlapped by rr to NA
  masked_rr <- mask(rr, pr, updatevalue = NA, maskvalue = 1) ## set cells in rr that are overlapped by pr to NA
  
  # plot(rr)
  # plot(pr, add = TRUE, col = "red")
  # plot(masked_pr, add = TRUE, col = "blue")
  # plot(masked_rr, add = TRUE, col = "orange")
  
  ## CALCULATE REALIZED RANGE OVERFILLING:
  #########################################
  ## cut masked realized range at potential range latitudinal midpoint
  ## calculate equatorward and poleward non-overlapping area separately
  if (pr_lmp > 0) {
    equ <- Polygon(matrix(c(-180,0,-180,pr_lmp,180,pr_lmp,180,0,-180,0),
                        ncol=2, byrow=TRUE))
    pol <- Polygon(matrix(c(-180,pr_lmp,-180,90,180,90,180,pr_lmp,-180,pr_lmp),
                          ncol=2, byrow=TRUE))
  }
  else {
    equ <- Polygon(matrix(c(-180,pr_lmp,-180,0,180,0,180,pr_lmp,-180,pr_lmp),
                          ncol=2, byrow=TRUE))
    pol <- Polygon(matrix(c(-180,-90,-180,pr_lmp,180,pr_lmp,180,-90,-180,-90),
                          ncol=2, byrow=TRUE))
  }
  rect_equ <- SpatialPolygons(list(Polygons(list(equ), "p1"))) 
  rect_pol <- SpatialPolygons(list(Polygons(list(pol), "p1"))) 
  
  rect_raster_equ <- rasterize(rect_equ, r, background = NA, getCover = TRUE)
  rect_raster_equ[rect_raster_equ == 0] <- NA
  rect_raster_pol <- rasterize(rect_pol, r, background = NA, getCover = TRUE)
  rect_raster_pol[rect_raster_pol == 0] <- NA
  
  equ_of <- mask(masked_rr, rect_raster_equ)
  pol_of <- mask(masked_rr, rect_raster_pol)
  
  range_filling$equ_of[i] <- freq(equ_of, value = 1) 
  range_filling$pol_of[i] <- freq(pol_of, value = 1)
  
  ## get total number of cells in realized range above and below latitudinal midpoint
  ncell_rr_equ <- freq(mask(rr, rect_raster_equ), value = 1)
  ncell_rr_pol <- freq(mask(rr, rect_raster_pol), value = 1)
  
  ## add proportion version:
  range_filling$equ_of_prop[i] <- ifelse(ncell_rr_equ == 0, 0, 
                                         freq(equ_of, value = 1) / ncell_rr_equ)
  range_filling$pol_of_prop[i] <- ifelse(ncell_rr_pol ==0, 0, 
                                         freq(pol_of, value = 1) / ncell_rr_pol)
  
  ## add weighted version:
  ## weight warm overfilling by abs(lmp_pr - min lat of rr) 
  range_filling$weighted_equ_of[i] <- range_filling$equ_of[i] * abs(pr_lmp - rr_min)
  
  ## weight cold overfilling by abs(lmp_pr - max lat of rr) 
  range_filling$weighted_pol_of[i] <- range_filling$pol_of[i] * abs(pr_lmp - rr_max)
  
  ## add proportion version:
  range_filling$weighted_equ_of_prop[i] <- ifelse(ncell_rr_equ == 0, 0, 
                                                  (range_filling$equ_of[i]/ ncell_rr_equ) 
                                                  * abs(pr_lmp - rr_min))
  range_filling$weighted_pol_of_prop[i] <- ifelse(ncell_rr_pol == 0, 0, 
                                                  (range_filling$pol_of[i]/ ncell_rr_pol)
                                                  * abs(pr_lmp - rr_max))
  
  
  ## CALCULATE POTNETIAL RANGE OVERPREDICTING
  #########################################
  ## cut masked potential range at realized range latitudinal midpoint
  ## calculate equatorward and poleward non-overlapping area separately
  if (rr_lmp > 0) {
    equ <- Polygon(matrix(c(-180,0,-180,rr_lmp,180,rr_lmp,180,0,-180,0),
                          ncol=2, byrow=TRUE))
    pol <- Polygon(matrix(c(-180,rr_lmp,-180,90,180,90,180,rr_lmp,-180,rr_lmp),
                          ncol=2, byrow=TRUE))
  }
  else {
    equ <- Polygon(matrix(c(-180,rr_lmp,-180,0,180,0,180,rr_lmp,-180,rr_lmp),
                          ncol=2, byrow=TRUE))
    pol <- Polygon(matrix(c(-180,-90,-180,rr_lmp,180,rr_lmp,180,-90,-180,-90),
                          ncol=2, byrow=TRUE))
  }
  rect_equ <- SpatialPolygons(list(Polygons(list(equ), "p1"))) 
  rect_pol <- SpatialPolygons(list(Polygons(list(pol), "p1"))) 
  
  rect_raster_equ <- rasterize(rect_equ, r, background = NA, getCover = TRUE)
  rect_raster_equ[rect_raster_equ == 0] <- NA
  rect_raster_pol <- rasterize(rect_pol, r, background = NA, getCover = TRUE)
  rect_raster_pol[rect_raster_pol == 0] <- NA
  
  equ_op <- mask(masked_pr, rect_raster_equ)
  pol_op <- mask(masked_pr, rect_raster_pol)
  
  range_filling$equ_op[i] <- freq(equ_op, value = 1) 
  range_filling$pol_op[i] <- freq(pol_op, value = 1) 
  
  ## get total number of cells in realized range above and below latitudinal midpoint
  ncell_pr_equ <- freq(mask(pr, rect_raster_equ), value = 1)
  ncell_pr_pol <- freq(mask(pr, rect_raster_pol), value = 1)
  
  ## add proportion version:
  range_filling$equ_op_prop[i] <- ifelse(ncell_pr_equ == 0, 0, 
                                         freq(equ_op, value = 1) / ncell_pr_equ)
  range_filling$pol_op_prop[i] <- ifelse(ncell_pr_pol == 0, 0, 
                                         freq(pol_op, value = 1) / ncell_pr_pol)
  
  ## add weighted version:
  ## weight warm overpredicting by abs(rr_lmp - min lat of pr) 
  range_filling$weighted_equ_op[i] <- range_filling$equ_op[i] * abs(rr_lmp - pr_min)
  
  ## weight cold overpredicting by abs(rr_lmp - max lat of pr) 
  range_filling$weighted_pol_op[i] <- range_filling$pol_op[i] * abs(rr_lmp - pr_max)
  
  ## add proportion version:
  range_filling$weighted_equ_op_prop[i] <- ifelse(ncell_pr_equ == 0, 0, 
                                                  (range_filling$equ_op[i]/ ncell_pr_equ)
                                                  * abs(rr_lmp - pr_min))
  range_filling$weighted_pol_op_prop[i] <- ifelse(ncell_pr_pol == 0, 0, 
                                                  (range_filling$pol_op[i]/ ncell_pr_pol) * 
                                                    abs(rr_lmp - pr_max))
  
  print(paste("Finished range number:", i))
  i = i + 1 
}


range_filling <- range_filling %>%
  mutate(species = str_replace_all(as.character(range), "\\.", " ")) %>%
  mutate(species = str_split_fixed(species, "_", n = 2)[,1]) %>%
  mutate(source = str_split_fixed(range, "_", n = 2)[,2]) 

write.csv(range_filling, 'data-processed/range-filling-quantifications.csv', row.names = FALSE)


##range_filling <- read.csv('data-processed/range-filling-quantifications.csv')

## inspect metrics:
## arrange 
arranged_equ_of <- range_filling %>% 
  arrange(equ_of, pol_of, equ_op, pol_op)

arranged_pol_of <- range_filling %>% 
  arrange(pol_of, equ_of, equ_op, pol_op)

arranged_equ_op <- range_filling %>% 
  arrange(equ_op, pol_op, pol_of, equ_of)

arranged_pol_op <- range_filling %>% 
  arrange(pol_op, equ_op, pol_of, equ_of)


countries <- ne_countries(returnclass = "sf") %>%
  st_transform(., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")

## most useful to look at cases where:
## 1. a species with a greater metric has a lower weighted metric than another and other way around 
## 2. cases where weighting doesn't change much 
## 3. cases where weighting dramatically changes metric such that two that were different are very similar or opposite


## 1. 
## equ_of: 
## Periophthalmus.kalolo_GBIF	and	Nacella.concinna_GBIF
## Nucella.lapillus_GBIF	and	Divaricella.irpex_GBIF
## Sphenomorphus.indicus_GBIF	and Macoma.balthica_GBIF
## Cerastoderma.edule_GBIF	and	Oncorhynchus.tshawytscha_GBIF

## 2. 
## equ_of 
## Plagiotremus.azaleus_IUCN and Thalassoma.lucasanum_IUCN
## Anolis.gundlachi_GBIF	and	Eleutherodactylus.coqui_IUCN	and	Sphaerodactylus.nicholsi_IUCN	and	Sphaerodactylus.roosevelti_GBIF	and Sphaerodactylus.townsendi_IUCN	and Anolis.cooki_GBIF	and	Sphaerodactylus.roosevelti_IUCN (all have same metric for weighted and unweighted)


## 3.
## equ_of:
## Littorina.saxatilis_GBIF	and	Takydromus.sexlineatus_IUCN
## Halichoeres.dispilus_IUCNand Ameiva.festiva_GBIF
## Fundulus.parvipinnis_GBIF	and	Gehyra.variegata_GBIF	and	Heteronotia.binoei_IUCN (initially the same)
## Lophognathus gilberti	IUCN and Chamaeleo.dilepis_IUCN	and	Engraulis.japonicus_IUCN


rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")
prs <- as.data.frame(potential_ranges, xy=TRUE)
colnames(prs)[1:2] <- c("longitude", "latitude")

## ggplots of all:
i = 1
while (i < ncol(prs) - 1) {
range <- colnames(prs)[i+2]
stats <- range_filling[which(as.character(range_filling$range) == range),]
  
  r <- rrs[,c(1:2, which(colnames(rrs) == stats$range))] %>%
    rename("rr" = stats$range) %>%
    left_join(., prs[,c(1:2, which(colnames(prs) == stats$range))]) %>%
    rename("pr" = stats$range) 
  
  r_gg <- r %>%
    ggplot(., aes(x = longitude, y = latitude)) +
    xlim(-180, 180) + ylim(-90,90) + coord_fixed(ratio = 1) +
    geom_raster(aes(fill=as.factor(rr))) + 
    scale_fill_manual(values = c("yellow"), aesthetics = 'fill', labels = ) +
    annotate(geom="raster", x=r$longitude, y=r$latitude, alpha=.6,
             fill = r$pr) +
    labs(title = paste(stats$range, "\nOverfilling: warm = ", 
                       stats$equ_of, " , cold = ", stats$pol_of,
                       "\nWeighted overfilling: warm = ", stats$weighted_equ_of,
                       ", cold = ", stats$weighted_pol_of, 
                       "\nOverpredicing: warm = ", 
                       stats$equ_op, " , cold = ", stats$pol_op,
                       "\nWeighted overpredicting: warm = ", stats$weighted_equ_op,
                       ", cold = ", stats$weighted_pol_op),
         y = "Latitude",
         x = "Longitude") +
    scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180), expand = c(0.01,0.01)) +
    scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90), expand = c(0.01,0.01)) +
    annotate(geom="text", label = stats$equ_op, x = -180, y = -100) +
    theme(legend.position = "none")
  
  
  ## write to file:
    ggsave(r_gg, path = "figures/plotting-with-metrics/", 
           filename = paste(stats$range, ".png",sep = "_"), 
           height = 6, width = 10, units = "in", device = "png")
    
    i = i + 1
}








# garbage:



  ## pick out ones with the same values, look at weighted vs normal metrics in these cases 
  ## make a function to plot the ranges:
  
  range <- arranged_equ_of$range[2]

rr <- rasterized_rrs[[which(names(rasterized_rrs) == range),]]
pr <- potential_ranges[[which(names(potential_ranges) == range)]]

plot(rr)
plot(pr, add = TRUE, col = "orange")

rrs <- as.data.frame(rasterized_rrs, xy=TRUE)
colnames(rrs)[1:2] <- c("longitude", "latitude")

prs <- as.data.frame(potential_ranges, xy=TRUE)
colnames(prs)[1:2] <- c("longitude", "latitude")

rr <- rrs %>%
  select(1:2, 55) %>%
  drop_na()

pr <- prs %>%
  select(1:2, 55) %>%
  drop_na()

pr %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = colnames(pr)[3])) + 
  ggtitle("Realized range") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") 
