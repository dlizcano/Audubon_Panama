
# load function
source("C:/CodigoR/AudubonPanama/R/matrix_creator.R")
############### load data

library(readr)
library(tidyverse)
parita_data_full <- read_delim("C:/CodigoR/AudubonPanama/data/all_records_audubon_panama_to_occu.csv", 
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(eventDate = col_date(format = "%Y-%m-%d")), 
                               trim_ws = TRUE)

## Nyctidromus albicollis
## Aramus guarauna
## Buteogallus anthracinus
## Dryocopus lineatus
## Piranga rubra

# filter by threshold
parita_data_5sp <-  parita_data_full %>% 
                      filter(confidence >= threshold) %>%
                      filter(scientificName == c("Nyctidromus albicollis", 
                                                 "Aramus guarauna",
                                                 "Buteogallus anthracinus",
                                                 "Dryocopus lineatus",
                                                 "Piranga rubra"
                                                 ))  
  


# insert Max and min 
# start_date
parita_data2 <- parita_data_5sp %>% 
  group_by(locationID) %>% 
  mutate(start_date = min(eventDate),
         end_date= max(eventDate))

parita <- f.matrix.creator2(data=parita_data2, year=2023)
covs <- read_csv("C:/CodigoR/AudubonPanama/shp/sites_covs_parita_nona.csv")
library(unmarked)
library(sf)
library(mapview)
library(ubms)
library(stocc)
library(terra)
library(raster)


# centroide for terra
# centoide <- centroids(puntos, TRUE)
centroide <- c(mean(as.numeric(parita_data_5sp$decimalLongitude)), mean(as.numeric(parita_data_5sp$decimalLatitude)))
clip_window <- extent(-80.67105, -80.05822, 7.713744, 8.489888) # 7.932311, -80.612233  8.360073, -80.180109
bb <- c(-80.612233, -80.180109, 7.932311, 8.360073)
srtm <- raster::getData('SRTM', lon=centroide[1], lat=centroide[2])

# crop the  raster using the vector extent
srtm_crop <- raster::crop(srtm, clip_window)
 #rast(srtm_crop) # convert to terra


# altitud <- elevation_3s(-72.893262, 7.664081007, path="data")
# altitud <- as.numeric(Camptostoma_obsoletum_occu_dia$Altitude)

# convierte a terra
puntos <- vect(covs, geom=c("Longitude", "Latitude"), crs="EPSG:4326")
# convierte a sf
puntos_sf <- sf::st_as_sf(puntos)
# extract elev from raster
cam_covs <- raster::extract(srtm_crop, puntos_sf)


AGB_Spawn <- rast("C:/CodigoR/AudubonPanama/raster/AGB Spawn.tif")
NDVI <- rast("C:/CodigoR/AudubonPanama/raster/S2_NDVI_median_v2.tif")
roads <- rast("C:/CodigoR/AudubonPanama/raster/roads_final_v2.tif")
carbon_stock <- rast("C:/CodigoR/AudubonPanama/raster/soil_organic_carbon_stock_0-30m.tif")

srtm_projected <- projectRaster(srtm_crop, crs=projection(NDVI), method="ngb") 
elev <- resample(rast(srtm_projected), NDVI)




# list of raster terras
many_rasters <- list(AGB_Spawn, NDVI, roads, carbon_stock, elev)
# terra stack
covs_raster <- rast(many_rasters)
covs_raster <- subst(covs_raster, NA, 0)# change NA to 0 
covs_raster_val <- extract(covs_raster, puntos_sf)

# Sp names
names(parita$dete[]) # dete es hora occur es presencia el dia

detach("package:raster", unload = TRUE)# elimina raster

umf_y_full_1<- unmarkedFrameOccu(y= as.data.frame(parita$occur[[1]])) #"Aramus guarauna"
siteCovs(umf_y_full_1) <- data.frame(elevation =  scale(as.numeric(covs_raster_val$srtm_20_11)),
                                   AGB_Spawn = scale(as.numeric(covs_raster_val$`AGB Spawn`)),
                                   roads = scale(as.numeric(covs_raster_val$roads_final_v2)),
                                   carbon = scale(as.numeric(covs_raster_val$`soil_organic_carbon_stock_0-30m`)),
                                   NDVI=scale(covs_raster_val$S2_NDVI_median_v2)
                                   ) # data.frame(Elev=full_covs$Elev) # Full

plot(umf_y_full_1)  #"Aramus guarauna" no converge

umf_y_full_2<- unmarkedFrameOccu(y= as.data.frame(parita$occur[[2]])) #"Aramus guarauna"
siteCovs(umf_y_full_2) <- data.frame(elevation =  covs_raster_val$srtm_20_11,
                                   AGB_Spawn = covs_raster_val$`AGB Spawn`,
                                   roads = covs_raster_val$roads_final_v2,
                                   carbon = covs_raster_val$`soil_organic_carbon_stock_0-30m`,
                                   NDVI=covs_raster_val$S2_NDVI_median_v2
) # data.frame(Elev=full_covs$Elev) # Full


plot(umf_y_full_2, title="Buteogallus anthracinus")  #""Buteogallus anthracinus no converge


obsCovs(umf_y_full_2) <- list(hour= matrix( parita$dete[[2]])) # hora

# Double formula: first part is for detection, second for occupancy
#  form <- ~detCov1 + detCov2 ~habCov1 + habCov2 + RSR(x, y, threshold=1)

# fit unmarked
fit_unm_sp2_0 <- unmarked::occu(~1~1, data=umf_y_full_2) # ok!
fit_unm_sp2_1 <- unmarked::occu(~1~scale(elevation), data=umf_y_full_2)


# fit stan
fit_stan_sp2_0 <- stan_occu(~1~1, data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_0h <- stan_occu(~scale(hour)~1, data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_1 <- stan_occu(~scale(hour)~scale(elevation), data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_2 <- stan_occu(~scale(hour)~scale(AGB_Spawn), data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_3 <- stan_occu(~scale(hour)~scale(roads), data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_4 <- stan_occu(~scale(hour)~scale(carbon), data=umf_y_full_2, chains=3, iter=500, cores=3)
fit_stan_sp2_5 <- stan_occu(~scale(hour)~scale(NDVI), data=umf_y_full_2, chains=3, iter=500, cores=3)


mods <- fitList(fit_stan_sp2_0,
                fit_stan_sp2_0h,
                fit_stan_sp2_1,
                fit_stan_sp2_2,
                fit_stan_sp2_3,
                fit_stan_sp2_4,
                fit_stan_sp2_5)

# model selection
round(modSel(mods), 3)

# plot top model
plot_effects(fit_stan_sp2_1, "det") # Detectiom
plot_effects(fit_stan_sp2_1, "state") # Ocupancy

# head(predict(fit_top, submodel="state"))

#### predict

library(RColorBrewer)
library(raster)
library(rasterVis)
library(mapview)

srtm_crop_s <- mask(elev, NDVI)#stack(raster(elev))#, 
#scale(roughness)) # scale altitud
names(srtm_crop_s) <- c("elevation")#, "roughness")
# crs(srtm_crop_s) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
nd <- srtm_crop_s #data.frame(srtm_crop_s, hour=0.5)
newd <- data.frame(srtm_crop_s, hour=0.5)
pred_psi_s_unm <-predict(fit_unm_sp2_1, type="state", newdata=nd)
#pred_psi_s_stan <-predict(fit_stan_sp2_1, submodel="state", newdata=nd, na.rm=FALSE)
#pred_psi_ras <- raster(srtm_crop_s) # make raster = to elev 
#pred_psi_ras[] <- pred_psi_s_stan$Predicted # drape on top


# pred_psi_r <- pred_psi_s # * sd(full_covs$elevation) + mean(full_covs$elevation)
# crs(pred_psi_ras) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
clr <- colorRampPalette(brewer.pal(9, "YlGn"))
# mapview (pred_psi_r[[1]], col.regions = clr,  legend = TRUE, alpha=0.7)
# plot(pred_psi_s[[1]], main="Occupancy")
levelplot(pred_psi_s[[1]], par.settings = YlOrRdTheme(), margin=FALSE, main="Ocupancy")

pal = mapviewPalette("mapviewSpectralColors")
mapview(raster(pred_psi_s[[1]]), map.types = c("Esri.WorldImagery"), alpha=0.5, col.regions = pal(10)) + mapview(puntos_sf)




