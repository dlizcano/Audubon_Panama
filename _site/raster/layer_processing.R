setwd("/Users/jorge/Documents/Projects/Panama/Environmental_layers")
#create mask
mgr <-st_read("./mangroves_parita.shp")
mgr.utm <- st_transform(mgr, crs=32617)
mgr.utm500 <- st_buffer(mgr.utm, dist=500)
mgr.500 <- st_transform(mgr.utm500, crs=4326)
st_write(mgr.500, "./mangroves_parita_500m.shp")

#Building density
bg225 <- raster("/Users/jorge/Documents/Data/buildings/EarthEngine/cell_225.tif")
bg224 <- raster("/Users/jorge/Documents/Data/buildings/EarthEngine/cell_224.tif")
bg_parita <- mosaic(bg225,bg224, fun=mean)
bg_parita <- bg_parita==0
out.rst <- crop(bg_parita, mgr.500)
writeRaster(out.rst, "/Users/jorge/Documents/Projects/Panama/Environmental_layers/building_density.tif",
            overwrite=T)

#finish mask
mgr_rst <- terra::rasterize(mgr.500, out.rst, cover=T)
mgr_rst <- mgr_rst >0
writeRaster(mgr_rst, "/Users/jorge/Documents/Projects/Panama/Environmental_layers/mangroves_parita.tif", overwrite=T)
rm(list=ls())

template <- raster("mangroves_parita.tif")

#Human footprint
hfp <- raster("/Users/jorge/Documents/Data/HumanFootprint2009/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif")
plot(hfp)
tmp <- drawExtent(show=TRUE, col="red")
hfp <- crop(hfp, tmp)
hfp2 <- projectRaster(hfp, crs=projection(template), method="ngb")
hfp2[hfp2==128] <- NA
out.rst <- resample(hfp2, template)
writeRaster(out.rst, "/Users/jorge/Documents/Projects/Panama/Environmental_layers/human_footprint.tif", overwrite=T)




#Forest integrity
fii <- raster("flii_NorthAmerica.tiff")
plot(fii)
tmp <- drawExtent(show=TRUE, col="red")
fii <- crop(fii, tmp)
fii[fii==-9999] <- NA
out.rst <- resample(fii, template)
out.rst <- out.rst/1000
writeRaster(out.rst, "/Users/jorge/Documents/Projects/Panama/Environmental_layers/forest_integrity_index.tif",
            overwrite=T)

#distance to buildings
bg <- raster("building_density.tif")
bg[bg==0] <- NA
bg_dist <- terra::distance(bg, grid=T)
writeRaster(bg_dist, "buildings_distance.tif", overwrite=T)

#Resample all layers to the same extent and resolution
setwd("/Users/jorge/Documents/Projects/Panama/Environmental_layers/")
dir.create("/Users/jorge/Documents/Projects/Panama/Environmental_layers/resample_30m")
aoi <- raster("mask_final_500m_v2.tiff")

#worldclim
rst.list <- list.files("/Users/jorge/Documents/Data/CuboTNC/Worldclim/NewExtent","*.tif$", full.names=T)
for(i in 1:length(rst.list)){
  in.rst <- raster(rst.list[i])
  out.rst <- resample(in.rst, aoi)
  out.rst <- mask(out.rst, aoi)
  writeRaster(out.rst, paste0("./resample_30m/", basename(rst.list[i]), ".tif"))
}

#soils
rst.list <- list.files("./soils","*.tiff$", full.names=T)
for(i in 1:length(rst.list)){
  in.rst <- raster(rst.list[i])
  out.rst <- resample(in.rst, aoi)
  out.rst <- mask(out.rst, aoi)
  writeRaster(out.rst, paste0("./resample_30m/", basename(rst.list[i]), ".tif"))
}

#redo soil groups which is categorical using ngb
in.rst <- raster("./soils/soil_group.tiff")
out.rst <- resample(in.rst, aoi, method="ngb")
out.rst <- mask(out.rst, aoi)
writeRaster(out.rst, paste0("./resample_30m/", basename(rst.list[i]), ".tif"), overwrite=T)

resample2 <- function(rst.filename, aoi, method, out.path){
  in.rst <- raster(rst.filename)
  out.rst <- resample(in.rst, aoi, method=method)
  out.rst <- mask(out.rst, aoi)
  writeRaster(out.rst, paste0(out.path, basename(rst.filename), ".tif"))
  plot(out.rst)
}

#carbon
resample2("AGB Spawn.tiff", aoi, "bilinear", "./resample_30m/")
resample2("BGB Spawn.tiff", aoi, "bilinear", "./resample_30m/")
resample2("buildings_distance.tif", aoi, "bilinear", "./resample_30m/")
resample2("canopy_height_jetz2.tiff", aoi, "bilinear", "./resample_30m/")
resample2("coast_final_v2.tiff", aoi, "bilinear", "./resample_30m/")
resample2("forest_integrity_index.tif", aoi, "bilinear", "./resample_30m/")
resample2("human_footprint.tif", aoi, "bilinear", "./resample_30m/")
resample2("rivers_final_v2.tiff", aoi, "bilinear", "./resample_30m/")
resample2("roads_final_v2.tiff", aoi, "bilinear", "./resample_30m/")
resample2("S2_NDVI_median_v2.tiff", aoi, "bilinear", "./resample_30m/")

#Sampled locations
rm(list=ls())
rst.list <- list.files("./resample_30m","*.tif$", full.names=T)
allvars <- stack(rst.list)
sites <- read.csv("sampled_sites.csv", as.is=T)
coordinates(sites) <- ~Longitude + Latitude
sites.covs <- extract(allvars, sites)
sites.covs <- cbind(sites, sites.covs)
rgdal::writeOGR(sites, getwd(), "sampled_sites", driver="ESRI Shapefile")
write.csv(sites.covs, "sites.covs.csv", row.names=FALSE)

i=47
plot(density(na.omit(values(allvars[[i]]))), main=names(allvars[[i]]))
points(x=sites.covs[,i], y=rep(1e-04,nrow(sites.covs)), col="red")
