# this is R code put together by Tim Meehan and Nicole Michel to do bird
# friendliness index with ebird s&t maps


# setup ------------------------------------------------------------------------
# set up a bunch of stuff. load libraries, set options, define spatial 
# reference systems, etc.

# libraries
library(doParallel)
library(mapview)
library(sf)
library(ebirdst)
library(raster)
library(dplyr)

# options
options(scipen=9999999)
options(max.print=99999)
options(stringsAsFactors = F)

# pointers
extract <- raster::extract
select <- dplyr::select

# directories
wd <- "~/GitHub/QuantitativeMetrics/colombia_migrants/data"
stemd <- "Z:/10_Data/eBird/STEM/"
wildc <- "abundance_seasonal_nonbreeding.tif"
outd <- "Z:/5_MigratoryBirdInitiative/ColombiaNASA/BFI"
setwd(wd)

# crs's
mbi_an_crs <- "+proj=aeqd +lat_0=15 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ebird_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
wgs84_crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# rescale 01 raster
rescale_ras_01 <- function(x) {
  x.min <- minValue(x)
  x.max <- maxValue(x)
  y <- (x - x.min) / (x.max - x.min)
  return(y)
}
# ------------------------------------------------------------------------------


# get maps and extra data ------------------------------------------------------
# colombia
nat1 <- st_read("col_mainland.shp") %>% st_transform(ebird_crs) %>% 
  st_simplify(dTolerance = 100)
plot(nat1$geometry)

# extent
ext1 <- st_buffer(st_as_sfc(st_bbox(nat1), crs=ebird_crs), 10000)
plot(ext1, add=T)

# parks
prk1 <- st_read("PNN_59_2021.shp") %>% st_transform(ebird_crs) %>% 
  st_simplify(dTolerance = 100)
plot(prk1$geometry, add=T)
# ------------------------------------------------------------------------------


# species info -----------------------------------------------------------------
# species names
dir(stemd)
# grep("2019", list.files(stemd, pattern=wildc, recursive=T, full.names=T), value=T)

# species table
# read.csv("BFI_SpeciesList_Enduser_v1_09192022.csv") %>%
#   rename(scientific_name=spp) %>% left_join(ebirdst_runs %>% select(1,2,3)) %>%
#   left_join(read.csv("ACAD Global 2021.02.05.csv") %>%
#   select(scientific_name=Scientific.Name, common_name=Common.Name, ccs_max=CCS.max)) %>%
#   write.csv("species_table.csv", na="", row.names=F)

spp_tab <- read.csv("species_table.csv")

# bring in traits
# read traits table
traits <- read.csv("EltonTraits_BirdFuncDat.csv")

# fix name mismatches
spp_tab %>% distinct(English=common_name) %>% 
  left_join(traits %>% select(Scientific, English, SpecID)) %>% 
  filter(is.na(Scientific))
traits$English[traits$Scientific=="Catharus minimus"] <- "Gray-cheeked Thrush"
traits$English[traits$Scientific=="Charadrius alexandrinus"] <- "Snowy Plover"
traits$English[traits$Scientific=="Contopus virens"] <- "Eastern Wood-Pewee"
traits$English[traits$Scientific=="Gallinago gallinago"] <- "Wilson's Snipe"
traits$English[traits$Scientific=="Pluvialis squatarola"] <- "Black-bellied Plover"
traits$English[traits$Scientific=="Riparia riparia"] <- "Bank Swallow"
traits$English[traits$Scientific=="Tyrannus dominicensis"] <- "Gray Kingbird"

# filter traits for particular species and traits
traits <- traits %>% 
  select(common_name=English, 
         Diet.Inv, Diet.Vend, Diet.Vect, Diet.Vfish, Diet.Vunk, Diet.Scav, 
         Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO,
         ForStrat.watbelowsurf, ForStrat.wataroundsurf, 
         ForStrat.ground, ForStrat.understory, ForStrat.midhigh, 
         ForStrat.canopy, ForStrat.aerial, BodyMass.Value) %>% 
  filter(common_name %in% unique(spp_tab$common_name)) %>% 
  arrange(common_name) 
summary(traits)  

# fill in missing body masses from BotW
traits %>% select(common_name, BodyMass.Value) %>% filter(is.na(BodyMass.Value))
traits$BodyMass.Value[traits$common_name=="Rose-breasted Grosbeak"] <- (39 + 49) / 2

# clean
traits <- traits %>% select(-Diet.Vunk) %>% 
  rename_with(~tolower(gsub(".", "_", .x, fixed=T))) %>% 
  arrange(common_name)

# diss matirx
library(cluster)
library(NbClust)
names(traits)
d_gower <- cluster::daisy(traits[,-1], metric="gower",
                 weights=c(rep(1/4/9, 9), 
                           rep(1/2/7, 7),
                           rep(1/4/1, 1)))

# clusters
clust <- NbClust(diss=d_gower, distance=NULL, method="ward.D", 
                index="silhouette")

# extract best clusters
inddf <- traits %>% select(common_name) %>% 
  mutate(func_spec=clust$Best.partition) %>% 
  arrange(func_spec)
species_table <- spp_tab %>% left_join(inddf) %>% 
  mutate(func_spec=factor(LETTERS[func_spec]))

# save
# write.csv(species_table, "species_characteristics_table.csv", na="", row.names=F)

# reread
species_table <- read.csv("species_characteristics_table.csv")
# ------------------------------------------------------------------------------


# make bfi stack for colombia for winter ---------------------------------------
species_table <- read.csv("species_characteristics_table.csv") %>% 
  mutate(func_spec=as.numeric(factor(func_spec))) %>% 
  distinct() %>% 
  arrange(habitat, scientific_name)
habitats <- unique(species_table$habitat)
dirs1 <- list.dirs("Z:/10_Data/eBird/STEM")
sp_ext1 <- as(nat1, "Spatial")

bfi_stk <- stack()
for(h in 1:length(habitats)){
  # get building blocks
  hab <- habitats[h]
  dat <- species_table %>% filter(habitat==hab)
  sci_names <- dat %>% pull(scientific_name)
  com_names <- dat %>% pull(common_name)
  spp_codes <- dat %>% pull(species_code)
  
  # get info on grouped taxa
  not_empid <- grepl("wilfly", spp_codes)
  spp_codes <- spp_codes[!not_empid]
  not_cont <- grepl("eawpew", spp_codes)
  spp_codes <- spp_codes[!not_cont]
  add_empid <- any(not_empid==TRUE)
  add_cont <- any(not_cont==TRUE)
  
  # make empty stacks
  oc_stk <- stack()
  ab_stk <- stack()
  cs_stk <- stack()
  fs_stk <- stack()
  out_nams <- c()

  # make stacks with single taxa
  for(i in 1:length(spp_codes)){
    # spp info
    spp_i <- com_names[i]
    scode_i <- spp_codes[i]
    out_nams <- c(out_nams, scode_i)
    cs_i <- dat$ccs_max[i]
    fs_i <- dat$func_spec[i]
  
    # try catch
    tryCatch({
    # get abundance data file path
    dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", scode_i, "-ERD2018-EBIRD_SCIENCE"), 
         dirs1, value=T)
    dirs3 <- grep(paste0("results/tifs"), 
         dirs2, value=T)
    files1 <- list.files(dirs3, 
              pattern="2018_abundance_seasonal_nonbreeding.tif$")
    files2 <- paste0(dirs3, "/", files1)
    # make ab, oc, cs rasters
    ab_ras <- raster(files2)
    ab_ras <- crop(ab_ras, as(nat1, "Spatial"))
    ab_ras <- mask(ab_ras, as(nat1, "Spatial")); names(ab_ras) <- scode_i
    oc_ras <- ab_ras > 0; names(ab_ras) <- scode_i
    cs_ras <- oc_ras * cs_i; names(ab_ras) <- scode_i
    fs_ras <- oc_ras * as.numeric(fs_i); names(ab_ras) <- scode_i
    # add to stacks
    oc_stk <- addLayer(oc_stk, oc_ras)
    ab_stk <- addLayer(ab_stk, ab_ras)
    cs_stk <- addLayer(cs_stk, cs_ras)
    fs_stk <- addLayer(fs_stk, fs_ras)
    # wrap up
    print(paste("Done with", spp_i))
    }, error=function(e){})
  }
  
  # add empids
  if(add_empid==TRUE){
    empids <- c("wilfly", "aldfly")
    ab_stk2 <- stack()
    for(emp in 1:length(empids)){
      dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", empids[emp], "-ERD2018-EBIRD_SCIENCE"), 
           dirs1, value=T)
      dirs3 <- grep(paste0("results/tifs"), 
           dirs2, value=T)
      files1 <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_nonbreeding.tif$")
      files2 <- paste0(dirs3, "/", files1)
      
      ab_ras <- raster(files2)
      ab_ras <- crop(ab_ras, as(nat1, "Spatial"))
      ab_ras <- mask(ab_ras, as(nat1, "Spatial"))
      ab_stk2 <- addLayer(ab_stk2, ab_ras)
    }
    ab_ras2 <- ab_stk2[[1]] + ab_stk2[[2]]; names(ab_ras2) <- "empids"
    oc_ras2 <- ab_ras2 > 0; names(ab_ras2) <- "empids"
    cs_ras2 <- oc_ras2 * 10; names(ab_ras2) <- "empids"
    fs_ras2 <- oc_ras2 * 5; names(ab_ras2) <- "empids"
    
    # add empids to stacks
    ab_stk <- addLayer(ab_stk, ab_ras2)
    oc_stk <- addLayer(oc_stk, oc_ras2)
    cs_stk <- addLayer(cs_stk, cs_ras2)
    fs_stk <- addLayer(fs_stk, fs_ras2)
  }
    
  # add contopus
  if(add_cont==TRUE){
    empids <- c("eawpew", "wewpew")
    ab_stk2 <- stack()
    for(emp in 1:length(empids)){
      dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", empids[emp], "-ERD2018-EBIRD_SCIENCE"), 
           dirs1, value=T)
      dirs3 <- grep(paste0("results/tifs"), 
           dirs2, value=T)
      files1 <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_nonbreeding.tif$")
      files2 <- paste0(dirs3, "/", files1)
      ab_ras <- raster(files2)
      ab_ras <- crop(ab_ras, as(nat1, "Spatial"))
      ab_ras <- mask(ab_ras, as(nat1, "Spatial"))
      ab_stk2 <- addLayer(ab_stk2, ab_ras)
    }
    ab_ras2 <- ab_stk2[[1]] + ab_stk2[[2]]; names(ab_ras2) <- "contopus"
    oc_ras2 <- ab_ras2 > 0; names(ab_ras2) <- "contopus"
    cs_ras2 <- oc_ras2 * 12; names(ab_ras2) <- "contopus"
    fs_ras2 <- oc_ras2 * 5; names(ab_ras2) <- "contopus"
    
    # add contopus to stacks
    ab_stk <- addLayer(ab_stk, ab_ras2)
    oc_stk <- addLayer(oc_stk, oc_ras2)
    cs_stk <- addLayer(cs_stk, cs_ras2)
    fs_stk <- addLayer(fs_stk, fs_ras2)
  }
  
  # conservation weighted abundance per pixel
  ab_x_cs <- ab_stk * cs_stk
  
  # richness per pixel
  rich <- (sum(ab_x_cs>0)/10) + 1
  
  # # list of spp per func spec
  # dat$species_code[dat$species_code=="wilfly OR aldfly"] <- "empids"
  # dat$species_code[dat$species_code=="eawpew OR wewpew"] <- "contopus"
  # func_lst <- dat %>% group_by(func_spec) %>% 
  #   summarise(fs=paste(species_code, collapse = ","))
  
  # # get tot abundance per func spec
  # func_sums <- stack()
  # func_nams <- c()
  # for(r in 1:nrow(func_lst)){
  #   print(r)
  #   tryCatch({
  #   new_stk <- ab_stk[[unlist(strsplit(as.character(func_lst[r,2]), split=","))]]
  #   if(nlayers(new_stk)>1) func_sum <- sum(new_stk)
  #   if(nlayers(new_stk)==1) func_sum <- new_stk
  #   names(func_sum) <- func_lst[r, 2]
  #   func_sums <- addLayer(func_sums, func_sum)
  #   }, error=function(e){})
  # }
  
  # # calc func div per cell
  # library(vegan)
  # func_div <- calc(func_sums, diversity) 
  # func_div[func_div==0] <- min(func_div[func_div!=0])

  # calc bfi
  bfi_lay <- sum(ab_x_cs) * rich
  # bfi_lay <- calc(bfi_lay, rescale_ras_01)
  names(bfi_lay) <- hab
  bfi_stk <- addLayer(bfi_stk, bfi_lay)
}  

# rename
bfi_nonbreed <- stack(lapply(1:nlayers(bfi_stk), 
                             function(i){rescale_ras_01(bfi_stk[[i]])}))

# check
plot(bfi_nonbreed)

# output grids
habs <- unique(species_table$habitat)
for(f in 1:length(habs)){
  spp <- habs[f]
  # clipr::write_clip(
  # species_table %>% filter(habitat==spp) %>% arrange(func_spec)
  # )
  plot(((bfi_nonbreed[[spp]])))
  raster::writeRaster(bfi_nonbreed[[spp]],
                      paste0(outd, "/", spp, "_bfi_nonbreeding.tif"), overwrite=T)  
}
# ------------------------------------------------------------------------------


# make bfi stack for colombia for migration ------------------------------------
species_table <- read.csv("species_characteristics_table.csv") %>% 
  mutate(func_spec=as.numeric(factor(func_spec))) %>% 
  arrange(habitat, scientific_name)
habitats <- unique(species_table$habitat)
dirs1 <- list.dirs("Z:/10_Data/eBird/STEM")
sp_ext1 <- as(ext1, "Spatial")

bfi_stk <- stack()
for(h in 1:length(habitats)){
  # get building blocks
  hab <- habitats[h]
  dat <- species_table %>% filter(habitat==hab)
  sci_names <- dat %>% pull(scientific_name)
  com_names <- dat %>% pull(common_name)
  spp_codes <- dat %>% pull(species_code)
  
  # get info on grouped taxa
  not_empid <- grepl("wilfly", spp_codes)
  spp_codes <- spp_codes[!not_empid]
  not_cont <- grepl("eawpew", spp_codes)
  spp_codes <- spp_codes[!not_cont]
  add_empid <- any(not_empid==TRUE)
  add_cont <- any(not_cont==TRUE)
  
  # make empty stacks
  oc_stk <- stack()
  ab_stk <- stack()
  cs_stk <- stack()
  fs_stk <- stack()
  out_nams <- c()

  # make stacks with single taxa
  for(i in 1:length(spp_codes)){
    # spp info
    spp_i <- com_names[i]
    scode_i <- spp_codes[i]
    out_nams <- c(out_nams, scode_i)
    cs_i <- dat$ccs_max[i]
    fs_i <- dat$func_spec[i]
  
    # try catch
    tryCatch({
    # get abundance data file path
    dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", scode_i, "-ERD2018-EBIRD_SCIENCE"), 
         dirs1, value=T)
    dirs3 <- grep(paste0("results/tifs"), 
         dirs2, value=T)
    files1a <- list.files(dirs3, 
              pattern="2018_abundance_seasonal_prebreeding_migration.tif$")
    files2a <- paste0(dirs3, "/", files1a)
    files1b <- list.files(dirs3, 
              pattern="2018_abundance_seasonal_postbreeding_migration.tif$")
    files2b <- paste0(dirs3, "/", files1b)
    # make ab, oc, cs rasters
    ab_ras_a <- raster(files2a)
    ab_ras_b <- raster(files2b)
    ab_ras <- max(ab_ras_a, ab_ras_b, na.rm=T)
    ab_ras <- crop(ab_ras, sp_ext1); names(ab_ras) <- scode_i
    ab_ras <- mask(ab_ras, as(nat1, "Spatial"))
    oc_ras <- ab_ras > 0; names(ab_ras) <- scode_i
    cs_ras <- oc_ras * cs_i; names(ab_ras) <- scode_i
    fs_ras <- oc_ras * as.numeric(fs_i); names(ab_ras) <- scode_i
    # add to stacks
    oc_stk <- addLayer(oc_stk, oc_ras)
    ab_stk <- addLayer(ab_stk, ab_ras)
    cs_stk <- addLayer(cs_stk, cs_ras)
    fs_stk <- addLayer(fs_stk, fs_ras)
    # wrap up
    print(paste("Done with", spp_i))
    }, error=function(e){})
  }
  
  # add empids
  if(add_empid==TRUE){
    empids <- c("wilfly", "aldfly")
    ab_stk2 <- stack()
    for(emp in 1:length(empids)){
      dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", empids[emp], "-ERD2018-EBIRD_SCIENCE"), 
           dirs1, value=T)
      dirs3 <- grep(paste0("results/tifs"), 
           dirs2, value=T)
      files1a <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_prebreeding_migration.tif$")
      files2a <- paste0(dirs3, "/", files1a)
      files1b <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_postbreeding_migration.tif$")
      files2b <- paste0(dirs3, "/", files1b)
      ab_ras_a <- raster(files2a)
      ab_ras_b <- raster(files2b)
      ab_ras <- max(ab_ras_a, ab_ras_b, na.rm=T)
      ab_ras <- crop(ab_ras, sp_ext1)
      ab_ras <- mask(ab_ras, as(nat1, "Spatial"))
      ab_stk2 <- addLayer(ab_stk2, ab_ras)
    }
    ab_ras2 <- ab_stk2[[1]] + ab_stk2[[2]]; names(ab_ras2) <- "empids"
    oc_ras2 <- ab_ras2 > 0; names(ab_ras2) <- "empids"
    cs_ras2 <- oc_ras2 * 10; names(ab_ras2) <- "empids"
    fs_ras2 <- oc_ras2 * 5; names(ab_ras2) <- "empids"
    
    # add empids to stacks
    ab_stk <- addLayer(ab_stk, ab_ras2)
    oc_stk <- addLayer(oc_stk, oc_ras2)
    cs_stk <- addLayer(cs_stk, cs_ras2)
    fs_stk <- addLayer(fs_stk, fs_ras2)
  }
    
  # add contopus
  if(add_cont==TRUE){
    empids <- c("eawpew", "wewpew")
    ab_stk2 <- stack()
    for(emp in 1:length(empids)){
      dirs2 <- grep(paste0("Z:/10_Data/eBird/STEM/", empids[emp], "-ERD2018-EBIRD_SCIENCE"), 
           dirs1, value=T)
      dirs3 <- grep(paste0("results/tifs"), 
           dirs2, value=T)
      files1a <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_prebreeding_migration.tif$")
      files2a <- paste0(dirs3, "/", files1a)
      files1b <- list.files(dirs3, 
                pattern="2018_abundance_seasonal_postbreeding_migration.tif$")
      files2b <- paste0(dirs3, "/", files1b)
      ab_ras_a <- raster(files2a)
      ab_ras_b <- raster(files2b)
      ab_ras <- max(ab_ras_a, ab_ras_b, na.rm=T)
      ab_ras <- crop(ab_ras, sp_ext1)
      ab_ras <- mask(ab_ras, as(nat1, "Spatial"))
      ab_stk2 <- addLayer(ab_stk2, ab_ras)
    }
    ab_ras2 <- ab_stk2[[1]] + ab_stk2[[2]]; names(ab_ras2) <- "contopus"
    oc_ras2 <- ab_ras2 > 0; names(ab_ras2) <- "contopus"
    cs_ras2 <- oc_ras2 * 12; names(ab_ras2) <- "contopus"
    fs_ras2 <- oc_ras2 * 5; names(ab_ras2) <- "contopus"
    
    # add contopus to stacks
    ab_stk <- addLayer(ab_stk, ab_ras2)
    oc_stk <- addLayer(oc_stk, oc_ras2)
    cs_stk <- addLayer(cs_stk, cs_ras2)
    fs_stk <- addLayer(fs_stk, fs_ras2)
  }
  
  # conservation weighted abundance
  ab_x_cs <- ab_stk * cs_stk
  
  # richness per pixel
  rich <- (sum(ab_x_cs>0)/10) + 1
  
  # # list of spp per func spec
  # dat$species_code[dat$species_code=="wilfly OR aldfly"] <- "empids"
  # dat$species_code[dat$species_code=="eawpew OR wewpew"] <- "contopus"
  # func_lst <- dat %>% group_by(func_spec) %>% 
  #   summarise(fs=paste(species_code, collapse = ","))
  
  # # get tot abundance per func spec
  # func_sums <- stack()
  # func_nams <- c()
  # for(r in 1:nrow(func_lst)){
  #   print(r)
  #   tryCatch({
  #   new_stk <- ab_stk[[unlist(strsplit(as.character(func_lst[r,2]), split=","))]]
  #   if(nlayers(new_stk)>1) func_sum <- sum(new_stk)
  #   if(nlayers(new_stk)==1) func_sum <- new_stk
  #   names(func_sum) <- func_lst[r, 2]
  #   func_sums <- addLayer(func_sums, func_sum)
  #   }, error=function(e){})
  # }
  
  # # calc func div per cell
  # library(vegan)
  # func_div <- calc(func_sums, diversity) 
  # func_div[func_div==0] <- min(func_div[func_div!=0])

  # calc bfi
  bfi_lay <- sum(ab_x_cs) * rich
  # bfi_lay <- calc(bfi_lay, rescale_ras_01)
  names(bfi_lay) <- hab
  bfi_stk <- addLayer(bfi_stk, bfi_lay)
}  

# rename
bfi_migration <- stack(lapply(1:nlayers(bfi_stk), 
                             function(i){rescale_ras_01(bfi_stk[[i]])}))

# check
plot(bfi_migration)

# output grids
habs <- unique(species_table$habitat)
for(f in 1:length(habs)){
  spp <- habs[f]
  # clipr::write_clip(
  # species_table %>% filter(habitat==spp) %>% arrange(func_spec)
  # )
  plot(((bfi_migration[[spp]])))
  raster::writeRaster(bfi_migration[[spp]],
                      paste0(outd, "/", spp, "_bfi_migration.tif"), overwrite=T)  
}
# ------------------------------------------------------------------------------



# max layers -------------------------------------------------------------------
# output grids
habs <- unique(species_table$habitat)
bfi_migration <- crop(bfi_migration, bfi_nonbreed)
for(f in 1:length(habs)){
  spp <- habs[f]
  plot(maxras <- max(stack(bfi_migration[[spp]], bfi_nonbreed[[spp]])))
  raster::writeRaster(maxras,
                    paste0(outd, "/", spp, "_bfi_max.tif"), overwrite=T)
}
# ------------------------------------------------------------------------------



# write shapes -----------------------------------------------------------------
st_write(nat1, paste0(outd, "/", "national_border.shp"))
st_write(ext1, paste0(outd, "/", "national_extent.shp"))
# ------------------------------------------------------------------------------