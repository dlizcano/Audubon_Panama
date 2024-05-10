# set up -----------------------------------------------------------------------
library(vegan)
library(unmarked)
library(mgcv)
library(mgcViz)
library(tidyverse)
library(lubridate)
library(stringr)

setwd("C:/Users/tmeehan/Box/Cauca Valley monitoring/r_work")
dir()
# ------------------------------------------------------------------------------


# get raw count data -----------------------------------------------------------
d0 <- read.csv("BASE DE DATOS FINAL AUDUBON.csv") %>% 
  rename_all(tolower) %>% 
  rename_all(~ gsub(".", "_", ., fixed=T)) %>% 
  filter(tipo_de_registro %in% c("Visual","Auditivo","visual","auditivo")) %>%
  mutate(distance=as.numeric(stringr::str_trim(distancia))) %>% 
  mutate(date=dmy(paste0(str_trim(fecha), "-2023"))) %>% 
  mutate(hora_inicio=hm(hora_inicio),
         hora_final=hm(hora_final)) %>% 
  mutate(start_time=hour(hora_inicio)+(minute(hora_inicio)/60),
         end_time=hour(hora_final)+(minute(hora_final)/60),
         month=month(date)) %>% 
  mutate(count=as.numeric(str_trim(individuos)),
         observador=str_trim(observador),
         cuenca=str_trim(cuenca),
         celda=str_trim(celda),
         punto=str_trim(punto)) %>% 
  filter(!is.na(count)) %>% 
  select(watershed=cuenca, site=celda, point=punto, date, month, 
         observer_id=observador, start_time, end_time, species=especie, 
         distance, count) %>% 
  uncount(count) %>% 
  mutate(month=ifelse(month==9, 8, month))
head(d0)

# check some sampling properties
d0 %>% mutate(site=as.numeric(site)) %>% distinct(site) %>% 
  arrange(site) # up to 8 sites
d0 %>% mutate(point=as.numeric(point)) %>% distinct(point) %>% 
  arrange(point) # up to 13 points per site
d0 %>% distinct(date) %>% arrange(date) # 40 dates
d0 %>% distinct(date, month) %>% arrange(month, date) %>% 
  group_by(month) %>% count() # 3 to 8 dates per month
d0 %>% distinct(site, month) %>% arrange(month, site) %>% 
  group_by(month) %>% count() # 3 to 8 sites per month
d0 %>% distinct(site, point, date, observer_id) %>% 
  mutate(point=as.numeric(point)) %>% 
  arrange(date, site, point, observer_id) %>% 
  group_by(date, site, point) %>% count() %>% 
  View() # almost always two observers per count

# just landbird species with trait data
spp1 <- read.csv("species_detect_diet_acad_forage_trait_data_022024.csv") %>% 
  filter(guild=="landbird")
d1 <- d0 %>% filter(species %in% spp1$species)
# ------------------------------------------------------------------------------



# distance analysis ------------------------------------------------------------
# filter out no distance data and distance over 100
max_d <- 80
d2 <- d1 %>% 
  filter(!is.na(distance),
         distance<=max_d)

# just species with large sample more than 20
thresh_ss <- 20
d2 <- d2 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
d2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(d2$species))
spp_nams <- unique(d2$species)

# make counts per bin
d3 <- d2 %>% 
  mutate(distance_bin=ifelse(distance<max_d/4, max_d/4, ifelse(
    distance>=max_d/4 & distance<max_d/2, max_d/2, ifelse(
      distance>=max_d/2 & distance<max_d*3/4, max_d*3/4, max_d)))) %>% 
  group_by(site, point, month, species, distance_bin) %>% 
  count() 

# zero fill
all_combos <- expand_grid(site=unique(d3$site),
                          point=unique(d3$point),
                          month=unique(d3$month),
                          species=unique(d3$species),
                          distance_bin=unique(d3$distance_bin))
d4 <- all_combos %>% left_join(d3) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(distance_bin, species, n)

# loop through species
dist_out <- c()
for(i in 1:length(spp_nams)){
  # data
  d_i <- d4 %>% filter(species==spp_nams[i]) %>% 
    pivot_wider(values_from=n, names_from=distance_bin)
  # umf
  umf_i <- unmarkedFrameDS(y=as.matrix(d_i[,5:8]),
                            #siteCovs=data.frame(),
                            dist.breaks=c(0,max_d/4,max_d/2,max_d*3/4, max_d), 
                            unitsIn="m", survey="point")
  # model
  dm_i <- distsamp(~1 ~1, umf_i)
  # effective area
  sig <- exp(coef(dm_i, type="det"))
  ea <- 2*pi * integrate(grhn, 0, max_d, sigma=sig)$value
  # detection probability
  dp_i <- ea / (pi*max_d^2)
  df_i <- data.frame(species=spp_nams[i], dist_det_p=dp_i)
  dist_out <- rbind(dist_out, df_i)
}
# ------------------------------------------------------------------------------



# occupancy analysis ----------------------------------------------------------
# just species with large sample more than 20
d2 <- d1 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
d2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(d2$species))
spp_nams <- unique(d2$species)

# make counts
d3 <- d2 %>%
  mutate(observer_id=str_trim(observer_id)) %>% 
  group_by(site, point, month, species, observer_id) %>% 
  count() 

# zero fill
all_combos <- expand_grid(site=unique(d3$site),
                          point=unique(d3$point),
                          month=unique(d3$month),
                          species=unique(d3$species),
                          observer_id=unique(d3$observer_id))
d4 <- all_combos %>% left_join(d3) %>% 
  mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(site, point, month, species) %>% 
  pivot_wider(values_from = n, names_from=observer_id) %>% 
  rename(y1=SANTIAGO, y2=JUAN) %>% 
  arrange(species, site, point, month) %>% 
  mutate(y2=ifelse(y2>0, 1, 0),
         y1=ifelse(y1>0, 1, 0))

# loop through species
occ_out <- c()
for(i in 1:length(spp_nams)){
  d_i <- d4 %>% filter(species==spp_nams[i]) %>% select(y1, y2) %>% as.matrix()
  umf_i <- unmarkedFrameOccu(y=d_i)
  om_i <- occu(~1 ~1, umf_i)
  p_i <- as.numeric(plogis(coef(om_i)[2]))
  df_i <- data.frame(species=spp_nams[i],
                      occ_det_p=p_i)
  occ_out <- rbind(occ_out, df_i)
}
# ------------------------------------------------------------------------------



# combine detection probabilities --------------------------------------------------------
detect_scores <- spp1 %>% select(species, com_name=Clements_2023_com_name,
                                detect_score=detectability_score) %>% 
  left_join(dist_out) %>% left_join(occ_out)
(dist_mean_coefs <- coef(lm(dist_det_p~-1+factor(detect_score), 
                            data=detect_scores)))
(occ_mean_coefs <- coef(lm(occ_det_p~-1+factor(detect_score), 
                           data=detect_scores)))
detect_scores <- detect_scores %>% 
  mutate(dist_det_p=ifelse(is.na(dist_det_p) & detect_score==1, 
                           dist_mean_coefs[1], 
                           ifelse(is.na(dist_det_p) & detect_score==2,
                                  dist_mean_coefs[2], 
                                  ifelse(is.na(dist_det_p) & 
                                           detect_score==3,
                                         dist_mean_coefs[3], dist_det_p)))) %>% 
  mutate(occ_det_p=ifelse(is.na(occ_det_p) & detect_score==1, 
                           occ_mean_coefs[1], 
                           ifelse(is.na(occ_det_p) & detect_score==2,
                                  occ_mean_coefs[2], 
                                  ifelse(is.na(occ_det_p) & 
                                           detect_score==3,
                                         occ_mean_coefs[3], occ_det_p))))
#-------------------------------------------------------------------------------



# define functional species ----------------------------------------------------
library(cluster)
library(NbClust)
head(spp1)

# functional data for clustering
t1 <- spp1 %>% filter(guild=="landbird") %>% 
  rename_with(~tolower(gsub(".", "_", .x, fixed=T))) %>% 
  select(-(2:7))

# make a distance matrix
d_gower <- daisy(t1[,-1], metric="gower", stand=TRUE,
                 weights=c(rep(3/9/10, 10), 
                           rep(3/9/8, 8),
                           rep(3/9/1, 1)))
# create clusters
clust <- NbClust(diss=d_gower, distance=NULL, method="ward.D", 
                 index="silhouette")

# best clusters
func_out <-  spp1 %>% filter(guild=="landbird") %>% 
  select(species, com_name=Clements_2023_com_name) %>% 
  mutate(func_spec=clust$Best.partition) %>% 
  arrange(func_spec)
# write.csv(inddf, "functional species v1.csv")
# ------------------------------------------------------------------------------





# get conservation scores ------------------------------------------------------
cscore_out <- spp1 %>% filter(guild=="landbird") %>% 
  rename_with(~tolower(gsub(".", "_", .x, fixed=T))) %>% 
  select(species, ccs_b)
# ------------------------------------------------------------------------------




# combine all species level info -----------------------------------------------
final_spp_tab <- detect_scores %>% 
  left_join(func_out %>% select(-2)) %>% 
  left_join(cscore_out)
write_csv(final_spp_tab, "final_species_table_030824.csv")
# ------------------------------------------------------------------------------



# combine counts and species level data ----------------------------------------
head(d1)
unique(d1$species)
d5 <- d1 %>% group_by(watershed, site, point, month, observer_id, species) %>% 
  count() %>% 
  ungroup %>% 
  group_by(watershed, site, point, month, species) %>% 
  summarise(mean_count=mean(n)) %>% 
  ungroup() %>% 
  left_join(final_spp_tab) %>% 
  mutate(occ_correct_mean_count=mean_count*(1/occ_det_p)) %>% 
  mutate(detect_correct_mean_density=(mean_count*
                                   (1/dist_det_p))*(pi*max_d*max_d/(100*100))) %>% 
  mutate(score_weighted_dens=detect_correct_mean_density*ccs_b)
# ------------------------------------------------------------------------------




# point level bfi and other metrics --------------------------------------------
point_bfi <- d5 %>% group_by(watershed, site, point, month) %>% 
  summarize(spp_rich=n(), 
            func_rich=length(unique(func_spec)),
            func_shannon=diversity(detect_correct_mean_density),
            sum_cw_abund=sum(score_weighted_dens)) %>% 
  ungroup() %>% 
  mutate(bfi_unscaled=func_shannon*sum_cw_abund,
         bfi_scaled=plogis(as.numeric(scale(bfi_unscaled))))
# ------------------------------------------------------------------------------
  


# join bfi with habitat --------------------------------------------------------
h1 <- read_csv("Caracterización_puntos_de_muestreo.csv") %>% 
  mutate(site_point=str_trim(site_point),
         site=str_split(site_point, pattern="-", simplify = T)[,1],
         point=str_split(site_point, pattern="-", simplify = T)[,2]) %>% 
  select(-site_point) %>% 
  select(watershed, site, point, x, y, starts_with("perc")) %>% 
  arrange(site, point)
point_bfi_habitat <- point_bfi %>% left_join(h1) %>% 
  mutate(watershed=factor(watershed),
         site=factor(site), 
         point=factor(point))
# ------------------------------------------------------------------------------





# quick model ------------------------------------------------------------------
names(point_bfi_habitat)
mod1 <- gam(bfi_scaled ~ 1 + site +
              s(perc_forest_no_fence, k=3) + 
              s(perc_savanna, k=3) +
              s(perc_pasture, k=3) +
              s(perc_water, k=3) + 
              s(month, k=3),
            data=point_bfi_habitat)
summary(mod1)
b <- getViz(mod1)
print(plot(b, allTerms = T), pages = 1)
# ------------------------------------------------------------------------------


