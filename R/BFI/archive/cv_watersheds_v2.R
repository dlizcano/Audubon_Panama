# set up -----------------------------------------------------------------------
library(unmarked)
library(spAbundance)
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

# filter data for distance sampling
names(d0)
dis1 <- d0 %>% filter(!is.na(distance))

# just species with trait data
spp1 <- read.csv("species_detect_diet_acad_forage_trait_data_022024.csv")
dis1 <- dis1 %>% filter(species %in% spp1$species)

# truncate the data to distance
hist(dis1$distance)
max_d <- round(ceiling(quantile(dis1$distance, probs=0.99)), -1); max_d
dis1 <- dis1 %>% filter(distance<=max_d)
hist(dis1$distance)

# just species with large sample more than 25
thresh_ss <- 25
dis2 <- dis1 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
dis2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(dis2$species))
spp_nams <- unique(dis2$species)
radius.km <- max_d / 1000
offset_km2_2_ha <- pi * radius.km^2 * 100; offset_km2_2_ha
# ------------------------------------------------------------------------------











# spAbundance distance ---------------------------------------------------------
# make y array for distance modeling
names(dis2)
dis3 <- dis2 %>% 
  mutate(distance_bin=ifelse(distance<20, .020, ifelse(
  distance>=20 & distance<40, .040, ifelse(
    distance>=40 & distance<60, .060, .080)))) %>% 
  group_by(site, species, distance_bin) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis3$site),
                          species=unique(dis3$species),
                          distance_bin=unique(dis3$distance_bin))
dis4 <- all_combos %>% left_join(dis3) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(distance_bin, species, n)
y1 <- array(dis4$n, dim=c(length(unique(dis4$species)),
                          length(unique(dis4$site)),
                          length(unique(dis4$distance_bin))))
dim(y1)
dimnames(y1)[[1]] <- sort(unique(dis4$species))
dimnames(y1)[[2]] <- sort(unique(dis4$site))
dimnames(y1)[[3]] <- c("0.02km", "0.04km", "0.06km", "0.08km")
head(y1)
head(dis4)

# make data list
dat_lst <- list()
dat_lst$y <- y1
dat_lst$dist.breaks <- c(0, 0.020, 0.040, 0.060, 0.080)
dat_lst$offset <-offset_km2_2_ha
str(dat_lst)

# formulas
abund.formula <- ~ 1
det.formula <- ~ 1

# Number of species and inits
n.sp <- dim(dat_lst$y)[1]
ms.inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 kappa = 1,
                 tau.sq.alpha = 1,
                 sigma.sq.mu = 0.5,
                 N = apply(dat_lst$y, c(1, 2), sum, na.rm = TRUE))

# priors
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 100),
                  alpha.comm.normal = list(mean = 0, var = 100),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                  kappa.unif = list(a = 0, b = 100))

# specify initial tuning values
ms.tuning <- list(beta = 0.3, alpha = 0.3, beta.star = 0.5, kappa = 0.5)

# run model
out.ms <- msDS(abund.formula = abund.formula,
               det.formula = det.formula,
               data = dat_lst,
               inits = ms.inits,
               n.batch = 200,
               tuning = ms.tuning,
               batch.length = 25,
               priors = ms.priors,
               n.omp.threads = 1,
               family = 'NB',
               det.func = 'halfnormal',
               transect = 'point',
               verbose = TRUE,
               n.report = 500,
               n.burn = 2000,
               n.thin = 30,
               n.chains = 3)
summary(out.ms)

# plot of species-specific detection probability
maxd_km <- max_d / 1000
sp.names <- unlist(dimnames(dat_lst$y)[[1]])
det.int.samples <- out.ms$alpha.samples[, 1:n.sp]
det.means <- apply(exp(det.int.samples), 2, mean)
det.comm.means <- mean(exp(out.ms$alpha.comm.samples[, 1]))
det.comm.quants <- quantile(exp(out.ms$alpha.comm.samples[, 1]), 
                            c(0.025, 0.975))
x.vals <- seq(0, maxd_km, length.out = 200)
n.vals <- length(x.vals)
p.plot.df <- data.frame(val = NA,
                        x.val = rep(x.vals, n.sp),
                        sp = rep(sp.names, each = n.vals))
for (i in 1:n.sp) {
  indx <- ((i - 1) * n.vals + 1):(i * n.vals)
  p.plot.df$val[indx] <- gxhn(x.vals, det.means[i])
}
comm.plot.df <- data.frame(mean = gxhn(x.vals, det.comm.means),
                           x.val = x.vals,
                           low = gxhn(x.vals, det.comm.quants[1]),
                           high = gxhn(x.vals, det.comm.quants[2]))
dist_plot <- ggplot(data = comm.plot.df) +
  geom_ribbon(aes(x = x.val*1000, ymin = low, ymax = high), fill = 'grey',
              alpha = 0.5) +
  geom_line(data = p.plot.df, aes(x = x.val*1000, y = val, col = sp), lwd = 0.8, 
            lty = 1, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  geom_line(aes(x = x.val*1000, y = mean), col = 'black', lwd = 1.3) +
  labs(x = 'Distance (m)', y = 'Detection Probability', col = 'Species')
dist_plot

# convert to integrated probabilities
sig <- det.means
dp <- c()
for(i in 1:length(sig)){
  ea_i <- 2*pi * integrate(grhn, 0, maxd_km, sigma=sig[i])$value # effective area in km
  er_i <- sqrt(ea_i / pi) # effective radius
  dp_i <- data.frame(species=sp.names[i], 
                     detect_prop=ea_i / (pi*maxd_km^2)) # detection probability
  dp <- rbind(dp, dp_i)
}
dist_out_detect_df <- dp
head(dist_out_detect_df)
dist_hist <- ggplot(dist_out_detect_df) + 
  geom_histogram(aes(x=detect_prop), col="white") +
  theme_bw(base_size = 12) +
  labs(x = 'Detection probability', y = 'Number of species')
dist_hist
# still need to figure out the units associated with these probabilities

dist_out_detect_df %>% arrange(detect_prop)
# ------------------------------------------------------------------------------




# spAbundance nmix ------------------------------------------------------------------
# make y array for n-mixture models
summary(dis2)
dis5 <- dis2 %>% mutate(distance_bin=ifelse(distance<20, .020, ifelse(
  distance>=20 & distance<40, .040, ifelse(
    distance>=40 & distance<60, .060, .080)))) %>% 
  mutate(observer_id=str_trim(observer_id)) %>% 
  group_by(site, date, species, observer_id) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis5$site),
                          date=unique(dis5$date),
                          species=unique(dis5$species),
                          observer_id=unique(dis5$observer_id))
dis6 <- all_combos %>% left_join(dis5) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  mutate(site_date=paste(site, date, sep="-")) %>% 
  arrange(site, date, species)
y1 <- array(dis6$n, dim=c(length(unique(dis6$species)),
                          length(unique(dis6$site_date)),
                          length(unique(dis6$observer_id))))
dim(y1)
dimnames(y1)[[1]] <- sort(unique(dis6$species))
dimnames(y1)[[2]] <- sort(unique(dis6$site_date))
dimnames(y1)[[3]] <- c("Juan", "Santiago")
head(y1)
head(dis6)

# make data list
abund.covs <- data.frame(abund.factor.1=as.numeric(factor(unique(dis6$site_date))))
det.covs <- data.frame(det.factor.1=as.numeric(factor(unique(dis6$site_date))))
dat_lst <- list()
dat_lst$y <- y1
dat_lst$abund.covs <- abund.covs
dat_lst$det.covs <- det.covs
str(dat_lst)

# formulas
abund.formula <- ~ 1 + (1|abund.factor.1)
det.formula <- ~ 1 + (1|det.factor.1)

# Number of species
n.sp <- dim(dat_lst$y)[1]
ms.inits <- list(alpha.comm = 0,
                 beta.comm = 0,
                 beta = 0,
                 alpha = 0,
                 tau.sq.beta = 1,
                 kappa = 1,
                 tau.sq.alpha = 1,
                 sigma.sq.mu = 0.5,
                 N = apply(dat_lst$y, c(1, 2), max, na.rm = TRUE))
# priors
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 100),
                  alpha.comm.normal = list(mean = 0, var = 2.72),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                  sigma.sq.mu.ig = list(a = 0.1, b = 0.1),
                  kappa.unif = list(a = 0, b = 100))

# Specify initial tuning values
ms.tuning <- list(beta = 0.3, alpha = 0.3, beta.star = 0.5, 
                  kappa = 0.5, alpha.star=0.3)

# Approx. run time:  11.5 min
out.ms.2 <- msNMix(abund.formula = abund.formula,
                 det.formula = det.formula,
                 data = dat_lst,
                 inits = ms.inits,
                 n.batch = 1000,
                 tuning = ms.tuning,
                 batch.length = 25,
                 priors = ms.priors,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 family = 'NB',
                 n.report = 400,
                 n.burn = 2000,
                 n.thin = 20,
                 n.chains = 3)
summary(out.ms.2)
str(out.ms.2)

# plot of species-specific detection probability
sp.names <- str_trim(unlist(dimnames(dat_lst$y)[[1]]))
det.int.samples <- out.ms.2$alpha.samples[, 1:n.sp]
det.means <- apply(plogis(det.int.samples), 2, mean)
det.quants <- apply(plogis(det.int.samples), 2, quantile, 
                    probs=c(0.025, 0.975)) %>% t() %>% 
  as.data.frame() %>% 
  mutate(species=str_trim(gsub("(Intercept)-", "", row.names(.), fixed=T))) %>% 
  remove_rownames()
det.comm.means <- mean(plogis(out.ms.2$alpha.comm.samples[, 1]))
det.comm.quants <- quantile(plogis(out.ms.2$alpha.comm.samples[, 1]), c(0.025, 0.975))
nmix_out_detect_df <- data.frame(species=sp.names,
                             detect_p=det.means) %>% 
  remove_rownames() %>% 
  left_join(det.quants)
head(nmix_out_detect_df)                           
 
nmix_hist <- ggplot(nmix_out_detect_df) + 
  geom_histogram(aes(x=detect_p), col="white") +
  theme_bw(base_size = 12) +
  labs(x = 'Detection probability', y = 'Number of species')
nmix_hist
cowplot::plot_grid(dist_plot, nmix_hist)
# ------------------------------------------------------------------------------


# combine probabilities --------------------------------------------------------
detect_probs <- dist_out_detect_df %>% rename(detect_p_dist=detect_prop) %>% 
  left_join(nmix_out_detect_df %>% select(species, detect_p_nmix=detect_p)) %>% 
  mutate(detect_p_comb=detect_p_dist * detect_p_nmix)
#-------------------------------------------------------------------------------


# export species table ---------------------------------------------------------
species_table <- d0 %>% group_by(species) %>% 
  summarise(n_counted=sum(count)) %>%
  arrange(desc(n_counted))
species_table <- detect_probs %>% left_join(species_table) %>%
  arrange(desc(n_counted))
write_csv(species_table, "species_table.csv")
# ------------------------------------------------------------------------------




# grab species traits ----------------------------------------------------------
# reopen species table
species_table <- read.csv("species_table.csv")
species_table %>% distinct(species) # 242

# diet traits data
traits1 <- read.csv("SAviTraits_1-0_1.csv") %>% 
  rename_all(tolower) %>% 
  rename_all(~ gsub(".", "_", ., fixed=T)) %>% 
  filter(species_scientific_name %in% species_table$species) %>% 
  select(1:3, 9:14) %>% 
  mutate(diet_cat=tolower(paste(diet_cat, diet_sub_cat, sep="-"))) %>% 
  select(-diet_sub_cat) %>% 
  rowwise() %>% 
  mutate(mean_percent=mean(c(jun, jul, aug, sep, oct, nov))) %>% 
  select(1,2,9) %>% 
  pivot_wider(names_from=diet_cat, values_from=mean_percent)
traits1 %>% distinct(species_scientific_name) # 232

# merge diets with species table
names(species_table)
names(traits1)
species_table <- species_table %>% left_join(traits1 %>% rename(species=species_scientific_name))
names(species_table)

# acad
acad1 <- read.csv("ACAD Global 2023.csv") %>% 
  select(species=Scientific.Name, com_name=Common.Name, guild=Group,
         ccs_b=CCS.b) %>% 
  rename_all(tolower) %>% 
  rename_all(~ gsub(".", "_", ., fixed=T)) 
names(acad1)

# merge acad and species table
species_table <- species_table %>% left_join(acad1)
summary(species_table)

# elton traits
et1 <- read.csv("EltonTraits_BirdFuncDat.csv") %>% 
  select(species=Scientific, com_name2=English, 10:20, 24:31, 
         body_mass=BodyMass.Value) %>% 
  rename_all(tolower) %>% 
  rename_all(~ gsub(".", "_", ., fixed=T)) 
names(et1)

# merge species table and elton table
species_table <- species_table %>% left_join(et1)

# export species diet, acad, trait data
write.csv(species_table, "species_detect_diet_acad_forage_trait_data_011124.csv", na = "", row.names=F)
# ------------------------------------------------------------------------------




# inla distance ----------------------------------------------------------------
library(inlabru)
library(INLA)

# n species
n_species <- length(unique(dis2$species_idx))

# distance function
hn <- function(distance, sigma) {
  exp(-0.5 * (distance / sigma)^2)
}
sigma <- function(x) {
  bru_forward_transformation(qexp, x, rate = 1 / (max_d/2))
}

# linear predictor
formula1 <- distance + species_idx ~ log(hn(distance, sigma)) + Intercept

# components
cmp1 <- ~
  sigma(species_idx,
        model = "iid",
        hyper=list(prec=list(initial=0, fixed = TRUE)),
        marginal = bru_mapper_marginal(qexp, rate = 1 / (max_d*2))) +
  Intercept(species_idx, model = "iid",
            hyper=list(prec=list(initial = -2, fixed = TRUE)))

# model object
dfit1 <- lgcp(
  components = cmp1,
  dis2,
  domain = list(
    distance = fm_mesh_1d(seq(0, max_d, length.out = 30)),
    species_idx = sort(unique(dis2$species_idx))),
  formula = formula1,
  options = list(bru_initial = list(sigma = rep(1, n_species),
                                    Intercept = rep(3, n_species))))
summary(dfit1)
sort((dfit1$summary.random$Intercept$mean))
sort((dfit1$summary.random$sigma$mean))

# model predictions
distdf <- expand.grid(distance = seq(0, max_d, length.out = 100),
                      species_idx=sort(unique(dis2$species_idx)))
detfun <- predict(dfit1, distdf, ~ hn(distance, sigma))

# plot and summarize probabilities
ggplot() + geom_line(data=detfun, 
                     aes(x=distance, y=mean, group=species_idx, 
                         color=(species_idx))) +
  scale_color_distiller(palette="Spectral")


# plot of species-specific detection probability
maxd_km <- max_d
sp.names <- sort(unique(dis2$species))
det.means <- dfit1$summary.random$sigma$mean
hist((dfit1$summary.random$sigma$mean))
x.vals <- seq(0, maxd_km, length.out = 200)
n.vals <- length(x.vals)
n.sp <- n_species
p.plot.df <- data.frame(val = NA,
                        x.val = rep(x.vals, n.sp),
                        sp = rep(sp.names, each = n.vals))
for (i in 1:n.sp) {
  indx <- ((i - 1) * n.vals + 1):(i * n.vals)
  p.plot.df$val[indx] <- gxhn(x.vals, det.means[i])
}
dist_plot <- ggplot() +
  geom_line(data = p.plot.df, aes(x = x.val, y = val, col = sp), lwd = 0.8, 
            lty = 1, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(x = 'Distance (m)', y = 'Detection Probability', col = 'Species')
dist_plot

# convert to integrated probabilities
sig <- det.means
dp <- c()
for(i in 1:length(sig)){
  ea_i <- 2*pi * integrate(grhn, 0, max_d, sigma=sig[i])$value # effective area in km
  er_i <- sqrt(ea_i / pi) # effective radius
  dp_i <- data.frame(species=sp.names[i], 
                     detect_prop=ea_i / (pi*maxd_km^2)) # detection probability
  dp <- rbind(dp, dp_i)
}
dist_out_detect_df <- dp
head(dist_out_detect_df)
dist_hist <- ggplot(dist_out_detect_df) + 
  geom_histogram(aes(x=detect_prop), col="white") +
  theme_bw(base_size = 12) +
  labs(x = 'Detection probability', y = 'Number of species')
dist_hist
# ------------------------------------------------------------------------------




# inla distance ----------------------------------------------------------------
library(inlabru)
library(INLA)

# data
names(dis2)
dis2$species_idx <- as.numeric(factor(dis2$species))

# distance function
hn <- function(distance, sigma) {exp(-0.5 * (distance / sigma)^2)}
sigma <- function(x) {bru_forward_transformation(qexp, x, rate = 1 / (max_d))}

# linear predictor
form1 <- distance + species_idx ~ Intercept + 
  log(hn(distance, sigma_common*exp(sigma_log_factor))) + 
  log(1e-6 + 2*pi*distance)

# components
cmp1 <- ~ 
  Intercept(species_idx, model = "iid",
            hyper=list(prec=list(initial = -2, fixed = TRUE))) +
  sigma_common(1, prec.linear = 1, marginal =
                 bru_mapper_marginal(qexp, rate = 1 / (max_d))) +
  sigma_log_factor(species_idx,
                   model = "iid",
                   hyper=list(prec=list(initial=0,
                                        fixed = FALSE)),
                   constr = TRUE)

# model object
mean_abund <- nrow(dis2) / max(dis2$species_idx)
dm1 <- lgcp(
  components = cmp1,
  dis2,
  domain = list(
    distance = fm_mesh_1d(seq(0, max_d, length.out = 50)),
    species_idx = sort(unique(dis2$species_idx))),
  formula = form1,
  options = list(bru_initial = list(sigma_common = -1,
                                    sigma_log_factor = rep(0, n_species),
                                    Intercept = rep(mean_abund, n_species))))

# model predictions
distdf <- expand.grid(distance = seq(1, max_d, length.out = 100),
                      species_idx=sort(unique(dis2$species_idx))) %>% 
  mutate(species=factor(paste0("sp", species_idx)))
detfun <- predict(dm1, distdf, ~hn(distance, sigma_common*exp(sigma_log_factor)))
summary(detfun)
# plot and summarize probabilities from inlabru
inlabru_plot <- ggplot() + geom_line(data=detfun, lwd=1,
                                     aes(x=distance, y=mean, group=species, 
                                         color=species)) +
  theme_bw(base_size = 10) +
  #scale_color_brewer(palette="Dark2") +
  labs(x = 'Distance (m)', y = 'inlabru detection probability', col = 'Species'); 

# get probs
sig1 <- predict(dm1, distdf %>% distinct(species_idx, species), 
                ~sigma_common*exp(sigma_log_factor))$mean
sp_names <- unique(dis2$species)
dp1 <- c()
for(i in 1:length(sig1)){
  ea_i <- 2*pi * integrate(grhn, 0, max_d, sigma=sig1[i])$value # effective area
  er_i <- sqrt(ea_i / pi) # effective radius
  dp_i <- data.frame(species=sp_names[i], 
                     detect_prob=ea_i / (pi*max_d^2)) # detection probability
  dp1 <- rbind(dp1, dp_i)
}

# show it
dp1 %>% mutate(sigma_val=sig1) %>% 
  left_join(spp1 %>% select(species, detectability_score), by="species") %>% 
  group_by(detectability_score) %>% summarise(mean(detect_prob))
summary(dp1)
inlabru_plot
# ------------------------------------------------------------------------------






# now nmix
# libraries
library(INLA)
library(unmarked)

# just species with trait data
spp1 <- read.csv("species_detect_diet_acad_forage_trait_data_022024.csv")
dis1 <- d0 %>% filter(species %in% spp1$species)

# just species with large sample more than 25
thresh_ss <- 25
dis2 <- dis1 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
dis2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(dis2$species))
spp_nams <- unique(dis2$species)

# with co data
summary(dis2)
dis5 <- dis2 %>%
  mutate(observer_id=str_trim(observer_id)) %>% 
  group_by(site, month, species, observer_id) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis5$site),
                          month=unique(dis5$month),
                          species=unique(dis5$species),
                          observer_id=unique(dis5$observer_id))
dis6 <- all_combos %>% left_join(dis5) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  mutate(site_month=paste(site, month, sep="-")) %>% 
  arrange(site, month, species)

dis7 <- dis6 %>% pivot_wider(values_from = n, names_from=observer_id) %>% 
  rename(y1=SANTIAGO, y2=JUAN) %>% 
  group_by(site_month) %>% 
  mutate(detected=ifelse(y1>0|y2>0, 1, 0)) %>% 
  ungroup() %>% 
  # filter(detected>0) %>% 
  mutate(species_idx=as.integer(factor(species)),
         site_month_idx=as.integer(factor(site_month)),
         mean_y = (y1 + y2) / 2) %>% 
  group_by(species) %>% 
  mutate(species_mean=mean(mean_y)) %>% 
  ungroup()

head(dis7)
ymat <- dis7 %>% select(y1, y2) %>% as.matrix()
detected <- dis7$detected
mean_y <- dis7$mean_y
species_mean <-dis7$species_mean
counts.and.count.covs <- inla.mdata(ymat, 1, species_mean)

# inla result, CRASHES
out.inla <- inla(counts.and.count.covs ~ -1 +  
                   f(species_idx, model="iid", constr=F) + 
                   f(site_month_idx, model="iid", constr=T),
                   data = list(counts.and.count.covs = counts.and.count.covs,
                               species_idx = dis7$species_idx,
                               site_month_idx = dis7$site_month_idx),
                   family = "nmixnb",
                   control.inla = list(cmin = 0))
summary(out.inla, digits = 3)
sort(plogis(out.inla$summary.random$species_idx$mean))
# show it
data.frame(species=spp_nams, 
           nmix_prob=plogis(out.inla$summary.random$species_idx$mean)) %>% 
  left_join(spp1 %>% select(species, detectability_score), by="species") %>% 
  group_by(detectability_score) %>% summarise(mean(nmix_prob))














# occupancy models -------------------------------------------------------------
# libraries
library(INLA)
library(unmarked)

# just species with trait data
spp1 <- read.csv("species_detect_diet_acad_forage_trait_data_022024.csv")
dis1 <- d0 %>% filter(species %in% spp1$species)

# just species with large sample more than 25
thresh_ss <- 30
dis2 <- dis1 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
dis2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(dis2$species))
spp_nams <- unique(dis2$species)

# with co data
summary(dis2)
dis5 <- dis2 %>%
  mutate(observer_id=str_trim(observer_id)) %>% 
  group_by(site, point, month, species, observer_id) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis5$site),
                          point=unique(dis5$point),
                          month=unique(dis5$month),
                          species=unique(dis5$species),
                          observer_id=unique(dis5$observer_id))
dis6 <- all_combos %>% left_join(dis5) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(site, point, month, species)

dis7 <- dis6 %>% pivot_wider(values_from = n, names_from=observer_id) %>% 
  rename(y1=SANTIAGO, y2=JUAN) %>% 
  mutate(species_idx=as.integer(factor(species))) %>% 
  arrange(species, site, point, month) %>% 
  mutate(y2=ifelse(y2>0, 1, 0),
         y1=ifelse(y1>0, 1, 0))

occ_out <- c()
for(i in 1:max(dis7$species_idx)){
  y <- dis7 %>% filter(species_idx==i) %>% select(y1, y2) %>% as.matrix()
  umf <- unmarkedFrameOccu(y=y)
  fm <- occu(~1 ~1, umf)
  p <- plogis(coef(fm)[2])
  out_i <- data.frame(species=unique(dis7$species)[i],
                      det_p=p)
  occ_out <- rbind(occ_out, out_i)
}


test1 <- data.frame(species=spp_nams, 
           nmix_prob=plogis(out.inla$summary.random$species_idx$mean)) %>% 
  left_join(occ_out) %>% left_join(spp1 %>% select(species, detectability_score), by="species")
# ------------------------------------------------------------------------------





# distance models --------------------------------------------------------------
# libraries
library(INLA)
library(unmarked)
max_d <- 100

# just species with trait data
spp1 <- read.csv("species_detect_diet_acad_forage_trait_data_022024.csv")
dis1 <- d0 %>% filter(species %in% spp1$species)

# filter out no distance data and distance over 100
dis1 <- dis1 %>% 
  filter(!is.na(distance),
         distance<=max_d)

# just species with large sample more than 30
thresh_ss <- 30
dis2 <- dis1 %>% group_by(species) %>% mutate(n_counted=n()) %>% 
  ungroup() %>% 
  filter(n_counted>=thresh_ss) %>% 
  arrange(species)
dis2 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# info
n_species <- length(unique(dis2$species))
spp_nams <- unique(dis2$species)

# with co data
dis3 <- dis2 %>% 
  mutate(distance_bin=ifelse(distance<max_d/4, max_d/4, ifelse(
    distance>=max_d/4 & distance<max_d/2, max_d/2, ifelse(
      distance>=max_d/2 & distance<max_d*3/4, max_d*3/4, max_d)))) %>% 
  group_by(site, point, month, species, distance_bin) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis3$site),
                          point=unique(dis3$point),
                          month=unique(dis3$month),
                          species=unique(dis3$species),
                          distance_bin=unique(dis3$distance_bin))

dis4 <- all_combos %>% left_join(dis3) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(distance_bin, species, n)



dis5 <- dis4 %>% filter(species=="Tyrannus melancholicus")
dis5 %>% group_by(distance_bin) %>% summarise(sum(n))

dis5 <- dis5 %>% pivot_wider(values_from=n, names_from=distance_bin)

umf <- unmarkedFrameDS(y=as.matrix(dis5[,5:8]),
                #siteCovs=data.frame(),
                dist.breaks=c(0,max_d/4,max_d/2,max_d*3/4, max_d), unitsIn="m", survey="point")
summary(umf)
dm1 <- distsamp(~1 ~1, umf)
summary(dm1)


hist(dm1)

# effective radius
sig <- exp(coef(dm1, type="det"))
ea <- 2*pi * integrate(grhn, 0, max_d, sigma=sig)$value # effective area
sqrt(ea / pi) # effective radius

# detection probability
ea / (pi*max_d^2)
