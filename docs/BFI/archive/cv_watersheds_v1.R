# set up -----------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(stringr)
setwd("C:/Users/tmeehan/Box/Cauca Valley monitoring/r_work")
dir()

# get raw count data
d0 <- read.csv("BASE DE DATOS FINAL AUDUBON.csv") %>% 
  rename_all(tolower) %>% 
  rename_all(~ gsub(".", "_", ., fixed=T)) %>% 
  filter(tipo_de_registro %in% c("Visual","Auditivo","visual","auditivo")) %>%
  mutate(distance=as.numeric(stringr::str_trim(distancia))) %>% 
  mutate(date=dmy(paste0(fecha, "-2023"))) %>% 
  mutate(hora_inicio=hm(hora_inicio),
         hora_final=hm(hora_final)) %>% 
  mutate(start_time=hour(hora_inicio)+(minute(hora_inicio)/60),
         end_time=hour(hora_final)+(minute(hora_final)/60),
         month=month(date)) %>% 
  mutate(count=as.numeric(individuos)) %>% 
  filter(!is.na(count)) %>% 
  select(watershed=cuenca, site=celda, point=punto, date, month, 
         observer_id=observador, start_time, end_time, species=especie, 
         distance, count)
head(d0)

# expand data table for distance sampling
names(d0)
dis1 <- d0 %>% uncount(count) %>% 
  filter(!is.na(distance)) 

# look at sample sizes
dis1 %>% group_by(species) %>% count() %>% ungroup() %>% arrange(desc(n)) %>% 
  View()

# truncate data
hist(dis1$distance)
max_d <- round(ceiling(quantile(dis1$distance, probs=0.95)), -1); max_d
max_d <- 80

# truncated data
dis2 <- dis1 %>% filter(distance<=max_d) %>% 
  mutate(species_idx=as.integer(factor(species)))

# offset
pi * max_d^2 / 10000

# get species list for distance analysis
# spp_lst <- d0 %>% group_by(especie) %>% summarise(tot_n=sum(individuos)) %>% 
#   filter(tot_n>=30) %>% pull(especie)
# ------------------------------------------------------------------------------



# simple Distance model --------------------------------------------------------
library(Distance)
head(dis2)
dis3 <- dis2 %>% group_by(species) %>% mutate(n_detects=n()) %>% 
  filter(n_detects>10) %>% select(-n_detects, -species_idx) %>% 
  ungroup()
dm1 <- ds(data=dis3,
          transect = "point", key="hn", 
          formula= ~ 1,
          adjustment = NULL)
summary(dm1)
plot(dm1)

# effective radius
sig <- exp(summary(dm1)$ds$coeff$key.scale$estimate)
sig_lb <- exp(summary(dm1)$ds$coeff$key.scale$estimate - 
                (1.97 * summary(dm1)$ds$coeff$key.scale$se))
sig_ub <- exp(summary(dm1)$ds$coeff$key.scale$estimate + 
                (1.97 * summary(dm1)$ds$coeff$key.scale$se))
ea <- 2*pi * integrate(grhn, 0, max_d, sigma=sig)$value # effective area
sqrt(ea / pi) # effective radius

# point count detection probability
ea / (pi*max_d^2)

# strange mean prob
x_vals <- 0:max_d
mean(gxhn(x_vals, sig))

# function plot
data.frame(med = gxhn(x_vals, sig), 
                        low = gxhn(x_vals, sig_lb), 
                        high = gxhn(x_vals, sig_ub),
                        x.val = x_vals) %>% 
  ggplot() + 
  geom_ribbon(aes(x = x.val, ymin = low, ymax = high), fill = 'grey', 
              alpha = 0.7) +
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val, y = med), col = 'black', linewidth = 1.3) + 
  labs(x = 'Distance (m)', y = 'Detection Probability')
# ------------------------------------------------------------------------------



# unmarked version -------------------------------------------------------------
library(unmarked)

# make y array
names(dis2)
dis4 <- dis2 %>% mutate(distance_bin=ifelse(distance<20, 10, ifelse(
  distance>=20 & distance<40, 30, ifelse(
    distance>=40 & distance<60, 50, 70)))) 
hist(dis4$distance_bin)
summary(dis4)
dis4 <- dis4 %>% group_by(species, distance_bin) %>% count() 
all_combos <- expand_grid(species=unique(dis4$species),
                          distance_bin=unique(dis4$distance_bin))
dis5 <- all_combos %>% left_join(dis4) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(distance_bin, species, n)
dis6 <- dis5 %>% pivot_wider(values_from=n, names_from=distance_bin) %>% 
  rename_if(is.numeric, ~paste0("dc", .))

# make a umf frame
umf1 <- with(dis6, {
  unmarkedFrameDS(y = cbind(dc10, dc30, dc50, dc70),
                  siteCovs = data.frame(species),
                  dist.breaks = c(0, 20, 40, 60, 80),
                  survey = "point", unitsIn = "m")
})
summary(umf1)

# make model
dm2 <- distsamp(~ 1 ~ 1, umf1)
summary(dm2)
predict(dm2, type="det")
hist(dm2, ylim=c(0, 0.03), xlab="Distance (m)")

# effective radius
sig <- exp(coef(dm2, type="det"))
ea <- 2*pi * integrate(grhn, 0, 80, sigma=sig[1])$value # effective area
sqrt(ea / pi) # effective radius

# detection probability
ea / (pi*max_d^2)
# ------------------------------------------------------------------------------



# inla distance ----------------------------------------------------------------
library(inlabru)
library(INLA)

# distance function
hn <- function(distance, sigma) {
  exp(-0.5*(distance / sigma)^2)
}
sigma <- function(x) {
  bru_forward_transformation(qexp, x, rate = 1 / (max_d/2))
}

# linear predictor
names(dis2)
formula2 <- distance ~ log(hn(distance, sigma(sigma_theta))) + Intercept

# components
cmp2 <- ~ sigma_theta(1, prec.linear = 1) + Intercept(1)
n_species <- length(unique(dis2$species))

# model object
in1 <- lgcp(
  components = cmp2,
  dis2,
  domain = list(distance = fm_mesh_1d(seq(0, max_d, length.out = 30))),
  formula = formula2,
  options = list(bru_initial = list(sigma_theta = 1, Intercept = 3))
)
summary(in1)

# model predictions
distdf <- expand.grid(distance = seq(0, max_d, length.out = 100))
detfun <- predict(in1, distdf, ~ hn(distance, sigma(sigma_theta)))
plot(detfun)
summary(detfun)
# ------------------------------------------------------------------------------




# inla distance ----------------------------------------------------------------
library(inlabru)
library(INLA)

# distance function
hn <- function(distance, sigma) {
  exp(-0.5 * (distance / sigma)^2)
}
sigma <- function(x) {
  bru_forward_transformation(qexp, x, rate = 1 / (max_d/2))
}


names(dis2)
n_species <- length(unique(dis2$species_idx))
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


# model predictions
distdf <- expand.grid(distance = seq(0, max_d, length.out = 100),
                      species_idx=sort(unique(dis2$species_idx)))
detfun <- predict(dfit1, distdf, ~ hn(distance, sigma))

# plot and summarize probabilities
ggplot() + geom_line(data=detfun, 
                     aes(x=distance, y=mean, group=species_idx, 
                         color=(species_idx))) +
  scale_color_distiller(palette="Spectral")
# ------------------------------------------------------------------------------




# inla distance ----------------------------------------------------------------
library(inlabru)
library(INLA)

# distance function
hn <- function(distance, sigma) {
  exp(-0.5 * (distance / sigma)^2)
}
sigma <- function(x) {
  bru_forward_transformation(qexp, x, rate = 1 / (max_d))
}

# linear predictor
names(dis2)
formula2 <- distance + species_idx ~ 
  log(hn(distance, sigma_common*exp(sigma_log_factor))) + Intercept

# components
cmp2 <- ~ 
  sigma_common(1, prec.linear = 1, marginal =
                 bru_mapper_marginal(qexp, rate = 1 / (max_d))) +
  sigma_log_factor(species_idx,
                   model = "iid",
                   hyper=list(prec=list(initial=0,
                                        fixed = FALSE)),
                   constr = TRUE) +
  Intercept(species_idx, model = "iid",
            hyper=list(prec=list(initial = -2, fixed = TRUE)))
n_species <- length(unique(dis2$species))

# model object
dm2 <- lgcp(
  components = cmp2,
  dis2,
  domain = list(
    distance = fm_mesh_1d(seq(0, max_d, length.out = 30)),
    species_idx = sort(unique(dis2$species_idx))),
  formula = formula2,
  options = list(bru_initial = list(sigma_common = 1,
                                    sigma_log_factor = rep(0, n_species),
                                    Intercept = rep(3, n_species))))
summary(dm2)
plot(dm2, "sigma_log_factor")
dm2$summary.fixed
dm2$summary.random
dm2$summary.hyperpar

# model predictions
distdf <- expand.grid(distance = seq(0, max_d, length.out = 100),
                      species_idx=sort(unique(dis2$species_idx)))
detfun <- predict(dm2, distdf, ~hn(distance,
                                     sigma_common*exp(sigma_log_factor)))

# plot and summarize probabilities
ggplot() + geom_line(data=detfun, 
                     aes(x=distance, y=mean, group=species_idx, 
                         color=(species_idx)), show.legend = FALSE) +
  scale_color_distiller(palette="Spectral")
# ------------------------------------------------------------------------------



# spAbundance distance ---------------------------------------------------------
# make y array
names(dis2)
dis7 <- dis2 %>% mutate(distance_bin=ifelse(distance<20, 20, ifelse(
  distance>=20 & distance<40, 40, ifelse(
    distance>=40 & distance<60, 60, 80)))) %>% 
  group_by(site, species, distance_bin) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis7$site),
                          species=unique(dis7$species),
                          distance_bin=unique(dis7$distance_bin))
dis8 <- all_combos %>% left_join(dis7) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(distance_bin, species, n)

y1 <- array(dis8$n, dim=c(length(unique(dis8$species)),
                          length(unique(dis8$site)),
                          length(unique(dis8$distance_bin))))
dim(y1)
dimnames(y1)[[1]] <- sort(unique(dis8$species))
dimnames(y1)[[2]] <- sort(unique(dis8$site))
dimnames(y1)[[3]] <- c("20m", "40m", "60m", "80m")
head(y1)
head(dis8)

# make data list
dat_lst <- list()
dat_lst$y <- y1
dat_lst$dist.breaks <- c(0, 0.020, 0.040, 0.060, 0.080)
radius.km <- 0.080
dat_lst$offset <- pi * radius.km^2 * 100
str(dat_lst)

# formulas
abund.formula <- ~ 1
det.formula <- ~ 1

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
                 N = apply(dat_lst$y, c(1, 2), sum, na.rm = TRUE))
# priors
ms.priors <- list(beta.comm.normal = list(mean = 0, var = 100),
                  alpha.comm.normal = list(mean = 0, var = 100),
                  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                  tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                  kappa.unif = list(a = 0, b = 100))

# Specify initial tuning values
ms.tuning <- list(beta = 0.3, alpha = 0.3, beta.star = 0.5, kappa = 0.5)

# Approx. run time:  11.5 min
library(spAbundance)
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
det.comm.quants <- quantile(exp(out.ms$alpha.comm.samples[, 1]), c(0.025, 0.975))
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
ggplot(data = comm.plot.df) +
  geom_ribbon(aes(x = x.val, ymin = low, ymax = high), fill = 'grey',
              alpha = 0.5) +
  geom_line(data = p.plot.df, aes(x = x.val, y = val, col = sp), lwd = 1, lty = 1, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val, y = mean), col = 'black', lwd = 1.3) +
  labs(x = 'Distance (m)', y = 'Detection Probability', col = 'Species')

# convert
sig <- det.means
dp <- c()
for(i in 1:length(sig)){
  ea_i <- 2*pi * integrate(grhn, 0, maxd_km, sigma=sig[i])$value # effective area
  er_i <- sqrt(ea_i / pi) # effective radius
  dp_i <- data.frame(species=sp.names[i], detect_prop=ea_i / (pi*maxd_km^2)) # detection probability
  dp <- rbind(dp, dp_i)
}
View(dp)
hist(dp$detect_prop)
# ------------------------------------------------------------------------------




# spAbundance nmix ------------------------------------------------------------------
# make y array
names(dis2)
dis7 <- dis2 %>% mutate(distance_bin=ifelse(distance<20, 20, ifelse(
  distance>=20 & distance<40, 40, ifelse(
    distance>=40 & distance<60, 60, 80)))) %>% 
  mutate(observer_id=str_trim(observer_id)) %>% 
  group_by(site, species, observer_id) %>% 
  count() 
all_combos <- expand_grid(site=unique(dis7$site),
                          species=unique(dis7$species),
                          observer_id=unique(dis7$observer_id))
dis8 <- all_combos %>% left_join(dis7) %>% mutate(n=ifelse(is.na(n), 0, n)) %>% 
  arrange(observer_id, species, n)

y1 <- array(dis8$n, dim=c(length(unique(dis8$species)),
                          length(unique(dis8$site)),
                          length(unique(dis8$observer_id))))
dim(y1)
dimnames(y1)[[1]] <- sort(unique(dis8$species))
dimnames(y1)[[2]] <- sort(unique(dis8$site))
dimnames(y1)[[3]] <- c("Juan", "Santiago")
head(y1)
head(dis8)

# make data list
dat_lst <- list()
dat_lst$y <- y1
str(dat_lst)

# formulas
abund.formula <- ~ 1
det.formula <- ~ 1

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
ms.tuning <- list(beta = 0.3, alpha = 0.3, beta.star = 0.5, kappa = 0.5)

# Approx. run time:  11.5 min
library(spAbundance)
out.ms <- msNMix(abund.formula = abund.formula,
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
summary(out.ms)

# ppc
ppc.ms.out <- ppcAbund(out.ms, fit.stat = 'chi-squared', group = 2)
summary(ppc.ms.out)

# plot of species-specific detection probability
sp.names <- unlist(dimnames(dat_lst$y)[[1]])
det.int.samples <- out.ms$alpha.samples[, 1:n.sp]
det.means <- apply(plogis(det.int.samples), 2, mean)
det.comm.means <- mean(plogis(out.ms$alpha.comm.samples[, 1]))
det.comm.quants <- quantile(plogis(out.ms$alpha.comm.samples[, 1]), c(0.025, 0.975))

# START HERE
comm.plot.df <- data.frame(mean = gxhn(x.vals, det.comm.means),
                           x.val = x.vals,
                           low = gxhn(x.vals, det.comm.quants[1]),
                           high = gxhn(x.vals, det.comm.quants[2]))
ggplot(data = comm.plot.df) +
  geom_ribbon(aes(x = x.val, ymin = low, ymax = high), fill = 'grey',
              alpha = 0.5) +
  geom_line(data = p.plot.df, aes(x = x.val, y = val, col = sp), lwd = 1, lty = 1, show.legend = FALSE) +
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val, y = mean), col = 'black', lwd = 1.3) +
  labs(x = 'Distance (m)', y = 'Detection Probability', col = 'Species')
# ------------------------------------------------------------------------------









# export species table ---------------------------------------------------------
species_table <- d0 %>% group_by(species=especie) %>% summarise(tot_n=sum(individuos)) %>%
  arrange(desc(tot_n))
write_csv(species_table, "species_table.csv")

# grab species traits
# diet traits data
traits1 <- read.csv("SAviTraits_1-0_1.csv") %>% 
  rename_all(tolower) %>% 
  filter(species_scientific_name %in% species_table$species) %>% 
  select(1:3, 9:14) %>% 
  mutate(diet_cat=tolower(paste(diet_cat, diet_sub_cat, sep="-"))) %>% 
  select(-diet_sub_cat) %>% 
  rowwise() %>% 
  mutate(mean_percent=mean(c(jun, jul, aug, sep, oct, nov))) %>% 
  select(1,2,9) %>% 
  pivot_wider(names_from=diet_cat, values_from=mean_percent)

species_table %>% distinct(species) # 257
traits1 %>% distinct(Species_Scientific_Name) # 246

traits1 %>% group_by(Species_Scientific_Name)

# figure months in dataset
d0 %>% distinct(month(fecha))

# pairs
d0 %>% arrange(fecha, celda, punto) %>% View()

# traits data
traits1 <- read.csv("SAviTraits_1-0_1.csv")

# distance
library(Distance)
spp_lst <- d0 %>% group_by(especie) %>% summarise(tot_n=sum(individuos)) %>% 
  filter(tot_n>=30) %>% pull(especie)
dis1 <- d0 %>% filter(!is.na(distancia)) %>% 
  filter(especie %in% spp_lst) %>% 
  mutate(distance=distancia,
         species=factor(especie))
head(dis1)
dm1 <- ds(data=dis1, transect = "point", truncation ="5%", key="hr", formula = ~ hora_inicio + I(hora_inicio^2) + species)
summary(dm1)
plot(dm1)




# -----------------------------------------------------------------------------