set.seed(210)
J.x <- 10
J.y <- 10 
J <- J.x * J.y
# Number of distance bins from which to simulate data. 
n.bins <- 5
# Length of each bin. This should be of length n.bins
bin.width <- c(.10, .10, .20, .3, .1)
# Number of species
n.sp <- 10
# Community-level abundance coefficients
beta.mean <- c(-1, 0, 0, 0)
p.abund <- length(beta.mean)
tau.sq.beta <- c(0.2, 0, 0, 0)
# Detection coefficients
alpha.mean <- c(-2.0, 0)
p.det <- length(alpha.mean)
tau.sq.alpha <- c(0.1, 0)
# Detection decay function
det.func <- 'halfnormal'
mu.RE <- list()
p.RE <- list()
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = n.sp, ncol = p.abund)
alpha <- matrix(NA, nrow = n.sp, ncol = p.det)
for (i in 1:p.abund) {
  beta[, i] <- rnorm(n.sp, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(n.sp, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
sp <- FALSE 
family <- 'Poisson'
kappa <- runif(n.sp, 0.3, 3) 
offset <- pi * .8^2
transect <- 'point'
factor.model <- FALSE

dat <- simMsDS(J.x = J.x, J.y = J.y, n.bins = n.bins, bin.width = bin.width,
               n.sp = n.sp, beta = beta, alpha = alpha, det.func = det.func, 
               mu.RE = mu.RE, p.RE = p.RE, sp = sp, 
               sigma.sq = sigma.sq, phi = phi, nu = nu, family = family, 
               offset = offset, transect = transect, factor.model = factor.model)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
dist.breaks <- dat$dist.breaks

covs <- cbind(X, X.p)
colnames(covs) <- c('int.abund', 'abund.cov.1', 'abund.cov.2', 'abund.cov.3', 
                    'int.det', 'det.cov.1')

data.list <- list(y = y, 
                  covs = covs,
                  dist.breaks = dist.breaks, 
                  offset = offset)

# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
                   alpha.comm.normal = list(mean = 0,
                                            var = 10), 
                   kappa.unif = list(0, 100), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1)) 
# Starting values
inits.list <- list(alpha.comm = 0, beta.comm = 0, beta = 0,
                   alpha = 0, kappa = 1)

tuning <- list(beta = 0.1, alpha = 0.1, beta.star = 0.3, alpha.star = 0.1, 
               kappa = 0.8) 

n.batch <- 10 
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.chains <- 1

out <- msDS(abund.formula = ~ 1,
            det.formula = ~ 1,
            data = data.list, 
            n.batch = n.batch, 
            batch.length = batch.length, 
            inits = inits.list, 
            family = 'Poisson',
            det.func = 'halfnormal', 
            transect = transect, 
            tuning = tuning,
            priors = prior.list, 
            accept.rate = 0.43, 
            n.omp.threads = 1, 
            verbose = TRUE, 
            n.report = 10,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains) 
summary(out, level = 'community')
