# Simulate some line-transect data

R <- 15 # number of transects
T <- 2  # number of replicates
strip.width <- 50
transect.length <- 100
breaks <- seq(0, 50, by=10)
lambda <- 1 # Abundance
phi <- 0.20  # Availability
sigma <- 20 # Half-normal shape parameter

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
  M <- rpois(1, lambda) # Individuals within the 1-ha strip
  for(t in 1:T) {
    # Distances from point
    d <- runif(M, 0, strip.width)
    # Detection process
    if(length(d)) {
      cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
      d <- d[rbinom(length(d), 1, cp) == 1]
      y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
    }
  }
}
y <- matrix(y, nrow=R) # convert array to matrix

# Organize data
umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                        dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)

# Fit the model
m1 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=50)
summary(m1)
backTransform(m1, type="lambda")
backTransform(m1, type="phi")
backTransform(m1, type="det")



