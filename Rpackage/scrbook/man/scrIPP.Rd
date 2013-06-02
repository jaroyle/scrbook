\name{scrIPP}
\alias{scrIPP}
\title{
  A Metropolis-within-Gibbs sampler to obtain posterior distributions
  for parameters of an SCR model with an inhomogeneous binomial point
  process.
}
\description{
  This is not meant to be generic. It assumes the state-space is rectangular.
}
\usage{
scrIPP(Z, X, M, niters, xlims, ylims, tune = rep(0.1, 4))
}
\arguments{
  \item{Z}{
    The 3-d array of encounter histories. The order of the dimensions is
    animal, trap, occassion.
}
  \item{X}{
    The matrix of trap coordinates. Should have 2 columns for easting
    and northing.
}
  \item{M}{
    The number of animals in the a
}
  \item{niters}{
    Number of MCMC iterations.
}
  \item{xlims}{
    The x-coordinates defining the state-space
}
  \item{ylims}{
    The y-coordinates defining the state-space
}
  \item{tune}{
    A vector of tuning parameters for the Metropolis steps in the
    algorithm. This should contain 3 values for the sigma, lambda_0,
    activity center proposals.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{out}{A matrix of posterior samples with one row per iteration
    and a column for the parameters: sigma, lambda_0, alpha, psi, N}
  \item{last}{The values of the latent parameters w and S for the last
    iteration. These can be used to update the model if it needs to be
    run longer}
}
\note{
  n
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{



# Spatial covariate (with mean 0)
elev.fn <- function(s) {
    s <- matrix(s, ncol=2)        # Force s to be a matrix
    (s[,1] + s[,2] - 100) / 40.8  # Returns (standardized) "elevation"
}

mu <- function(s, beta0, beta1) exp(beta0 + beta1*elev.fn(s=s))

library(R2Cuba)
xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(100,100), beta0=0, beta1=2)

xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta0=0, beta1=2,
            flags=list("verbose"=0))


# Simulate PP using rejection sampling
set.seed(31025)
beta0 <- -6 # intercept of intensity function
beta1 <- 1  # effect of elevation on intensity
# Next line computes integral, which is expected value of N
EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
            lower=c(0,0), upper=c(100,100))$value
EN
N <- rpois(1, EN) # Realized N
s <- matrix(NA, N, 2) # This matrix will hold the coordinates
elev.min <- elev.fn(c(0,0))
elev.max <- elev.fn(c(100, 100))
Q <- max(c(exp(beta0 + beta1*elev.min),
           exp(beta0 + beta1*elev.max)))
counter <- 1
while(counter <= N) {
  x.c <- runif(1, 0, 100); y.c <- runif(1, 0, 100)
  s.cand <- c(x.c,y.c)
  pr <- mu(s.cand, beta0, beta1) #/ EN
  if(runif(1) < pr/Q) {
    s[counter,] <- s.cand
    counter <- counter+1
    }
  }

plot(s)


xsp <- seq(20, 80, by=10); len <- length(xsp)
X <- cbind(rep(xsp, each=len), rep(xsp, times=len)) # traps
ntraps <- nrow(X); noccasions <- 5
y <- array(NA, c(N, ntraps, noccasions)) # capture data
sigma <- 5  # scale parameter
lam0 <- 1   # basal encounter rate
lam <- matrix(NA, N, ntraps)
set.seed(5588)
for(i in 1:N) {
    for(j in 1:ntraps) {
        # The object "s" was simulated in previous section
        distSq <- (s[i,1]-X[j,1])^2 + (s[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(noccasions, lam[i,j])
    }
}
# data augmentation
nz <- 80
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, noccasions))
yz[1:nrow(y),,] <- y # Fill data augmentation array





# Fit the model using MCMC
# Sample the parameters: "sigma", "lam0", "beta0", "beta1", "N", "EN"

\dontrun{
set.seed(3434)
system.time({
fm1 <- scrIPP(yz, X, M, 10000, xlims=c(0,100), ylims=c(0,100),
              space.cov=elev.fn,
              tune=c(0.4, 0.2, 0.3, 0.3, 7))
}) # 328s



library(coda)

mc1 <- mcmc(fm1$out)
plot(mc1)
summary(mc1)
rejectionRate(mc1)
}


}
