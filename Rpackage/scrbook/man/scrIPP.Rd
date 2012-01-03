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

# spatial covariate (with mean 0)
#elev.fn <- function(x) x[,1]+x[,2]-1
elev.fn <- function(x) x[1]+x[2]-1

mu <- function(x, beta) exp(beta*elev.fn(x=x))

library(R2Cuba)

# Simulate PP using rejection sampling
set.seed(300225)
N <- 100
count <- 1
s <- matrix(NA, N, 2)
beta <- 2 # parameter of interest
int.mu <- cuhre(2, 1, mu, beta=beta)$value
elev.min <- elev.fn(c(0,0)) #elev.fn(cbind(0,0))
elev.max <- elev.fn(c(1,1)) #elev.fn(cbind(1,1))
Q <- max(c(exp(beta*elev.min) / int.mu,   #2d(beta),
           exp(beta*elev.max) / int.mu))   #2d(beta)))
while(count <= 100) {
  x.c <- runif(1, 0, 1); y.c <- runif(1, 0, 1)
  s.cand <- c(x.c,y.c)
  pr <- exp(beta*elev.fn(s.cand)) / int.mu #2d(beta)
  if(runif(1) < pr/Q) {
    s[count,] <- s.cand
    count <- count+1
    }
  }


# Maximum likelihood
nll <- function(beta) {
    int.mu <- cuhre(2, 1, mu, beta=beta)$value
    -sum(beta*elev.fn(s) - log(int.mu))
}
starting.value <- 0
fm <- optim(starting.value, nll, method="Brent",
            lower=-5, upper=5, hessian=TRUE)
c(Est=fm$par, SE=sqrt(1/fm$hessian)) # estimates and SEs


# Counts in discrete space... not used for anything
n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))











# Create trap locations
xsp <- seq(0.2, 0.8, by=0.1)
len <- length(xsp)
X <- cbind(rep(xsp, each=len), rep(xsp, times=len))

# Simulate capture histories, and augment the data
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps, T))

nz <- 50 # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, T))

sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.5   # basal encounter rate
lam <- matrix(NA, N, ntraps)

set.seed(5588)
for(i in 1:N) {
    for(j in 1:ntraps) {
        distSq <- (s[i,1]-X[j,1])^2 + (s[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(T, lam[i,j])
    }
}
yz[1:nrow(y),,] <- y # Fill






# MCMC

set.seed(3434)
fm1 <- scrIPP(yz, X, M, 6000, xlims=c(0,1), ylims=c(0,1),
            tune=c(0.003, 0.08, 0.3, 0.07) )

plot(mcmc(fm1$out))
rejectionRate(mcmc(fm1$out))



}
