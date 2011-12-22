\name{ch9secrYjags}
\alias{ch9jags}
\title{
  Example of fitting inhomogeneous point process using JAGS
}
\description{
  The analysis of is of data simulated on a unrealistically small 5x5
  raster. See Chapter 9 for an application using real data
}
\details{
  See Chapter 9 or is it 10?
}
\seealso{
  \code{\link{scrIPP}}
}
\examples{

# spatial covariate (with mean 0)
elev.fn <- function(x) x[,1]+x[,2]-1

# Elevation on a coarse grid
Sx <- seq(0.1, 0.9, 0.2)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn(S)

# 2-dimensional integration over unit square
int2d <- function(alpha, delta=0.02) {
  z <- seq(delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(alpha*elev.fn(S)) * cell.area)
  }

# Simulate IPP using rejection sampling
set.seed(300225)
N <- 100
count <- 1
s <- matrix(NA, N, 2)
alpha <- 2 # parameter of interest
while(count <= 100) {
  x.c <- runif(1, 0, 1); y.c <- runif(1, 0, 1)
  s.cand <- cbind(x.c,y.c)
  elev.min <- elev.fn(cbind(0,0))
  elev.max <- elev.fn(cbind(1,1))
  pr <- exp(alpha*elev.fn(s.cand)) / int2d(alpha)
  Q <- max(c(exp(alpha*elev.min) / int2d(alpha),
             exp(alpha*elev.max) / int2d(alpha)))
  if(runif(1) < pr/Q) {
    s[count,] <- s.cand
    count <- count+1
    }
  }

# Counts in each grid cell
n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))

# 3x3 grid of traps centered
xc <- seq(0.3, 0.7, 0.2)
X <- cbind(rep(xc, each=3), rep(xc, times=3))

# Parameters
lam0 <- 2
sigma <- 0.15

# Generate encounter histories
set.seed(32235)
y <- matrix(NA, nrow(s), nrow(X))
for(i in 1:nrow(s)) {
    for(j in 1:nrow(X)) {
        dist <- sqrt((s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2)
        lam <- lam0 * exp(-dist^2/(2*sigma*sigma))
        y[i,j] <- rpois(1, lam)
    }
}

# Augment data
nz <- 50
yz <- matrix(0, nrow(y)+nz, ncol(y))
yz[1:nrow(y),] <- y


# BUGS model
sink("ippDiscrete.txt")
cat("
model{
sigma ~ dunif(0, 1)
lam0 ~ dunif(0, 5)
beta ~ dnorm(0,0.1)
psi ~ dbeta(1,1)

for(j in 1:nPix) {
  theta[j] <- exp(beta*elevation[j])
  probs[j] <- theta[j]/sum(theta[])
}

for(i in 1:M) {
  w[i] ~ dbern(psi)
  s[i] ~ dcat(probs[])
  x0g[i] <- Sgrid[s[i],1]
  y0g[i] <- Sgrid[s[i],2]
  for(j in 1:ntraps) {
    dist[i,j] <- sqrt(pow(x0g[i]-grid[j,1],2) + pow(y0g[i]-grid[j,2],2))
    lambda[i,j] <- lam0*exp(-dist[i,j]*dist[i,j]/(2*sigma*sigma)) * w[i]
    y[i,j] ~ dpois(lambda[i,j])
    }
  }

N <- sum(w[])
D <- N/1 # unit square
}

", fill=TRUE)
sink()


\dontrun{
library(rjags)

modfile <- "ippDiscrete.txt"
file.show(modfile)
dat <- list(y=yz, elevation=elev, nPix=prod(dim(n.k)),
            M=nrow(yz), ntraps=nrow(X), Sgrid=S, grid=X)
set.seed(434)
init <- function() {
    list(sigma=runif(1), lam0=runif(1), beta=rnorm(1),
         s=sample.int(dat$nPix, dat$M, replace=TRUE),
         w=c(rep(1,100), rep(1,nz)), psi=1)
}
pars <- c("sigma", "lam0", "beta", "N")

# Compile and adapt
jm <- jags.model(modfile, dat, init, n.chains=2, n.adapt=500)
# MCMC
jc <- coda.samples(jm, pars, n.iter=2000)

plot(jc)
summary(window(jc, start=1001))
}

unlink(modfile)


}
