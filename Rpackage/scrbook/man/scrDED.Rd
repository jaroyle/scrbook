\name{scrDED}
\alias{scrDED}
\title{
Fit SCR model with state-space covariates and covariates of ecological distance.
}
\description{
Likelihood-based method.
}
\usage{
scrDED(y = y, traplocs = traplocs, den.formula = ~1, dist.formula = ~1, rasters, start, ...)
}
\arguments{
  \item{y}{
    An nAnimal by nTrap encounter history matrix
}
  \item{traplocs}{
    An nTrap by 2 matrix of trap coordinates (x then y)
}
  \item{den.formula}{
    Right-hand sided formula
}
  \item{dist.formula}{
    Righ-hand sided formula
}
  \item{rasters}{
    An object of class RasterLayer or RasterStack containing the
    covariate data. This determines the extent of the state-space and
    the resolution used for the discrete sum approximation of the
    integral over space.
}
  \item{start}{
    An optional vector of starting values. The order is lam0, sigma, n0,
    density covs, and then ecological distance covs
}
  \item{\dots}{
    Optional arguments to pass to \code{\link{optim}}
}
}
\details{
  Right now, this assumes a Poisson encounter history model. Need to
  allow for binomial model too.
}
\value{
  Currently, it just returns the object spit out by optim.
}
\note{
  Note
}

\examples{

library(raster)
library(gdistance)

set.seed(353)
pix <- 0.05
dat <- spcov(pix=pix)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","elev")
cell <- seq(pix/2, 1-pix/2, pix)
image(cell, cell, t(matrix(dat$elev, 1/pix, 1/pix)), ann=FALSE)

head(dat)

# Simulate IPP
set.seed(30275)
N <- 50
alpha <- 2
dat$cp <- exp(alpha*dat$elev) / sum(exp(alpha*dat$elev))
s.tmp <- rmultinom(1, N, dat$cp) # a single realization to be ignored later

# Trap locations
xsp <- seq(0.275, 0.725, by=0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)

# Elevation covariate as a matrix, then as a raster
elevMat <- t(matrix(dat$elev, 1/pix, 1/pix))
elev <- flip(raster(t(elevMat)), direction="y")
layerNames(elev) <- "elev"

# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps))

# Parameters
sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.8   # basal encounter rate
lam <- matrix(NA, N, ntraps)
theta <- 1

# Activity centers
s <- matrix(NA, N, 3)
colnames(s) <- c("pixID", "x", "y")

set.seed(557828)
elev.tran <- elev-cellStats(elev, min)
elev.tran <- elev.tran/cellStats(elev.tran, max)
cost <- exp(theta*elev.tran)
tr1 <- transition(cost, transitionFunction = function(x) 1/mean(x),
                  directions=8)
tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE)
for(i in 1:N) {
    s.i <- sample(1:npix, 1, prob=dat$cp)
    sx <- dat[s.i, "x"]
    sy <- dat[s.i, "y"]
    s[i,] <- c(s.i, sx, sy)
    distSq <- costDistance(tr1CorrC, X, s[i,2:3])^2
    lam[i,] <- exp(-distSq/(2*sigma*sigma)) * lam0
    y[i,] <- rpois(nrow(X), lam[i,])
}

sum(y)



y.ded <- y[rowSums(y)>0,]
str(y.ded)


\dontrun{

(fm1 <- scrDED(y.ded, X, ~1, ~1, rasters=elev,
#               start=c(-1, -1, 1),
               method="BFGS",
               control=list(trace=TRUE, REPORT=1, maxit=50)))

exp(fm1$par[1:2])            # 0.8, 0.1
exp(fm1$par[3])+nrow(y.ded)  # 50


(fm2 <- scrDED(y.ded, X, ~elev, ~1, rasters=elev,
#               start=c(log(0.8), log(0.1), log(10), 2),
               method="BFGS",
               control=list(trace=TRUE, REPORT=1, maxit=500)))

exp(fm2$par[1:2])               # 0.8, 0.1
exp(fm2$par[3])+nrow(y.ded)     # 50



(fm3 <- scrDED(y.ded, X, ~elev, ~elev, rasters=elev,
#               start=c(log(0.8), log(0.1), log(10), 2, 0),
               method="BFGS",
               control=list(trace=TRUE, REPORT=1, maxit=500)))

exp(fm3$par[1:2])                       # 0.8, 0.1
c(N=exp(fm3$par[3])+nrow(y.ded))        # 50
fm3$par[4:5]                            # 2, 1


}

}
% R documentation directory.
