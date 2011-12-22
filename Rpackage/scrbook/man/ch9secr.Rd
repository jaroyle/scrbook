\name{ch9secr}
\alias{ch9secr}
\title{
  Example of fitting inhomogeneous point process using secr
}
\description{
  Analysis of simulated data using secr
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



spcov(B=1, v=30)








library(secr)


# Create a "traps" object
Xs <- data.frame(X)
colnames(Xs) <- c("x","y")
secr.traps <- read.traps(data=Xs, detector="count")

summary(secr.traps)

# Huh?
plot(secr.traps)

plot.default(secr.traps, xlim=c(0,1), asp=1, pch="+")

# Create a "capthist" object
secr.caps <- matrix(NA, sum(y), 5)
colnames(secr.caps) <- c("Session", "ID", "Occasion", "X", "Y")
counter <- 0
for(i in 1:nrow(y)) {
    for(j in 1:ncol(y)) {
        for(k in 1:dim(y)[3]) {
            y.ijk <- y[i,j,k]
            if(y.ijk==0)
                next
            for(v in 1:y.ijk) {
                counter <- counter+1
                secr.caps[counter,] <- c(1, i, k, X[j,1], X[j,2])
            }
        }
    }
}
ch <- make.capthist(secr.caps, secr.traps, fmt="XY")
plot(ch, tol=0.0005) # ouch

# Make mask

msk <- make.mask(secr.traps, buffer=0.2, spacing=.1, nx=100)
summary(msk)
plot(msk)

summary(elev.fn(msk))

covariates(msk) <- data.frame(elev=elev.fn(msk))

m0 <- secr.fit(ch, mask=msk, start=c(log(N), log(lam0), log(sigma)))

m1 <- secr.fit(ch, model=D~elev,
               mask=msk, start=c(log(N), 0, log(lam0), log(sigma)))


m1


expected.n(m1)
region.N(m1, se.N=TRUE)






}
