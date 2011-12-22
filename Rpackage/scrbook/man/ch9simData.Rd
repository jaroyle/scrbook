\name{ch9simData}
\alias{ch9simData}
\docType{data}
\title{
  Simulated data from an inhomogeneous point process model
}
\description{
  Data simulated from an inhomogeneous point process model. Used to
  compare secr and JAGS in Chapter 9.
}
\usage{data(ch9simData)}
\format{
  The format is:
List of 5
 $ ch.secr   : capthist [1:31, 1, 1:81] 0 0 1 1 0 0 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : chr [1:31] "4" "5" "6" "7" ...
  .. ..$ : chr "1"
  .. ..$ : chr [1:81] "1" "2" "3" "4" ...
  ..- attr(*, "covariates")='data.frame':	0 obs. of  0 variables
  ..- attr(*, "traps")=Classes 'traps' and 'data.frame':	81 obs. of  2 variables:
  .. ..$ x: num [1:81] 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.35 ...
  .. ..$ y: num [1:81] 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.3 ...
  .. ..- attr(*, "detector")= chr "count"
  .. ..- attr(*, "spacex")= num 0.05
  .. ..- attr(*, "spacey")= num 0.05
  .. ..- attr(*, "spacing")= num 0.05
  ..- attr(*, "session")= chr "1"
 $ ch.jags   : num [1:100, 1:81] 0 0 0 0 0 1 1 0 0 0 ...
 $ spcov.jags:'data.frame':	441 obs. of  4 variables:
  ..$ x   : num [1:441] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ y   : num [1:441] 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...
  ..$ elev: num [1:441, 1] 0.643 0.441 0.641 0.573 0.317 ...
  ..$ cp  : num [1:441, 1] 0.00545 0.00364 0.00543 0.00473 0.00284 ...
 $ spcov.secr:Classes 'mask' and 'data.frame':	441 obs. of  2 variables:
  ..$ x: num [1:441] -2.08e-17 5.00e-02 1.00e-01 1.50e-01 2.00e-01 ...
  ..$ y: num [1:441] -2.08e-17 -2.08e-17 -2.08e-17 -2.08e-17 -2.08e-17 ...
  ..- attr(*, "type")= chr "traprect"
  ..- attr(*, "meanSD")='data.frame':	2 obs. of  2 variables:
  .. ..$ x: num [1:2] 0.5 0.303
  .. ..$ y: num [1:2] 0.5 0.303
  ..- attr(*, "area")= num 2.5e-07
  ..- attr(*, "spacing")= num 0.05
  ..- attr(*, "boundingbox")='data.frame':	4 obs. of  2 variables:
  .. ..$ x: num [1:4] -0.025 1.025 1.025 -0.025
  .. ..$ y: num [1:4] -0.025 -0.025 1.025 1.025
  .. ..- attr(*, "out.attrs")=List of 2
  .. .. ..$ dim     : Named int [1:2] 2 2
  .. .. .. ..- attr(*, "names")= chr [1:2] "x" "y"
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ x: chr [1:2] "x=-0.025" "x= 1.025"
  .. .. .. ..$ y: chr [1:2] "y=-0.025" "y= 1.025"
  ..- attr(*, "covariates")='data.frame':	441 obs. of  1 variable:
  .. ..$ elev: num [1:441] 0.643 0.561 0.973 1.064 0.756 ...
 $ traps     : num [1:81, 1:2] 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.35 ...
}
\examples{
data(ch9simData)

# Here is how the data were generated

# state-space covariate
set.seed(3453)
v <- 21
dat <- spcov(v=v)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","elev")
image(t(matrix(dat$elev, v, v)))

# Multinomial cell probs
set.seed(300225)
N <- 50
alpha <- 2
dat$cp <- exp(alpha*dat$elev) / sum(exp(alpha*dat$elev))
s.tmp <- rmultinom(1, N, dat$cp) # a single realization to be ignored later
image(t(matrix(dat$elev, v, v)))
points(dat[s.tmp>0,c("x","y")])

# Trap locations
xsp <- seq(0.3, 0.7, 0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))


# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps))

nz <- 50 # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps))

sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.8   # basal encounter rate
lam <- matrix(NA, N, ntraps)

s <- matrix(NA, N, 3)
colnames(s) <- c("pixID", "x", "y")

set.seed(5588)
for(i in 1:N) {
    s.i <- sample(1:npix, 1, prob=dat$cp)
    sx <- dat[s.i, "x"]
    sy <- dat[s.i, "y"]
    s[i,] <- c(s.i, sx, sy)
    for(j in 1:ntraps) {
        distSq <- (sx-X[j,1])^2 + (sy - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j] <- rpois(1, lam[i,j])
    }
}
yz[1:nrow(y),] <- y # Fill

sum(y)





# Format using secr

# Create a "traps" object
Xs <- data.frame(X)
colnames(Xs) <- c("x","y")
secr.traps <- read.traps(data=Xs, detector="count")

# Create a "capthist" object
secr.caps <- matrix(NA, sum(y), 5)
colnames(secr.caps) <- c("Session", "ID", "Occasion", "X", "Y")
counter <- 0
for(i in 1:nrow(y)) {
    for(j in 1:ncol(y)) {
        y.ij <- y[i,j]
        if(y.ij==0)
            next
        for(v in 1:y.ij) {
            counter <- counter+1
            secr.caps[counter,] <- c(1, i, 1, X[j,1], X[j,2])
        }
    }
}
ch <- make.capthist(secr.caps, secr.traps, fmt="XY")
plot(ch, tol=0.0005) # ouch

# Make mask
msk <- make.mask(secr.traps, buffer=0.325, spacing=.05, nx=v)
summary(msk)
plot(msk)

ssArea <- attr(msk, "area")*nrow(msk)

covariates(msk) <- data.frame(elev=dat$elev[order(dat$y, dat$x)])

ch9simData <- list(ch.secr=ch, ch.jags=yz, spcov.jags=dat, spcov.secr=msk,
                   traps=X)
}
\keyword{datasets}
