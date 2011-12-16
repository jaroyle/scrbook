




# spatial covariate
elev.fn <- function(x) x[,1]+x[,2]


# 2-dimensional integration over [-1, 1] square
int2d <- function(alpha, delta=0.02) {
  z <- seq(-1+delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(alpha*elev.fn(S)) * cell.area)
  }

# Simulate PP using rejection sampling
set.seed(395)
N <- 100
count <- 1
s <- matrix(NA, N, 2)
alpha <- 2 # parameter of interest
while(count <= 100) {
  x.c <- runif(1, -1, 1)
  y.c <- runif(1, -1, 1)
  s.cand <- cbind(x.c,y.c)
  elev.min <- elev.fn(cbind(-1,-1))
  elev.max <- elev.fn(cbind(1,1))
  pr <- exp(alpha*elev.fn(s.cand)) / int2d(alpha)
  Q <- max(c(exp(alpha*elev.min) / int2d(alpha),
             exp(alpha*elev.max) / int2d(alpha)))
  if(runif(1) < pr/Q) {
    s[count,] <- s.cand
    count <- count+1
    cat("accept\n")
    }
  else
    cat("  reject\n")
  }

# plot the simulated data
png("../figs/elevMap.png", width=7, height=7, units="in", res=400)
Sx <- seq(-1, 1, length=50)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn(S)
image(Sx, Sx, matrix(elev, len), col=rgb(0,seq(0.1,1,0.01),0,0.8),
      xlab="Easting", ylab="Northing", main="Elevation map")
box()
points(s, pch=16)
dev.off()


# Negative log-likelihood
nll <- function(beta) {
  -sum(beta*elev.fn(S[,1], S[,2]) - log(int2d(beta)))
  }


starting.value <- 0
fm <- optim(starting.value, nll, method="Brent",
            lower=-5, upper=5, hessian=TRUE)
c(Est=fm$par, SE=sqrt(solve(fm$hessian)))
















# Create trap locations
xsp <- seq(-0.8, 0.8, by=0.2)
len <- length(xsp)
X <- cbind(rep(xsp, each=len),
           rep(xsp, times=len))
str(X)


plot(X, pch="+", xlim=c(-1,1), ylim=c(-1,1))
points(S, pch=16)


# Simulate capture histories
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps, T))

nz <- 50 # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, T))

sigma <- 0.1
lam0 <- 0.5
lam <- matrix(NA, N, ntraps)

for(i in 1:N) {
    for(j in 1:ntraps) {
        distSq <- (S[i,1]-X[j,1])^2 + (S[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(T, lam[i,j])
    }
}
yz[1:nrow(y),,] <- y # Fill

table(y)
summary(c(yz))
table(yz)
sum(rowSums(y)>0)


















# Create trap locations
xsp <- seq(-0.8, 0.8, by=0.2)
len <- length(xsp)
X <- cbind(rep(xsp, each=len),
           rep(xsp, times=len))
str(X)


plot(X, pch="+", xlim=c(-1,1), ylim=c(-1,1))
points(s, pch=16)


# Simulate capture histories
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps, T))

nz <- 50 # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, T))

sigma <- 0.1
lam0 <- 0.5
lam <- matrix(NA, N, ntraps)

for(i in 1:N) {
    for(j in 1:ntraps) {
        distSq <- (s[i,1]-X[j,1])^2 + (s[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(T, lam[i,j])
    }
}
yz[1:nrow(y),,] <- y # Fill

table(y)
summary(c(yz))
table(yz)
sum(rowSums(y)>0)






# MCCM code

source("thinSmcmc_exp.R")
ls()

fm1 <- scrIPP(yz, X, M, 3000, xlims=c(-1,1), ylims=c(-1,1),
            tune=c(0.002, 0.1, 0.25, 0.07) )

library(coda)
plot(mcmc(fm1$out))

rejectionRate(mcmc(fm1$out))















