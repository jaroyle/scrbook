
# Simulation study with both SS covs and ED covs

library(scrbook)
library(raster)
library(gdistance)
library(rjags)


set.seed(9373)
pix <- 5
pixArea <- (pix*pix) / 10000
B <- 100  # Length of square side
dat <- spcov(B=B, pix=pix, cor=10)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","CANHT")
cell <- seq(pix/2, B-pix/2, pix)

canhtMat <- t(matrix(dat$CANHT, B/pix, B/pix))
canht <- flip(raster(t(canhtMat)), direction="y")
names(canht) <- "canht"



# Trap locations
xsp <- seq(27.5, 72.5, by=5)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)




# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps))

sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.8   # basal encounter rate
lam <- matrix(NA, N, ntraps)
theta <- 1

# s defined earlier
#s <- matrix(NA, N, 3)
#colnames(s) <- c("pixID", "x", "y")

set.seed(557828)
cost <- exp(theta*canht)
tr1 <- transition(cost, transitionFunction = function(x) 1/mean(x),
                  directions=8)
tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE)
D <- t(costDistance(tr1CorrC, X, s[,2:3]))
for(i in 1:N) {
#    s defined previously
#    s.i <- sample(1:npix, 1, prob=dat$cp)
#    sx <- dat[s.i, "x"]
#    sy <- dat[s.i, "y"]
#    s[i,] <- c(s.i, sx, sy)
    for(j in 1:ntraps) {
#        distSq <- (sx-X[j,1])^2 + (sy - X[j,2])^2
        lam[i,j] <- exp(-D[i,j]^2/(2*sigma^2)) * lam0
        y[i,j] <- rpois(1, lam[i,j])
    }
}

sum(y)



y.ded <- y[rowSums(y)>0,]
str(y.ded)

(fm1 <- scrDED(y.ded, X, ~1, ~1, rasters=canht,
               start=c(-1, -1, 1),
               method="BFGS",
#               lower=c(-3,-3,-3), upper=c(5,1,5),
               control=list(trace=TRUE, REPORT=1, maxit=50)))

exp(fm1$par[1:2])            # 0.8, 0.1
exp(fm1$par[3])+nrow(y.ded)  # 50


(fm2 <- scrDED(y.ded, X, ~canht, ~1, rasters=canht,
#               start=c(log(0.8), log(0.1), log(10), 2),
               method="BFGS",
               control=list(trace=TRUE, REPORT=1, maxit=500)))

exp(fm2$par[1:2])               # 0.8, 0.1
exp(fm2$par[3])+nrow(y.ded)     # 50



(fm3 <- scrDED(y.ded, X, ~canht, ~canht, rasters=canht,
#               start=c(log(0.8), log(0.1), log(10), 2, 0),
               method="BFGS",
               control=list(trace=TRUE, REPORT=1, maxit=500)))

exp(fm3$par[1:2])                       # 0.8, 0.1
c(N=exp(fm3$par[3])+nrow(y.ded))        # 50
fm3$par[4:5]                            # 2, 0


















# Simulation study



set.seed(353)
pix <- 0.05
dat <- spcov(pix=pix)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","canht")
cell <- seq(pix/2, 1-pix/2, pix)
image(cell, cell, t(matrix(dat$canht, 1/pix, 1/pix)), ann=FALSE)

head(dat)

# Trap locations
xsp <- seq(0.275, 0.725, by=0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)

# Canhtation covariate as a matrix, then as a raster
canhtMat <- t(matrix(dat$canht, 1/pix, 1/pix))
canht <- flip(raster(t(canhtMat)), direction="y")
names(canht) <- "canht"

sim.data <- function(N=50, sigma=0.1, lam0=0.8, beta=1, theta=1, X,
                     covar) {
    ntraps <- nrow(X)
    y <- lam <- matrix(NA, N, ntraps)

    # Activity centers
    s <- matrix(NA, N, 3)
    colnames(s) <- c("pixID", "x", "y")

    npix <- ncell(covar)

    den <- exp(beta*dat$canht)
    den <- den/sum(den)

    covar.tran <- covar-cellStats(covar, min)
    covar.tran <- covar.tran/cellStats(covar.tran, max)

    cost <- exp(theta*covar.tran)
    tr1 <- transition(cost, transitionFunction = function(x) 1/mean(x),
                      directions=8)
    tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE)
    for(i in 1:N) {
        s.i <- sample(1:npix, 1, prob=den)
        sx <- dat[s.i, "x"]
        sy <- dat[s.i, "y"]
        s[i,] <- c(s.i, sx, sy)
        distSq <- costDistance(tr1CorrC, X, s[i,2:3])^2
        lam[i,] <- exp(-distSq/(2*sigma*sigma)) * lam0
        y[i,] <- rpois(nrow(X), lam[i,])
    }
    y <- y[rowSums(y)>0,]
    return(y)
}


sum(sim.data(X=X, covar=canht))



nsim <- 2
simout <- matrix(NA, nsim, 5)
colnames(simout) <- c("lam0", "sigma", "N", "beta", "theta")
set.seed(343)
for(i in 1:nsim) {
    cat("doing", i, format(Sys.time(), "%H:%M:%S"), "\n")
    lam0 <- 2
    sigma <- 0.1
    N <- 50
    beta <- 1
    theta <- 1
    y.i <- sim.data(N=N, sigma=sigma, lam0=lam0, beta=beta, theta=theta,
                    X=X, covar=canht)
    cat("  ncaught =", nrow(y.i), "\n")
    fm.i <- scrDED(y=y.i, traplocs=X, ~canht, ~canht, rasters=canht,
                   start=c(0, -2, 3, 1, 1),
                   method="BFGS", control=list(trace=TRUE, REPORT=5))
    mle <- fm.i$par
    simout[i,] <- c(exp(mle[1:2]), exp(mle[3])+nrow(y.i), mle[4:5])
    cat("  mle =", simout[i,], "\n\n")
}


png("figs/scrDEDsim.png", width=6, height=3, units="in", res=400)
op <- par(mfrow=c(1,3), mai=c(0.6,0.3,0.3,0.2))
#hist(simout[,1], main="", xlab="lam0", freq=FALSE, cex.lab=1.2)
#abline(v=lam0, lwd=3, col=4, lty=1)
#hist(simout[,2], main="", xlab="sigma", freq=FALSE, cex.lab=1.2)
#abline(v=sigma, lwd=3, col=4, lty=1)
hist(simout[,3], main="", xlab="population size (N)", freq=FALSE,
     cex.lab=1.2)
abline(v=N, lwd=3, col=4, lty=1)
hist(simout[,4], main="", xlab="Density effect (beta)",
     freq=FALSE, cex.lab=1.2)
abline(v=beta, lwd=3, col=4, lty=1)
hist(simout[,5], main="",
     xlab="Ecological distance effect (theta)",
     freq=FALSE, cex.lab=1.2)
abline(v=theta, lwd=3, col=4, lty=1)
par(op)
dev.off()

ytest <- sim.data(X=X, covar=canht, beta=0)
fm.test <- scrDED(y=y.i, traplocs=X, ~1, ~canht, rasters=canht,
                  start=c(0, -2, 2, 1),
                  method="BFGS", control=list(trace=TRUE, REPORT=1))








































combn(4,2)

choose(4,2)
choose(100, 3)

