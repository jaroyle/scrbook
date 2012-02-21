#test

# Homogeneous BPP


set.seed(4234)
N <- 50
s <- cbind(runif(N), runif(N)) # points
# counts in each pixel
n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))


# plot continuous space and discrete space
png("figs/homoPlots.png", width=5, height=2.5, units="in", res=400)
op <- par(mfrow=c(1, 2), mai=c(0.1, 0.1, 0.1, 0.1))
plot(s, frame=T, ann=FALSE, axes=FALSE, asp=1, cex=0.5)
#segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
#segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
plot(s, frame=T, ann=FALSE, axes=F, type="n", asp=1)
segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
y <- 0.1
for(i in 1:nrow(n.k)) {
    text(seq(0.1, 1, by=0.2), y, labels=n.k[,i], cex=0.6)
    y <- y+0.2
}
par(op)
dev.off()



















# Heterogeneous BPP





# Create a spatial covariate





# spatial covariate (with mean 0)
#elev.fn <- function(x) x[,1]+x[,2]-1
elev.fn <- function(x) x[1]+x[2]-1

mu <- function(x, beta) exp(beta*elev.fn(x=x))

library(R2Cuba)
xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta=2)

xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta=2,
            flags=list("verbose"=0))


# 2-dimensional integration over unit square
int2d <- function(beta, delta=0.02) {
  z <- seq(delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
#  sum(exp(beta*elev.fn(S)) * cell.area)
  sum(exp(beta*(S[,1]+S[,2]-1)) * cell.area)
  }

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
#  int.mu <- cuhre(2, 1, mu, beta=beta)$value
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






n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))


# plot continuous space and discrete space
png("figs/heteroPlots.png", width=5, height=2.5, units="in", res=400)
op <- par(mfrow=c(1, 2), mai=c(0.1, 0.1, 0.1, 0.1))
Sx <- seq(0.01, 0.99, 0.01)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev.fn2 <- function(x) x[,1]+x[,2]-1
elev <- elev.fn2(S)
image(Sx, Sx, matrix(elev, len), col=rgb(0,seq(0.1,1,0.01),0,0.8),
      ann=FALSE, axes=FALSE, asp=1)
points(s, cex=0.5)
#segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
#segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
box(col=gray(0.5))
Sx <- seq(0.1, 0.9, 0.2)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn2(S)
image(Sx, Sx,
      matrix(elev, len), xlim=c(0,1), ylim=c(0,1),
      col=rgb(0,seq(0.1,1,0.01),0,0.8),
      ann=FALSE, axes=FALSE, asp=1)
segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
box(col=gray(0.5))
y <- 0.1
for(i in 1:nrow(n.k)) {
    text(seq(0.1, 1, by=0.2), y, labels=n.k[,i], cex=0.6)
    y <- y+0.2
}
par(op)
dev.off()




























# Analysis using custom MCMC





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








library(scrbook)
set.seed(3434)
fm1 <- scrIPP(yz, X, M, 6000, xlims=c(0,1), ylims=c(0,1),
            tune=c(0.003, 0.08, 0.3, 0.07) )

plot(mcmc(fm1$out))
rejectionRate(mcmc(fm1$out))


fm1.s <- summary(window(mcmc(fm1$out), start=1001))
fm1.r <- cbind(fm1.s$stat[,1:2], fm1.s$quant[,c(1,3,5)])
colnames(fm1.r) <- c("& Mean", "SD", "2.5\\%", "50\\%", "97.5\\%")
rownames(fm1.r) <- c("$\\sigma=0.1$", "$\\lambda_0=0.5$", "$\\psi=0.66$",
                     "$\\beta=2$", "$N=100$")
round(fm1.r,2)




sink("fm1.r.tex")
cat("
\\begin{table}
\\centering
\\caption{Posterior summaries from inhomogeneous point proces model}
\\begin{tabular}{lrrrrr}
\\hline
")
write.table(format(fm1.r, digits=2, nsmall=4, scientific=FALSE),
            quote=FALSE, sep=" & ", eol=" \\\\\n ")
cat("
\\hline
\\end{tabular}
\\label{ch9:tab:simIPP}
\\end{table}
")
sink()
















# Discrete space


library(scrbook)
library(secr)
library(raster) # raster after secr b/c of flip
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


elevMat <- t(matrix(dat$elev, 1/pix, 1/pix))
library(raster)
elev <- flip(raster(t(elevMat)), direction="y")

windows(width=3, height=6)
op <- par(mfrow=c(2,1), mai=rep(0.2,4))
#png("figs/discrete.png", width=7, height=7, units="in", res=400)
image(cell, cell, elevMat, ann=FALSE)
points(dat[s.tmp>0,c("x","y")], cex=s.tmp[s.tmp>0])
points(X, pch="+")
box()
#dev.off()
plot(elev)
par(op)




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

set.seed(557828)
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





# Analysis using secr
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
#plot(ch, tol=0.0005) # ouch

# Make mask

msk <- make.mask(secr.traps, buffer=0.275, spacing=.05, nx=v)
summary(msk)
#plot(msk)

ssArea <- attr(msk, "area")*nrow(msk)

covariates(msk) <- data.frame(elev=dat$elev[order(dat$y, dat$x)])




ch9simData <- list(ch.secr=ch, ch.jags=yz, spcov.jags=dat, spcov.secr=msk,
                   traps=X)


#save(ch9simData, file="../Rpackage/scrbook/data/ch9simData.rda")
#promptData(ch9simData, filename="../Rpackage/scrbook/man/ch9simData.Rd")




library(scrbook)
library(secr)
library(rjags)

data(ch9simData)

ch <- ch9simData$ch.secr
msk <- ch9simData$spcov.secr


# SECR analysis

secr1 <- secr.fit(ch, model=D~elev, mask=msk)

region.N(secr1, se.N=TRUE)





# JAGS analysis

# JAGS model
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



modfile <- "ippDiscrete.txt"

jags.data <- with(ch9simData, {
    list(y=ch.jags, elevation=drop(spcov.jags$elev),
            nPix=nrow(spcov.jags),
            M=nrow(ch.jags), ntraps=nrow(traps),
            Sgrid=as.matrix(spcov.jags[,1:2]),
            grid=traps)
    })
str(jags.data)

all(matrix(jags.data$elevation, 20, byrow=T) ==
    matrix(covariates(msk)$elev, 20))

init <- function() {
    list(sigma=runif(1), lam0=runif(1), beta=rnorm(1),
         s=sample.int(jags.data$nPix, jags.data$M, replace=TRUE),
         w=rep(1, jags.data$M), psi=1)
}
str(init())

pars <- c("sigma", "lam0", "beta", "N")

# Obtain posterior samples. This takes a few minutes
# Compile and adapt
system.time({
set.seed(03453)
jm <- jags.model(modfile, jags.data, init, n.chains=2, n.adapt=1000)
jags1 <- coda.samples(jm, pars, n.iter=10000)
})

plot(jags1)
summary(jags1)

summary(window(jags1, start=3501))
plot(window(jags1, start=3501))



unlink(modfile)




save.image("scratch.RData")







# compare results

library(scrbook)
#example(ch9secrYjags)


plot(window(jags1, start=5001))

summary(window(jags1, start=5001))
gelman.diag(window(jags1, start=5001))


jags.est <- summary(window(jags1, start=5001))
jags.r <- cbind(jags.est$stat[,1:2], jags.est$quant[,c(1,5)])

secr.est <- predict(secr1)
secr.r <- cbind(secr.est[2:3,2:5])
secr.r <- rbind(beta=as.numeric(coef(secr1)[2,]), secr.r)
secr.r <- rbind(region.N(secr1)[2,1:4], secr.r)


comp.out <- data.frame(rbind(as.matrix(secr.r), jags.r))
comp.out <- cbind(Software=c(rep("secr",4), rep("JAGS",4)),
      Parameter=c("$N$", "$\\beta$", "$\\lambda_0$", "$\\sigma$",
                  "$N$", "$\\beta$", "$\\lambda_0$", "$\\sigma$"),
      comp.out)

colnames(comp.out) <- c("Software", "Par", "Est.", "SD", "lower", "upper")


sink("comp.out.tex")
cat("
\\begin{table}
\\centering
\\caption{Comparision of \\jags~and \\secr~results}
\\begin{tabular}{llrrrr}
\\hline
")
write.table(format(comp.out, digits=2, nsmall=4, scientific=FALSE),
            quote=FALSE, sep=" & ", eol=" \\\\\n ", row.names=FALSE)
cat("
\\hline
\\end{tabular}
\\label{ch9:tab:simIPP}
\\end{table}
")
sink()



















# Consider both SS covs and ED covs














elevMat <- t(matrix(dat$elev, 1/pix, 1/pix))
elev <- flip(raster(t(elevMat)), direction="y")
layerNames(elev) <- "elev"







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
cost <- exp(theta*elev)
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

(fm1 <- scrDED(y.ded, X, ~1, ~1, rasters=elev,
               start=c(-1, -1, 1),
               method="BFGS",
#               lower=c(-3,-3,-3), upper=c(5,1,5),
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
fm3$par[4:5]                            # 2, 0


















# Simulation study


library(scrbook)
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

# Trap locations
xsp <- seq(0.275, 0.725, by=0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)

# Elevation covariate as a matrix, then as a raster
elevMat <- t(matrix(dat$elev, 1/pix, 1/pix))
elev <- flip(raster(t(elevMat)), direction="y")
layerNames(elev) <- "elev"

sim.data <- function(N=50, sigma=0.1, lam0=0.8, beta=1, theta=1, X,
                     covar) {
    ntraps <- nrow(X)
    y <- lam <- matrix(NA, N, ntraps)

    # Activity centers
    s <- matrix(NA, N, 3)
    colnames(s) <- c("pixID", "x", "y")

    npix <- ncell(covar)

    den <- exp(beta*dat$elev)
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


sum(sim.data(X=X, covar=elev))



nsim <- 500
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
                    X=X, covar=elev)
    cat("  ncaught =", nrow(y.i), "\n")
    fm.i <- scrDED(y=y.i, traplocs=X, ~elev, ~elev, rasters=elev,
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

ytest <- sim.data(X=X, covar=elev, beta=0)
fm.test <- scrDED(y=y.i, traplocs=X, ~1, ~elev, rasters=elev,
                  start=c(0, -2, 2, 1),
                  method="BFGS", control=list(trace=TRUE, REPORT=1))















































