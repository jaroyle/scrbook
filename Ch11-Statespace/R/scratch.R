
# PART I. Misc code shown in chapter


set.seed(4234)
N <- 50
s <- cbind(runif(N), runif(N)) # points
# counts in each pixel
n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))


# plot continuous space and discrete space
png("../figs/homoPlots.png", width=5, height=2.5, units="in", res=400)
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
system("open ../figs/homoPlots.png")



set.seed(454500)
n.Bk <- rmultinom(1, size=50, prob=rep(1/25, 25))
matrix(n.Bk, 5, 5)






set.seed(340)
Area <- 1                  # Area of state-space
M <- 100                   # Data augmentation size
mu <- 10                   # Intensity (points per area)
psi <- (mu*Area)/M         # Data augmentation parameter (thinning rate)
N <- rbinom(1, M, psi)     # Realized value of N under binomial prior
cbind(runif(N), runif(N))  # Point pattern from thinned binomial model











# PART II. Fitting IPP models when points are observed


elev.fn <- function(x) {          # spatial covriate
    x <- matrix(x, ncol=2)        # Force x to be a matrix
    (x[,1] + x[,2] - 100) / 40.8  # Returns (standardized) "elevation"
}
# intensity function
mu <- function(x, beta0, beta1) exp(beta0 + beta1*elev.fn(x=x))
beta0 <- -6 # intercept of intensity function
beta1 <- 1  # effect of elevation on intensity
# Next line computes integral
EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
            lower=c(0,0), upper=c(100,100))$value



# Spatial covariate (with mean 0)
#elev.fn <- function(s) {
#    s <- matrix(s, ncol=2)        # Force s to be a matrix
#    (s[,1] + s[,2] - 100) / 40.8  # Returns (standardized) "elevation"
#}

#mu <- function(s, beta0, beta1) exp(beta0 + beta1*elev.fn(s=s))

library(R2Cuba)
xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(100,100), beta0=0, beta1=2)

xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta0=0, beta1=2,
            flags=list("verbose"=0))

# 2-dimensional integration over unit square
int2d <- function(beta0=0, beta1=2, delta=0.02) {
  z <- seq(delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(beta0 + beta1*elev.fn(S)) * cell.area)
  }

int2d(2, delta=0.001)
int2d(2, delta=0.01)
int2d(2, delta=0.02)
int2d(2, delta=0.03)

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

# Maximum likelihood
nll <- function(beta) {
    beta0 <- beta[1]
    beta1 <- beta[2]
    EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
                lower=c(0,0), upper=c(100,100))$value
    -(sum(beta0 + beta1*elev.fn(s)) - EN)
}
starting.values <- c(-10, 0)
fm <- optim(starting.values, nll, hessian=TRUE)
cbind(Est=fm$par, SE=sqrt(diag(solve(fm$hessian)))) # estimates and SEs


# Now with binomial instead of Poisson prior
nllBin <- function(beta, M=10000) {
    beta0 <- beta[1]
    beta1 <- beta[2]
    EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
                flags=list(verbose=0),
                lower=c(0,0), upper=c(100,100))$value
    N <- nrow(s)
    psi <- EN/M
    -(sum(beta0 + beta1*elev.fn(s) - log(EN)) +
      dbinom(N, M, psi, log=TRUE))
}
starting.values <- c(-10, 0)
fmBin <- optim(starting.values, nllBin, hessian=TRUE)
cbind(Est=fmBin$par, SE=sqrt(diag(solve(fmBin$hessian)))) # est and SE






n.k <- table(cut(s[,1], seq(0, 100, 20)),
             cut(s[,2], seq(0, 100, 20)))


# plot continuous space and discrete space
png("../figs/heteroPlots.png", width=5, height=2.5, units="in", res=400)
op <- par(mfrow=c(1, 2), mai=c(0.1, 0.1, 0.1, 0.1))
Sx <- seq(1, 99, 1)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn(S)
plot(0, type="n", xlim=c(0, 100), ylim=c(0, 100), asp=1, axes=FALSE)
image(Sx, Sx, matrix(elev, len),
      col=gray(seq(0.2, 0.8, 0.01)),
#      col=rgb(0,seq(0.1,1,0.01),0,0.8),
      add=TRUE) #,
#      ann=FALSE, axes=FALSE, asp=1)
points(s, cex=0.4)
#segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
#segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
box(col=gray(0))
rect(0, 0, 100, 100, lwd=2)
plot(0, type="n", xlim=c(0, 100), ylim=c(0, 100), asp=1, axes=FALSE)
image(Sx, Sx, matrix(elev, len), xlim=c(0,1), ylim=c(0,1),
      col=gray(seq(0.2, 0.8, 0.01)), add=TRUE)
#      col=rgb(0,seq(0.1,1,0.01),0,0.8),
#      ann=FALSE, axes=FALSE, asp=1)
segments(seq(0, 100, 20), 0, seq(0, 100, 20), 100, col=gray(0))
segments(0, seq(0, 100, 20), 100, seq(0, 100, 20), col=gray(0))
box(col=gray(0))
rect(0, 0, 100, 100, lwd=2)
y <- 10
for(i in 1:nrow(n.k)) {
    text(seq(10, 100, by=20), y, labels=n.k[,i], cex=0.6)
    y <- y+20
}
par(op)
dev.off()
system("open ../figs/heteroPlots.png")




















# PART III. Analysis using custom MCMC


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
#nz <- 80
#M <- nz+nrow(y)
#yz <- array(0, c(M, ntraps, noccasions))
#yz[1:nrow(y),,] <- y # Fill data augmentation array





library(scrbook)
library(coda)

#source("../../Rpackage/scrbook/R/Ch11.R")
#source("../../Rpackage/scrbook/R/scrIPP.R")


set.seed(3434)
system.time({
fm1 <- scrIPP(y, X, M=200, 10000, xlims=c(0,100), ylims=c(0,100),
              space.cov=elev.fn,
              tune=c(0.4, 0.2, 0.3, 0.25, 9))
}) # 360s

plot(mcmc(fm1$out))
summary(mcmc(fm1$out))
summary(window(mcmc(fm1$out), start=5001))#$q

which.max(table(fm1$out[,"N"]))
HPDinterval(window(mcmc(fm1$out, start=5001)))

rejectionRate(mcmc(fm1$out))

c(N=N, n=sum(apply(y>0, 1, any)), M=M)

plot(fm1$last$S)



png("../figs/fm1p.png", width=7, height=7, units="in", res=400)
par(mfrow=c(4,2), mai=c(0.3, 0.4, 0.5, 0.2), cex.main=1.8, cex.axis=1.8)
plot(mcmc(fm1$out[,c(3,4,5)]))
dev.off()
system("open ../figs/fm1p.png")





# Analysis using secr
library(secr)


# Create a "traps" object
Xs <- data.frame(X)
colnames(Xs) <- c("x","y")
secr.traps <- read.traps(data=Xs, detector="count")

summary(secr.traps)

plot(secr.traps)
plot.default(secr.traps, xlim=c(0,100), asp=1, pch="+")

# Create a "capthist" object
secr.caps <- matrix(NA, sum(y), 5)
colnames(secr.caps) <- c("Session", "ID", "Occasion", "X", "Y")
counter <- 0
for(i in 1:nrow(y)) {
    for(j in 1:dim(y)[2]) {
        for(k in 1:dim(y)[3]) {
            y.ij <- y[i,j,k]
            if(y.ij==0)
                next
            for(v in 1:y.ij) {
                counter <- counter+1
                secr.caps[counter,] <- c(1, i, k, X[j,1], X[j,2])
            }
        }
    }
}
ch <- make.capthist(secr.caps, secr.traps, fmt="XY")
#plot(ch, tol=0.0005) # ouch

# Make mask

msk <- make.mask(secr.traps, buffer=20, spacing=1) #, nx=v)
str(msk)
summary(msk)
#plot(msk)

ssArea <- attr(msk, "area")*nrow(msk)

covariates(msk) <- data.frame(elev=apply(as.matrix(msk), 1,
                              function(x) elev.fn(x)))


system.time(secr1.0 <- secr.fit(ch, model=D~1, mask=msk)) # 78s
predict(secr1.0)
region.N(secr1.0, se.N=TRUE)

system.time(secr1.1 <- secr.fit(ch, model=D~elev, mask=msk)) # 328s
predict(secr1.1)
region.N(secr1.1, se.N=TRUE)

system.time(secr1.2 <- secr.fit(ch, model=D~xy, mask=msk)) # 129s
predict(secr1.2)
region.N(secr1.2, se.N=TRUE)

AIC(secr1.0, secr1.1, secr1.2)







fm1.s <- summary(window(mcmc(fm1$out), start=5001))
fm1.r <- cbind(fm1.s$stat[,1:2], fm1.s$quant[,c(1,5)])[c(1,2,4:6),]
fm1.r

secr1.1s <- rbind(data.matrix(predict(secr1.1)[c(3,2), 2:5]),
                  beta1=as.numeric(coef(secr1.1)[2,]),
                  region.N(secr1.1)[,1:4])[c(1,2,3,5,4),]
secr1.1s

(scrVsecr <- cbind(par=rep(c("sigma", "lam0", "beta1", "N", "EN"),each=2),
                   Method=c("MCMC", "ML"),
                   rbind(fm1.r[1,], secr1.1s[1,],
                   fm1.r[2,], secr1.1s[2,],
                   fm1.r[3,], secr1.1s[3,],
                   fm1.r[4,], secr1.1s[4,],
                   fm1.r[5,], secr1.1s[5,])))

format(scrVsecr, digits=2, nsmall=3, scientific=FALSE)

write.table(format(scrVsecr, digits=2, nsmall=3, scientific=FALSE),
            file="scrVsecr1.txt", row.names=FALSE,
            quote=FALSE, sep=" \t& ", eol=" \\\\\n ")


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













# PART IV. Discrete space


library(scrbook)
library(secr)
library(raster) # raster after secr b/c of flip
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

image(cell, cell, t(matrix(dat$CANHT, B/pix, B/pix)), ann=FALSE)

head(dat)

# Simulate IPP
set.seed(3216)
beta0 <- 3
beta1 <- 1
(EN <- sum(exp(beta0 + beta1*dat$CANHT)*pixArea))
# N <- rpois(1, EN) # 45
M <- 100
(N <- rbinom(1, M, EN/M))
dat$cp <- exp(beta0 + beta1*dat$CANHT) / EN
s.tmp <- rmultinom(1, N, dat$cp) # a single realization to be ignored later

# Trap locations
xsp <- seq(27.5, 72.5, by=5)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)


canhtMat <- t(matrix(dat$CANHT, 100/pix, 100/pix))
canht <- raster:::flip(raster(t(canhtMat)), direction="y")



png("../figs/discrete.png", width=6, height=6, units="in", res=400)
par(mai=c(0.4, 0.4, 0.2, 0.2))
image(cell, cell, canhtMat, ann=FALSE, col=gray(seq(0,1,len=50)))
points(dat[s.tmp>0,c("x","y")], cex=s.tmp[s.tmp>0])
points(X, pch="+")
box()
dev.off()
system("open ../figs/discrete.png")



# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
y <- array(NA, c(N, ntraps))

sigma <- 10  # half-normal scale parameter
lam0 <- 1   # basal encounter rate
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
y <- y[rowSums(y)>0,]

sum(y)
dim(y)
table(y)

nz <- M - nrow(y) # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps))
yz[1:nrow(y),] <- y # Fill




# Analysis using secr
library(secr)


# Create a "traps" object
Xs <- data.frame(X)
colnames(Xs) <- c("x","y")
secr.traps <- read.traps(data=Xs, detector="count")

summary(secr.traps)


plot(secr.traps)

plot.default(secr.traps, xlim=c(0, B), asp=1, pch="+")

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

msk <- make.mask(secr.traps, buffer=27.5, spacing=5) #, nx=v)
summary(msk)
#plot(msk)

ssArea <- attr(msk, "area")*nrow(msk)

covariates(msk) <- data.frame(canht=dat$CANHT[order(dat$y, dat$x)])




ch11simData <- list(ch.secr=ch, ch.jags=yz, spcov.jags=dat, spcov.secr=msk,
                    traps=X)


# save(ch11simData, file="../../Rpackage/scrbook/data/ch11simData.rda")
#promptData(ch9simData, filename="../Rpackage/scrbook/man/ch9simData.Rd")




#data(ch9simData)

#ch <- ch9simData$ch.secr
#msk <- ch9simData$spcov.secr


# SECR analysis

secr2 <- secr.fit(ch, model=D~canht, mask=msk)

coef(secr2)
beta0; beta1

region.N(secr2, se.N=TRUE)
N




# JAGS analysis

# JAGS model
sink("ippDiscrete.jag")
cat("
model{
sigma ~ dunif(0, 20)
lam0 ~ dunif(0, 5)
beta0 ~ dunif(-10, 10)
beta1 ~ dunif(-10, 10)
for(j in 1:nPix) {
  mu[j] <- exp(beta0 + beta1*CANHT[j])*pixArea
  probs[j] <- mu[j]/EN
}
EN <- sum(mu[])
psi <- EN/M
for(i in 1:M) {
  w[i] ~ dbern(psi)
  s[i] ~ dcat(probs[])
  x0g[i] <- grid[s[i],1]
  y0g[i] <- grid[s[i],2]
  for(j in 1:ntraps) {
    dist[i,j] <- sqrt(pow(x0g[i]-traps[j,1],2) + pow(y0g[i]-traps[j,2],2))
    lambda[i,j] <- lam0*exp(-dist[i,j]*dist[i,j]/(2*sigma*sigma)) * w[i]
    y[i,j] ~ dpois(lambda[i,j])
    }
  }
N <- sum(w[])
D <- N/1 # 1ha state-space
}
", fill=TRUE)
sink()



modfile <- "ippDiscrete.jag"

jags.data <- list(y=yz, CANHT=drop(dat$CANHT),
                  nPix=nrow(dat), pixArea=pixArea,
                  M=nrow(yz), ntraps=nrow(X),
                  grid=as.matrix(dat[,1:2]),
                  traps=X)
str(jags.data)

all(matrix(jags.data$CANHT, 20, byrow=TRUE) ==
    matrix(covariates(msk)$canht, 20))

init1 <- function() {
    list(sigma=runif(1, 2, 8), lam0=runif(1),
         beta0=rnorm(1, 3, .2), beta1=rnorm(1, 1, 0.3),
#         w=ifelse(rowSums(jags.data$y>0), 1, 0),
         w=rep(1, M),
         s=c(s[,"pixID"], sample(1:400, M-nrow(s))))
}
str(init1())

pars1 <- c("sigma", "lam0", "beta0", "beta1", "N", "EN")

# Obtain posterior samples.
# Compile and adapt
system.time({
    set.seed(453)
    jm <- jags.model(modfile, jags.data, init1, n.chains=2, n.adapt=1000)
    jags1 <- coda.samples(jm, pars1, n.iter=10000)
}) # 1.6hr

jags2 <- coda.samples(jm, pars1, n.iter=10000)


plot(jags1, ask=TRUE)
summary(jags1)

summary(window(jags1, start=6001))
plot(window(jags1, start=90001), ask=TRUE)


library(R2WinBUGS)
bugs1 <- bugs(jags.data, init1, pars1, modfile, n.chains=2, n.iter=2000,
              n.burnin=1000, n.thin=1, DIC=FALSE,
              bugs.directory="C:/Winbugs/WinBUGS14/")


save.image("scratch.RData")





library(parallel)
cl <- makeCluster(3)
clusterExport(cl, c("jags.data", "init1", "pars1", "modfile", "M", "s"))

fm <- clusterEvalQ(cl, {
    library(rjags)
    jm <- jags.model(modfile, jags.data, init1, n.chains=1, n.adapt=500)
    jags1 <- coda.samples(jm, pars1, n.iter=2000)
    jags1
})

str(fml <- mcmc.list(sapply(fm, mcmc)))
plot(fml, ask=TRUE)


fm2 <- clusterEvalQ(cl, {
    jags2 <- coda.samples(jm, pars1, n.iter=10000)
    jags2
})

str(fml2 <- mcmc.list(sapply(fm2, mcmc)))
plot(fml2, ask=TRUE)


stopCluster(cl)




# compare results

library(scrbook)
#example(ch9secrYjags)


plot(window(jags1, start=5001))

summary(window(jags1, start=5001))
gelman.diag(window(jags1, start=5001))


#jags.est <- summary(window(fml, start=1001))
jags.est <- summary(jags2)
jags.r <- cbind(jags.est$stat[,1:2], jags.est$quant[,c(1,5)])
jags.r

secr.est <- predict(secr2)
secr.r <- cbind(secr.est[2:3,2:5])
secr.r <- rbind(beta=as.numeric(coef(secr2)[2,]), secr.r)
secr.r <- data.matrix(rbind(region.N(secr2)[,1:4], secr.r))
secr.r

jagsVsecr <-
data.frame(Par=rep(c("$\\lambda_0$", "$\\sigma$", "$\\beta_1$", "$N$",
                   "$\\mathbb{E}[N]$"), each=2),
           Truth=rep(c(1, 10, 1, 30, 32.3), each=2),
           Software = rep(c("\\textbf{JAGS}", "\\texttt{secr}"), 5),
           rbind(jags.r[5,], secr.r[4,],
                 jags.r[6,], secr.r[5,],
                 jags.r[4,], secr.r[3,],
                 jags.r[2,], secr.r[2,],
                 jags.r[1,], secr.r[1,]))


format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE)

write.table(format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE),
            file="jagsVsecr2.txt", row.names=FALSE,
            quote=FALSE, sep=" \t& ", eol=" \\\\\n ")


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










