
tr <- seq(0.3, 0.7, length=10)
X <- cbind(rep(tr, each=length(tr)),
           rep(tr, times=length(tr)))    # trap coords
set.seed(10)
xlim <- c(0, 1); ylim <- c(0, 1)         # S is the unit square
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1]) # area of S
mu <- 50                                 # density (animals/unit area)
(N <- rpois(1, mu*A))                    # Generate N=50 as Poisson deviate
s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
plot(X, xlim=xlim, ylim=ylim, pch="+")
points(s, col=gray(0.5), pch=16)

sigma <- 0.04
lam0 <- 0.3
J <- nrow(X)
K <- 5
y <- array(NA, c(N, J, K))
for(j in 1:J) {
    dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2] - s[,2])^2)
    lambda <- lam0*exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        y[,j,k] <- rpois(N, lambda)
    }
}
table(y)

table(apply(apply(y, c(1,2), sum)>0, 1, sum))


n <- apply(y, c(2,3), sum)
dimnames(n) <- list(paste("trap", 1:J, sep=""),
                    paste("night", 1:K, sep=""))
n[1:5,]


plot(X, cex=rowSums(n), asp=1, xlim=c(0,1))

# Analyze in JAGS


library(rjags)
dat1 <- list(n=n, X=X, J=J, K=K, M=200, xlim=xlim, ylim=ylim)
init1 <- function() {
    yi <- array(0, c(dat1$M, dat1$J, dat1$K))
    for(j in 1:dat1$J) {
        for(k in 1:dat1$K) {
            yi[sample(1:dat1$M, dat1$n[j,k]),j,k] <- 1
        }
    }
    list(sigma=runif(1, 1, 2), lam0=runif(1),
         y=yi, z=rep(1, dat1$M))
}
pars1 <- c("lam0", "sigma", "N", "mu")

system.time({
jm <- jags.model("SCmod1.jag", data=dat1, inits=init1, n.chain=1,
                 n.adapt=1000)
samples1 <- coda.samples(jm, pars1, n.iter=10000)
})

samples2 <- coda.samples(jm, pars1, n.iter=5000)

plot(samples1)







# Analyze using scrUN()

library(scrbook)
library(coda)

source("../../Rpackage/scrbook/R/scrUN.R")

fm1 <- scrUN(n=n, X=X, M=200, niter=50000, xlims=xlim, ylims=ylim,
             inits=list(lam0=0.3, sigma=0.01),
             updateY=TRUE,
             tune=c(0.004, 0.07, 0.3))

mc1 <- mcmc(fm1)
plot(mc1)
summary(window(mc1, start=10001))

rejectionRate(mc1)
rejectionRate(window(mc1, start=10001))

save(mc1, file="scrUNmc1.gzip")
