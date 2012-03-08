
# Simulate quadrat count data
# guys have activity centers and they move
# Only 2 parameters! tau = scale parameter of bivariate normal move model
#                      p = detection prob

library(raster)


# Dimensions of the state-space
xlims <- ylims <- c(0, 100)

# As a raster
r <- raster(xmn=xlims[1], xmx=xlims[2], ymn=ylims[1], ymx=ylims[2])
res(r) <- c(5,5)
values(r) <- 1
r
plot(r)


# Traps
buff <- 25
xc <- seq(xlims[1]+buff, xlims[2]-buff, by=5)
lx <- length(xc)
X <- cbind(rep(xc, times=lx), rep(xc, each=lx))
dim(X)


# guys have activity centers (s) and coordinates at time k (u)
N <- 50
K <- 5
tau <- 2
s <- cbind(runif(N, xlims[1], ylims[2]), runif(N, ylims[1], ylims[2]))
u <- array(NA, c(N, 2, K))
#plot(r)
plot(X, pch="+", xlim=xlims, ylim=ylims, asp=1)
points(s, pch=16, col=4)
for(k in 1:K) {
    u[,,k] <- cbind(rnorm(N, s[,1], tau), rnorm(N, s[,2], tau))
    points(u[,,k], pch=16, col=3, cex=0.5)
}

# tabulate the number of guys in each trap pixel (quadrat)
# True abundance at trap pixels during each occasion
Nj <- matrix(0, nrow(X), K)
rownames(Nj) <- cellFromXY(r, X)
for(k in 1:K) {
    cells <- cellFromXY(r, u[,,k])
    counts <- table(cells)
    counts.in <- counts[names(counts) %in% rownames(Nj)]
    Nj[names(counts.in),k] <- counts.in
}

Nj


# Imperfect detection
p <- 0.5
n <- Nj
n[] <- rbinom(prod(dim(Nj)), Nj, p)
n



source("sampler.R")


fm1 <- scrQUAD(n, X, 50, r, 10000, xlims, ylims, tune=c(0.1, 0.5, 2))

mc1 <- mcmc(fm1)
plot(mc1)
rejectionRate(mc1)
