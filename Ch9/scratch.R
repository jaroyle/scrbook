

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
segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
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





# spatial covariate
elev.fn <- function(x) x[,1]+x[,2]-1


# 2-dimensional integration over unit square
int2d <- function(alpha, delta=0.02) {
  z <- seq(delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(alpha*elev.fn(S)) * cell.area)
  }

# Simulate PP using rejection sampling
set.seed(3005)
N <- 100
count <- 1
s <- matrix(NA, N, 2)
alpha <- 2 # parameter of interest
while(count <= 100) {
  x.c <- runif(1, 0, 1)
  y.c <- runif(1, 0, 1)
  s.cand <- cbind(x.c,y.c)
  elev.min <- elev.fn(cbind(0,0))
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






n.k <- table(cut(s[,1], seq(0, 1, 0.2)),
             cut(s[,2], seq(0, 1, 0.2)))


# plot continuous space and discrete space
png("figs/heteroPlots.png", width=5, height=2.5, units="in", res=400)
op <- par(mfrow=c(1, 2), mai=c(0.1, 0.1, 0.1, 0.1))
Sx <- seq(0.01, 0.99, 0.01)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn(S)
image(Sx, Sx, matrix(elev, len), col=rgb(0,seq(0.1,1,0.01),0,0.8),
      ann=FALSE, axes=FALSE, asp=1)
points(s, cex=0.5)
segments(seq(0, 1, 0.2), 0, seq(0, 1, 0.2), 1, col=gray(0.5))
segments(0, seq(0, 1, 0.2), 1, seq(0, 1, 0.2), col=gray(0.5))
box(col=gray(0.5))
Sx <- seq(0.1, 0.9, 0.2)
len <- length(Sx)
S <- cbind(rep(Sx, each=len), rep(Sx, times=len))
elev <- elev.fn(S)
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















# Analysis in JAGS


# 3x3 grid of traps centered
xc <- seq(0.3, 0.7, 0.2)
X <- cbind(rep(xc, each=3), rep(xc, times=3))

lam0 <- 2
sigma <- 0.15

set.seed(32235)
y <- matrix(NA, nrow(s), nrow(X))
for(i in 1:nrow(s)) {
    for(j in 1:nrow(X)) {
        dist <- sqrt((s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2)
        lam <- lam0 * exp(-dist^2/(2*sigma*sigma))
        y[i,j] <- rpois(1, lam)
    }
}

sum(y)

nz <- 50
yz <- matrix(0, nrow(y)+nz, ncol(y))
yz[1:nrow(y),] <- y


sink("ippDiscrete.txt")
cat("
model{
sigma ~ dunif(0, 1)
lam0 ~ dunif(0, 5)
beta ~ dnorm(0,0.1)
psi ~ dbeta(1,1)

for(j in 1:nPix) {
  theta[j] <- exp(beta*elevation[j])
}

for(j in 1:nPix) {
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





library(rjags)

modfile <- "ippDiscrete.txt"
file.show(modfile)
dat <- list(y=yz, elevation=elev, nPix=prod(dim(n.k)),
            M=nrow(yz), ntraps=nrow(X), Sgrid=S, grid=X)
init <- function() {
    list(sigma=runif(1), lam0=runif(1), beta=rnorm(1),
         s=
         w=c(rep(1,100), rep(0,nz)), psi=1)
}
pars <- c("sigma", "lam0", "beta", "N")


jm <- jags.model(modfile, dat, init, n.chains=2, n.adapt=500)
jc <- coda.samples(jm, pars, n.iter=2000)

plot(jc)

summary(window(jc, start=1001))
