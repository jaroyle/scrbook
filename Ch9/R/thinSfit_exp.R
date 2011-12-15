
# spatial covariate
cov <- function(x, y) x+y


# 2-dimensional integration over [-1, 1] square
int2d <- function(beta, delta=0.02) {
  z <- seq(-1+delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  zz <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(beta*cov(zz[,1], zz[,2])) * cell.area)
  }


# Assess influence of grid resolution on accuracy
int2d(1, delta=0.1)
int2d(1, delta=0.05)
int2d(1, delta=0.005)
int2d(1, delta=0.001)
int2d(1, delta=0.0005)

beta <- 2

exp(beta*cov(1,1)) / int2d(beta)
exp(beta*cov(0,0)) / int2d(beta)
exp(beta*cov(-1,-1)) / int2d(beta)



# Simulate IPP using rejection sampling

N <- 100
count <- 1
S <- matrix(NA, N, 2)
beta <- -1
while(count <= 100) {
  x.c <- runif(1, -1, 1)
  y.c <- runif(1, -1, 1)
  pr <- exp(beta*cov(x.c, y.c)) / int2d(beta)
  M <- max(c(exp(beta*cov(1,1)) / int2d(beta),
             exp(beta*cov(-1,-1)) / int2d(beta)))
  if(runif(1) < pr/M) {
    S[count,] <- c(x.c, y.c)
    count <- count+1
    cat("accept\n")
    }
  else
    cat("  reject\n")
  }

plot(S, xlim=c(-1, 1), ylim=c(-1, 1))





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


# MCCM code

source("thinSmcmc_exp.R")
ls()

fm1 <- scrIPP(yz, X, M, 3000, xlims=c(-1,1), ylims=c(-1,1),
            tune=c(0.002, 0.1, 0.25, 0.07) )


plot(mcmc(fm1$out))

rejectionRate(mcmc(fm1$out))






















# Simulation study


source("thinSmcmc.R")
ls()



plot(function(x) plogis(0 + 5*x), -1, 1)


simIPP <- function(N=100, sigma=0.05, lam0=5, #beta0=0,
                   beta1=2, T=5, M=250, xlims=c(0,1), ylims=c(0,1)) {
    S <- matrix(NA, N, 2)
    count <- 1
#    elev <- function(x, y) x+y
    elev <- function(x, y) x-mean(xlims) + y-mean(ylims)
    while(count<=N) {
        sx <- runif(1, xlims[1], xlims[2])
        sy <- runif(1, ylims[1], ylims[2])
        psi <- plogis(beta1*elev(sx,sy))
        if(runif(1) < psi) {
            S[count,] <- c(sx, sy)
            count <- count+1
        }
    }
    X <- cbind(rep(seq(0.2, 0.8, by=0.08), each=8),
               rep(seq(0.2, 0.8, by=0.08), times=8))
    ntraps <- nrow(X)
    y <- array(NA, c(N, ntraps, T))
    yz <- array(0, c(M, ntraps, T))
    lam <- matrix(NA, N, ntraps)
    for(i in 1:N) {
        for(j in 1:ntraps) {
            distSq <- (S[i,1]-X[j,1])^2 + (S[i,2] - X[j,2])^2
            lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
            y[i,j,] <- rpois(T, lam[i,j])
        }
    }
    yz[1:nrow(y),,] <- y
    return(list(z=yz, X=X, S=S, M=M, xlims=xlims, ylims=ylims))
}

set.seed(34503)
sim1 <- simIPP(beta0=0, beta1=2, M=150)
str(sim1)

plot(sim1$X, xlim=c(0,1), ylim=c(0,1))
points(sim1$S, col=4, pch=16)

sum(rowSums(sim1$z)>0)


fm150 <- scrIPP(sim1$z, sim1$X, 6000, sim1$xlims, sim1$ylims,
              tune=c(0.001, 0.2, 0.01, 0.9, 0.02))


plot(mcmc(fm150$out))

hist(fm150$out[1001:6000,"N"], breaks=50:150)

rejectionRate(mcmc(fm150$out))




set.seed(34503)
sim2 <- simIPP(beta0=0, beta1=2, M=300)
str(sim1)

plot(sim2$X, xlim=c(0,1), ylim=c(0,1))
points(sim2$S, col=4, pch=16)

sum(rowSums(sim2$z)>0)

fm300.2 <- scrIPP(sim2$z, sim2$X, 6000, sim2$xlims, sim2$ylims,
              tune=c(0.001, 0.2, 0.01, 0.9, 0.02))


plot(mcmc(fm300.2$out))



hist(fm150$out[5001:6000,"N"], breaks=50:150, border="white",
     col=gray(0.9), xlab="N", main="")
hist(fm300.2$out[5001:6000,"N"], breaks=50:150,
#     col=rgb(0,0,0,0.1),
     add=TRUE)



source("thinSmcmc.R")
ls()


nsim <- 20
niter <- 27000
simout <- array(NA, c(niter, 5, nsim))
colnames(simout) <- c("sigma", "lam0", "beta0", "beta1", "N")
for(i in 1:nsim) {
    cat("\nsim", i, "\n")
    sim.i <- simIPP(beta1=2, M=300)
    fm.i <- scrIPP(sim.i$z, sim.i$X, niter,
                   sim.i$xlims, sim.i$ylims,
                   tune=c(0.001, 0.2, 0.4, 0.9, 0.02))
    simout[,,i] <- fm.i$out
}



stats1 <- t(colMeans(simout[2001:niter,,]))

hist(stats1[,1])
hist(stats1[,2])
hist(stats1[,3])
hist(stats1[,4])
hist(stats1[,5])

plot(mcmc(simout[,,1]))
plot(mcmc(simout[,,2]))
plot(mcmc(simout[,,3]))
plot(mcmc(simout[,,4]))


summary(mcmc(simout[,,1]))

crosscorr(mcmc(simout[1001:niter,,1]))
crosscorr(mcmc(simout[1001:niter,,2]))
crosscorr(mcmc(simout[1001:niter,,3]))
crosscorr(mcmc(simout[1001:niter,,4]))

rejectionRate(mcmc(simout[1001:niter,,1]))

























# Andy's stuff



# Create spatial covariate
B <- 1
v <- 20
R <- data.frame(x=rep(seq(0, B, length=v), each=v),
                y=rep(seq(0, B, length=v), times=v))
elev <- function(x, y) (x-.5)+(y-.5) # Elevation is a function of x,y
D<-e2dist1(R,R)
V<-exp(-D/2)
Vi<-solve(V)

cov1<-t(chol(V))%*%rnorm(nrow(R))
image(matrix(cov1,20,20))
R$elev <- apply(R[,1:2], 1, function(x) elev(x[1], x[2]))
R$cov1<-cov1-mean(cov1)
cov1.fn<-function(newpt,cov1,cov1.coords=cbind(R$x,R$y),Vi){
newpt<-matrix(newpt,ncol=2)
k<- exp(-e2dist1(newpt,cov1.coords)/2)
pred<-k%*%Vi%*%cov1
as.numeric(pred)
}
