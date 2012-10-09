
# There are two ways to fit inhomogeneous point process SCR model
# (1) Change the prior from uniform to weighted
# (2) Thin uniform point process by modeling psi

# This script does method (2)








# Trap locations in [0,1] x [0,1] state-space
X <- cbind(rep(seq(.3,.7,,5), 5), rep(seq(.3,.7,,5), times=5))
X



# Spatial covariate. Elevation as function of location
elev <- function(x) x[1]+x[2]
grid <- expand.grid(seq(0,1,.1),seq(0,1,0.1))
library(lattice)
levelplot(apply(grid, 1, elev) ~ grid[,1] + grid[,2])


# Simulate data
M <- 400 # Data augmentation size. Could be anything, but we need it
s <- cbind(runif(M), runif(M)) # activity centers
beta0 <- -4
beta1 <- 2
pr <- plogis(beta0 + beta1*apply(s, 1, elev))  # thinning prob
mean(pr)*M
w <- rbinom(M, 1, pr) # Thinned guys
N <- sum(w) # Number of real guys
N

plot(s)
points(s[w==1,], col=rgb(0,0,1,0.5), pch=16)

p0 <- 0.7       # Basal capture prob
sigma <- 0.1    # Half-normal scale parameters
J <- nrow(X)    # nTraps
K <- 5          # nOcc
ym <- matrix(NA, M, J)
for(i in 1:M) {
    for(j in 1:J) {
        d2 <- (s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2
        p <- w[i]*p0*exp(-d2/(2*sigma^2))
        ym[i,j] <- rbinom(1, K, p)
    }
}


y <- ym[!apply(ym==0,1,all),]
table(ym)

M2 <- 150
yz <- rbind(y, matrix(0, M2-nrow(y), ncol(y)))
rowSums(yz)





# JAGS

library(rjags)

dat1 <- list(y=yz, X=X, J=J, K=K, M=M2)
par1 <- c("beta0", "beta1", "sigma", "p0", "N")
init1 <- function() list(w=as.integer(apply(dat1$y>0, 1, any)))

jm1 <- jags.model("modelPsi.txt", data=dat1, init1,
                  n.chains=2, n.adapt=500)
jc1 <- coda.samples(jm1, par1, n.iter=1500)
jc2 <- coda.samples(jm1, par1, n.iter=3000)


plot(jc1, ask=TRUE)
summary(jc1)

plot(jc2, ask=TRUE)
summary(jc2)




