











X <- cbind(rep(seq(.3,.7,,5), 5), rep(seq(.3,.7,,5), times=5))
X

x <- seq(-2, 2, .01)
plogis(-2 + 2*.1)/sum(plogis(-2 + 2*x))


M <- 200
elev <- function(x) x[1]+x[2]
grid <- expand.grid(seq(0,1,.1),seq(0,1,0.1))
library(lattice)
levelplot(apply(grid, 1, elev) ~ grid[,1] + grid[,2])

sM <- cbind(runif(M), runif(M))
pr <- plogis(-2 + 2*apply(sM, 1, elev))
w <- rbinom(M, 1, pr)

plot(sM)
points(sM[w==1,], fg="blue", pch=16)

p0 <- 0.7
sigma <- 0.2
J <- nrow(X)
K <- 5
y <- matrix(NA, M, J)
for(i in 1:M) {
    for(j in 1:J) {
        d2 <- (sM[i,1]-X[j,1])^2 + (sM[i,2]-X[j,2])^2
        p <- w[i]*p0*exp(-d2/(2*sigma^2))
        y[i,j] <- rbinom(1, K, p)
    }
}

table(y)






# JAGS

library(rjags)

dat1 <- list(y=y, X=X, J=J, K=K, M=M)
par1 <- c("beta0", "beta1", "sigma", "p0", "N")
init1 <- function() list(w=rep(1, dat1$M))
jm1 <- jags.model("modelPsi.txt", data=dat1, init1,
                  n.chains=2, n.adapt=500)
jc1 <- coda.samples(jm1, par1, n.iter=1500)


plot(jc1, ask=TRUE)
summary(jc1)




