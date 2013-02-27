
tr <- seq(0.3, 0.7, length=5)
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

sigma <- 0.1
lam0 <- 0.5
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

n <- apply(y, c(2,3), sum)
dimnames(n) <- list(paste("trap", 1:J, sep=""),
                    paste("night", 1:K, sep=""))
n


# Analyze in JAGS


