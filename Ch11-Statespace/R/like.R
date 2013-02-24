
choose(10, 5)


allH <- function(N, n) {
    nH <- choose(N, n)
    h1 <- combn(N, n)
    H <- matrix(0L, N, nH)
    for(h in 1:nH) {
        H[h1[,h], h] <- 1
    }
    return(H)
}

allH(4, 2)








# Simulate data

N <- 10
sigma <- 0.1
p0 <- 0.8

co <- seq(0.3, 0.7, length=5)
X <- cbind(rep(co, each=length(co)), rep(co, times=length(co)))
J <- nrow(X)
K <- 5

s <- cbind(runif(N), runif(N))

y <- array(NA, c(N, J, K))
for(i in 1:N) {
    for(j in 1:J) {
        dist <- sqrt((X[j,1]-s[i,1])^2 + (X[j,2]-s[i,2])^2)
        p <- p0*exp(-dist^2/(2*sigma^2))
        y[i,j,] <- rbinom(K, 1, p)
    }
}

n <- apply(y, 2:3, sum)
n

plot(s)




cond.like <- function(s, ...) {
    dist <- sqrt((X[j,1]-s[1])^2 + (X[j,2]-s[2])^2)
    p <- p0 * exp(-dist^2/(2*sigma^2))
    dbinom(H[i,h], 1, p)/A
}

int2d <- function(delta=0.05) {
    sco <- seq(delta/2, 1-delta/2, by=delta)
    smat <- cbind(rep(sco, each=length(sco)),
                  rep(sco, times=length(sco)))
    dist <- apply(smat, 1, function(s) sqrt((X[j,1]-s[1])^2 +
                                            (X[j,2]-s[2])^2))
    p <- p0*exp(-dist^2/(2*sigma^2))
    mean(dbinom(H[i,h]), 1, p)


nll <- function(pars, n) {
    sigma <- pars[1]
    p0 <- pars[2]
    density <- pars[3]
    A <- 1 # area of unit square
    J <- nrow(n)
    K <- ncol(n)
    L <- matrix(NA, J, K)
    for(j in 1:J) {
        for(k in 1:K) {
            LcN <- rep(1, Nmax)
            for(N in max(n):Nmax) {
                H <- allH(N, n[j,k])
                nH <- ncol(H)
                LcH <- rep(NA, nrow(H), nH)
                for(h in 1:ncol(H)) {
                    for(i in 1:nrow(H)) {
                        LcH[i,h] <- cuhre(2, 1, cond.like,
#                                         p0=p0, sigma=sigma,
#                                          y.ijk=H[i,h], A=A,
                                         lower=c(0,0), upper=c(1,1))$value
                    }
                }
                LcN[N+1] <- sum(log(LcH))
            }
            L[j,k] <- sum(exp(LcN + dpois(0:N, density*A, log=TRUE)))
        }
    }
    -sum(log(L))
}


Nmax <- 15
fm <- optim(c(0.1,0.5,5), nll, n=n, control=list(trace=TRUE, REPORT=1))


debugonce(nll)
