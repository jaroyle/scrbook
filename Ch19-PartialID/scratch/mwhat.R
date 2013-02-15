


# Simulate data for "m unknown" case

M <- 100
psi <- 0.5
EN <- M*psi
omega <- 0.3
Em <- M*omega


# You can be in one of three states: marked, unmarked, or faker
# Here are categorical probs:
pi <- c(Em/M, (EN-Em)/M, (M-EN)/M)
sum(pi)


# Which state is a guy in
set.seed(545)
h <- rmultinom(M, 1, pi)
w <- h[1,]
u <- h[2,]
z <- w+u

(N <- sum(z))   # pop size
(m <- sum(w))   # marked guys
(U <- sum(u))   # unmarked guys

# Constraints. Should be true
N == m + U
M == m + U + sum(h[3,])

# set.seed(5459)
#z <- w <- rep(0, M)
#z1 <- sample(1:M, N)
#z[z1] <- 1
#w1 <- sample(z1, m)
#w[w1] <- 1

cbind(z, w, u)

#z <- rbinom(M, 1, psi)
#q <- rbinom(M, 1, omega*z)
#(N <- sum(z))
#(m <- sum(q))

# Activity centers
s <- cbind(runif(M), runif(M))

co <- seq(0.3, 0.7, len=5)
X <- cbind(rep(co, each=5), rep(co, times=5))

plot(s, xlim=c(0,1), ylim=c(0,1), cex=1.3)
points(s[z==1,], col="grey", pch=16, cex=1.3)
points(s[w==1,], col="blue", pch=16, cex=0.8)
points(X, pch="+", cex=2)


J <- nrow(X)
K <- 5

lam0 <- 0.5
sigma <- 0.1


yM <- yU <- array(NA, c(M, J, K)) # Capture data
lambda <- dist <- matrix(NA, M, J)
for(i in 1:M) {
    for(j in 1:J) {
        dist[i,j] <- sqrt((s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2)
        lambda[i,j] <- lam0*exp(-dist[i,j]^2/(2*sigma^2))
        yM[i,j,] <- rpois(K, lambda[i,j] * w[i])
        yU[i,j,] <- rpois(K, lambda[i,j] * u[i])
    }
}


# Observed data
y <- yM[rowSums(yM)>0,,]
dim(y)
nind <- nrow(y)

nU <- apply(yU, c(2,3), sum)

# Augment data

nz <- 100

yz <- array(0, c(nind+nz, J, K))
yz[1:nind,,] <- y




# JAGS

paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)


dat1 <- list(y=yz, nU=nU, X=X, M=nind+nz, J=J, K=K,
             xlim=c(0, 1), ylim=c(0,1))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((nind+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

wi <- ifelse(rowSums(dat1$y)>0, 1, 0)
ui <- 1-wi
hi <- cbind(wi, ui, 0)

init1 <- function() list(omega=0.2,
#                         h=hi,
                         H=apply(hi==1, 1, which),
                         yu=yui,
                         psi=0.3)
#z=rep(1, dat1$M),
#                         u=c(),
#                         yu=yui,
#                         w=c(rep(1, nind), rep(0,dat1$M-nind)))
#                         w=c(rep(NA, nind), rep(0, nz)))




str(dat1)
str(init1())


pars1 <- c("N", "m", "U", "sigma", "lam0")

jm1 <- jags.model("munknown.jag", dat1, init1, n.chains=1,
                  n.adapt=100)

mc1 <- coda.samples(jm1, pars1, n.iter=500)
mc2 <- coda.samples(jm1, pars1, n.iter=500)


plot(mc1, ask=TRUE)
summary(mc1)

plot(mc2, ask=TRUE)
summary(mc2)


N
m

