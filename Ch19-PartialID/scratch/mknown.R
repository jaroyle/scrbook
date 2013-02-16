







# Simulate data for "m known" case

#M <- 100
#psi <- 0.5
#EN <- M*psi
#omega <- 0.3
#Em <- M*omega

M <- 125
D <- 50
A <- 1
N <- D*1
m <- 40



# Which state is a guy in
set.seed(545)
z <- w <- rep(0, M)
z1 <- sample(1:M, N)
z[z1] <- 1            # guy is real
w1 <- sample(z1, m)
w[w1] <- 1            # real and marked
u <- z*(1-w)          # real and unmarked

cbind(z, w, u)
colSums(cbind(z, w, u))

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
y <- yM[w==1,,]
dim(y)
nind <- nrow(y)

nU <- apply(yU, c(2,3), sum)





# JAGS


# Augment data

nz <- 50

yz <- array(0, c(nind+nz, J, K))
yz[1:nind,,] <- y

rowSums(yz)



paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)


dat1 <- list(y=yz, nU=nU, X=X, M=nind+nz, J=J, K=K,
             z=c(rep(0, nind), rep(NA, nz)),
             w=c(rep(1, nind), rep(0, nz)),
             xlim=c(0, 1), ylim=c(0,1))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((nind+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

init1 <- function() list(omega=0.2,
#                         z=rep(1, dat1$M),
                         z=c(rep(NA, nind), rep(1, nz)),
#                         u=c(rep(NA, nind), rep(1, nz)),
                         yu=yui,
                         psi=0.3)




str(dat1)
str(init1())


pars1 <- c("N", "m", "U", "sigma", "lam0", "D")

jm1 <- jags.model("mknown.jag", dat1, init1, n.chains=1,
                  n.adapt=100)

mc1 <- coda.samples(jm1, pars1, n.iter=100)
mc2 <- coda.samples(jm1, pars1, n.iter=500)


plot(mc1, ask=TRUE)
summary(mc1)

tail(mc1)

plot(mc2, ask=TRUE)
summary(mc2)


N
m
sum(u)












# Try again with larger state-space

nz <- 200

A2 <- 4  # [-0.5,1.5]x[-0.5,1.5] square
D*A2

yz <- array(0, c(nind+nz, J, K))
yz[1:nind,,] <- y

rowSums(yz)



paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)


dat1 <- list(y=yz, nU=nU, X=X, M=nind+nz, J=J, K=K,
#             z=c(rep(1, nind), rep(NA, nz)),
             u=c(rep(0, nind), rep(NA, nz)),
             w=c(rep(1, nind), rep(0, nz)),
             xlim=c(-0.5, 1.5), ylim=c(-0.5,1.5))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((nind+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

init1 <- function() list(omega=0.2,
#                         z=c(rep(NA, nind), rep(1, nz)),
                         u=c(rep(NA, nind), rep(1, nz)),
                         yu=yui,
                         psi=0.3)




str(dat1)
str(init1())


pars1 <- c("N", "m", "U", "sigma", "lam0", "D")

jm1 <- jags.model("mknown.jag", dat1, init1, n.chains=1,
                  n.adapt=100)

mc1 <- coda.samples(jm1, pars1, n.iter=100)
mc2 <- coda.samples(jm1, pars1, n.iter=500)


plot(mc1, ask=TRUE)
summary(mc1)

plot(mc2, ask=TRUE)
summary(mc2)

