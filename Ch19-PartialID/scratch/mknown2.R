
# Andy's idea with nested squares






# Simulate data for "m known" case

#M <- 125
#D <- 100
#A <- 1
#N <- D*1
N <- 100  # Guys in S
m <- 20   # Marked guys



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


# Dimensions of S and B (square in which marked guys are sampled)
xlimS <- c(0, 4)
ylimS <- c(0, 4)
xlimB <- c(1, 3)
ylimB <- c(1, 3)

# Activity centers in S
set.seed(43)
s <- cbind(runif(N, xlimS[1], xlimS[2]),
           runif(N, ylimS[1], ylimS[2]))

# Randomly sample m guys in B
# NOTE: Not really possilbe to randomly sample activity centers,
#       but you could sample based on proportion of home range in B
#       This is just a simplification

sInB <- apply(s, 1, function(x) (x[1] > xlimB[1]) & (x[1] < xlimB[2]) &
                                (x[2] > ylimB[1]) & (x[2] < ylimB[2]))
sum(sInB) # If less than m, start over

w <- rep(0, N)
w[sample(which(sInB), m)] <- 1
w  # Indicator of marked guys


co <- seq(1.5, 2.5, len=5)
X <- cbind(rep(co, each=5), rep(co, times=5))

plot(s, xlim=xlimS, ylim=ylimS, cex=1.3)
rect(xlimB[1], ylimB[1], xlimB[2], ylimB[2])
points(s, col="grey", pch=16, cex=1.3)
points(s[w==1,], col="blue", pch=16, cex=0.8)
points(X, pch="+", cex=2)

# Andy's theta parameter:
sum(w)/sum(sInB)



# Generate capture histories for the marked guys (yM)
# and latent cap histories for the unmarked guys (yU)

J <- nrow(X)
K <- 5
lam0 <- 0.5
sigma <- 0.2

yM <- yU <- array(NA, c(N, J, K)) # Capture data
lambda <- dist <- matrix(NA, N, J)
for(i in 1:N) {
    for(j in 1:J) {
        dist[i,j] <- sqrt((s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2)
        lambda[i,j] <- lam0*exp(-dist[i,j]^2/(2*sigma^2))
        yM[i,j,] <- rpois(K, lambda[i,j] * w[i])
        yU[i,j,] <- rpois(K, lambda[i,j] * (1-w[i]))
    }
}

sum(rowSums(yM)>0) # 16 of m guys actually detected

# Observed data
y <- yM[w==1,,]
dim(y)


nU <- apply(yU, c(2,3), sum) # Trap counts, excluding marked guys
nU




# JAGS


# Augment data

nz <- 100

yz <- array(0, c(m+nz, J, K))
yz[1:m,,] <- y

rowSums(yz)



paste("yu[", (m+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)

A1 <- (xlimB[2]-xlimB[1])*(ylimB[2]-ylimB[1])
AS <- (xlimS[2]-xlimS[1])*(ylimS[2]-ylimS[1])
A2 <- AS-A1
pi <- c(A1/AS, A2/AS)


dat1 <- list(y=yz, nU=nU, X=X, M=nrow(yz), J=J, K=K,
             h=c(rep(1, m), rep(NA, nz)),
             pi=pi,
             xlimB=c(xlimB[1], xlimB[2]), ylimB=c(ylimB[1], ylimB[2]),
             xlimS=c(xlimS[1], xlimS[2]), ylimS=c(xlimS[1],xlimS[2]))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((m+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

init1 <- function() list(h=c(rep(NA, m), rep(2, nz)),
                         yu=yui)




str(dat1)
str(init1())


pars1 <- c("N", "sigma", "lam0", "D", "ED", "EN")

jm1 <- jags.model("mknown2.jag", dat1, init1, n.chains=1,
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

nz <- 100

xlim <- ylim <- c(-.1, 1.1)
(A2 <- (ylim[2]-ylim[1])*(xlim[2]-xlim[1]))
D*A2

yz <- array(0, c(nind+nz, J, K))
yz[1:nind,,] <- y

rowSums(yz)



paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)


dat2 <- list(y=yz, nU=nU, X=X, M=nind+nz, J=J, K=K,
#             z=c(rep(0, nind), rep(NA, nz)),
#             h=c(rep(1, nind), rep(NA, nz)),
#             u=c(rep(0, nind), rep(NA, nz)),
             w=c(rep(1, nind), rep(0, nz)),
             xlim=xlim, ylim=ylim)

yui <- array(0, c(dat2$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((nind+1):dat2$M, dat2$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

init1 <- function() list(#omega=0.2,
                         z=c(rep(NA, nind), rep(1, nz)),
#                         u=c(rep(NA, nind), rep(1, nz)),
                         yu=yui,
                         psi=0.3)




str(dat2)
str(init1())


jm2 <- jags.model("mknown.jag", dat2, init1, n.chains=1,
                  n.adapt=100)

mc2.1 <- coda.samples(jm2, pars1, n.iter=100)
mc2.2 <- coda.samples(jm2, pars1, n.iter=10000)


plot(mc2.1, ask=TRUE)
summary(mc2.1)

plot(mc2.2, ask=TRUE)
summary(mc2.2)



D*A2
