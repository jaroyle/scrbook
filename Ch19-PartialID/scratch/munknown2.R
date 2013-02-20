
# Andy's idea with nested squares






# Simulate data for "m known" case

#M <- 125
#D <- 100
#A <- 1
#N <- D*1
N <- 100  # Guys in S
m <- 15   # Marked guys




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


co <- seq(0.8, 3.2, len=6)
X <- cbind(rep(co, each=length(co)), rep(co, times=length(co)))

plot(s, xlim=xlimS, ylim=ylimS, cex=1.3)
rect(xlimB[1], ylimB[1], xlimB[2], ylimB[2])
points(s, col="grey", pch=16, cex=1.3)
points(s[w==1,], col="blue", pch=16, cex=0.8)
points(X, pch="+", cex=2)

# Andy's theta parameter:
sum(w)/sum(sInB)


A1 <- (xlimB[2]-xlimB[1])*(ylimB[2]-ylimB[1])
AS <- (xlimS[2]-xlimS[1])*(ylimS[2]-ylimS[1])
A2 <- AS-A1
(pi <- c(A1/AS, A2/AS))



# Generate capture histories for the marked guys (yM)
# and latent cap histories for the unmarked guys (yU)

J <- nrow(X)
K <- 5
lam0 <- 0.2
sigma <- 0.25

set.seed(5060)
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

sum(rowSums(yM)>0) # 12 of 15 guys actually detected
rowSums(y>0)

# Observed data
y <- yM[rowSums(yM)>0,,]
dim(y)

rowSums(y>0)
apply(y>0, c(1,2), sum)

(nind <- nrow(y))

nU <- apply(yU, c(2,3), sum) # Trap counts, excluding marked guys
nU




# JAGS


# Augment data

nz <- 100

yz <- array(0, c(nind+nz, J, K))
yz[1:nind,,] <- y

rowSums(yz)



paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)



dat1 <- list(y=yz, nU=nU, X=X, M=nrow(yz), J=J, K=K,
#             m=m,
             h=c(rep(1, nind), rep(NA, nz)),
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

init1 <- function() list(h=c(rep(NA, nind), rep(2, nz)), # start guys in B
                         # need to start all guys in B!
                         s=cbind(runif(dat1$M, xlimB[1], ylimB[2]),
                                 runif(dat1$M, ylimB[1], ylimB[2])),
                         yu=yui)




str(dat1)
str(init1())


pars1 <- c("N", "sigma", "lam0", "D", "ED", "EN", "uInB", "uOutB",
           "theta", "m")

jm1 <- jags.model("munknown2.jag", dat1, init1, n.chains=1,
                  n.adapt=100)

mc1 <- coda.samples(jm1, pars1, n.iter=100)
mc2 <- coda.samples(jm1, pars1, n.iter=500)


plot(mc1, ask=TRUE)
summary(mc1)

tail(mc1)

plot(mc2, ask=TRUE)
summary(mc2)


sigma
lam0
N
m
sum(sInB) - m # uInB
N - sum(sInB) # uOutB










