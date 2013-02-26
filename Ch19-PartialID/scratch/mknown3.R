
# Pr(marked) decreases from trap array centroid





# Simulate data for "m known" case

#M <- 125
#D <- 100
#A <- 1
#N <- D*1
N <- 80  # Guys in S
#m <- 20



# Dimensions of S and B (square in which marked guys are sampled)
xlimS <- c(0, 4)
ylimS <- c(0, 4)
cen <- c(mean(xlimS), mean(ylimS))

# Activity centers in S
set.seed(4830)
s <- cbind(runif(N, xlimS[1], xlimS[2]),
           runif(N, ylimS[1], ylimS[2]))

dc <- apply(s, 1, function(x) sqrt((x[1]-cen[1])^2 + (x[2]-cen[2])^2))
tau <- 0.6
w0 <- 0.7
PrMark <- exp(-dc^2/(2*tau^2)) # Change this to "net capture" prob
hist(PrMark)
#w <- rep(0, N)
#w[sample(1:N, m, prob=PrMark)] <- 1
w <- rbinom(N, 1, PrMark)
(m <- sum(w))

# Randomly sample m guys in B


co <- seq(0.8, 3.2, len=6)
X <- cbind(rep(co, each=length(co)), rep(co, times=length(co)))

plot(s, xlim=xlimS, ylim=ylimS, cex=1.3)
points(s, col="grey", pch=16, cex=1.3)
points(s[w==1,], col="blue", pch=16, cex=0.8)
points(X, pch="+", cex=2)



# Generate capture histories for the marked guys (yM)
# and latent cap histories for the unmarked guys (yU)

J <- nrow(X)
K <- 5
lam0 <- 0.5
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

sum(rowSums(yM)>0) # 16 of m guys actually detected
rowSums(y>0)

# Observed data
y <- yM[w==1,,]
dim(y)

sum(y)

rowSums(y>0)
apply(y>0, c(1,2), sum)


nU <- apply(yU, c(2,3), sum) # Trap counts, excluding marked guys
nU




# JAGS


# Augment data

nz <- 80

yz <- array(0, c(m+nz, J, K))
yz[1:m,,] <- y

rowSums(yz)



paste("yu[", (m+1):(m+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)



dat1 <- list(y=yz, nU=nU, X=X, M=nrow(yz), J=J, K=K,
             w=c(rep(1, m), rep(0, nz)),
             h=c(rep(1, m), rep(NA, nz)),
             cen=cen,
             xlimS=c(xlimS[1], xlimS[2]), ylimS=c(xlimS[1],xlimS[2]))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((m+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:m,,] <- 0

init1 <- function() list(h=c(rep(NA, m), rep(2, nz)),
                         tau=0.5,
                         yu=yui)




str(dat1)
str(init1())


pars1 <- c("N", "sigma", "lam0", "D", "ED", "EN",
           "tau", "w0", "m")

jm1 <- jags.model("mknown3.jag", dat1, init1, n.chains=1,
                  n.adapt=100)

mc1 <- coda.samples(jm1, pars1, n.iter=1000)
mc2 <- coda.samples(jm1, pars1, n.iter=5000)

#system.time({
#mc1 <- coda.samples(jm1, pars1, n.iter=10000)
#mc2 <- coda.samples(jm1, pars1, n.iter=10000)
#})

plot(mc1, ask=TRUE)
summary(mc1)

tail(mc1)

plot(mc2, ask=TRUE)
summary(mc2)


sigma
lam0
w0
tau
N
m
N / 16










