
# Simulate data for "m unknown" case

M <- 100
N <- 50
m <- 15

(psi <- N/M)      # Pr(z=1)
(omega <- m/N)    # Pr(marked)
(pi <- omega*psi) #

# set.seed(5459)
z <- w <- rep(0, M)
z1 <- sample(1:M, N)
z[z1] <- 1
w1 <- sample(z1, m)
w[w1] <- 1

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

yM <- array(NA, c(M, J, K))
lambda <- dist <- matrix(NA, M, J)
q <- matrix(NA, M, K)
for(i in 1:M) {
    for(j in 1:J) {
        dist[i,j] <- sqrt((s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2)
        lambda[i,j] <- lam0*exp(-dist[i,j]^2/(2*sigma^2)) * z[i]
        yM[i,j,] <- rpois(K, lambda[i,j])
    }
}

# Matrix indicating if a guy is known to be marked
qM <- apply(yM, c(1,3), sum)
qM[qM>0] <- 1
qM <- qM*w   # These guys are marked
for(i in 1:M) {
    if(all(qM[i,]==0))
        next
    first1 <- min(which(qM[i,]==1))
    qM[i,first1:K] <- 1
}
colSums(qM)


y1 <- apply(yM>0, 1, any)
sum(y1)

# This would be data if everyone was marked
yL <- yM[y1,,]  # All capture histories (observed and latent)
qL <- qM[y1,]   # Marked status (some of this is latent too)

# Sort so marked guys are first 1:nMarked rows
markedfirst <- order(rowSums(qL), decreasing=TRUE)

yL <- yL[markedfirst,,]
qL <- qL[markedfirst,]

# Observed data
y <- yL[rowSums(qL)>0,,]
q <- qL[rowSums(qL)>0,]
w <- ifelse(rowSums(q)>0, 1, 0)
n <- apply(yL, c(2,3), sum)
nU <- apply(yL[rowSums(qL)==0,,], c(2,3), sum)  # Counts of unmarked guys

dim(y)
dim(q)
nind <- dim(y)[1]

# Augment data

nz <- 100

yz <- array(NA, c(nind+nz, J, K))
yz[1:nind,,] <- y

wz <- c(w, rep(NA,nz))



# JAGS

paste("yw[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)


dat1 <- list(y=yz, w=wz, nU=nU, X=X, M=nind+nz, J=J, K=K,
             xlim=c(0, 1), ylim=c(0,1))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample(1:M, dat1$nU[j,k]),j,k] <- 1
    }
}
ywi[1:nind,,] <- NA
init1 <- function() list(z=rep(1, dat1$M),
                         yu=yui,
                         w=c(rep(NA, nind), rep(0, nz)))



str(dat1)
str(init1())



jm1 <- jags.model("munknown.jag", dat1, init1, n.chains=1, n.adapt=100)

