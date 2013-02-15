
# Simulate data for "m unknown" case
M <- 200
N <- 100
m <- 50

(psi <- N/M)      # Pr(z=1)
(omega <- m/N)    # Pr(marked)
(pi <- omega*psi) #

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

plot(s, xlim=c(0,1), ylim=c(0,1))
points(s[z==1,], col="grey", pch=16)
points(s[q==1,], col="blue", pch=16)
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
qM <- qM*w
for(i in 1:M) {
    if(all(qM[i,]==0))
        next
    first1 <- min(which(qM[i,]==1))
    qM[i,first1:K] <- 1
}
colSums(qM)


y1 <- apply(yM>0, 1, any)
sum(y1)

# Observed data
y <- yM[y1,,]
q <- qM[y1,]
