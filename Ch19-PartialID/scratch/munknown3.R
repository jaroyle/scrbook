
# Pr(marked) decreases from trap array centroid





# Simulate data for "m known" case

#M <- 125
#D <- 100
#A <- 1
#N <- D*1
N <- 80  # Guys in S
m <- 20   # Marked guys




# Dimensions of S and B (square in which marked guys are sampled)
xlimS <- c(0, 4)
ylimS <- c(0, 4)
cen <- c(mean(xlimS), mean(ylimS))

# Activity centers in S
set.seed(43)
s <- cbind(runif(N, xlimS[1], xlimS[2]),
           runif(N, ylimS[1], ylimS[2]))

dc <- apply(s, 1, function(x) sqrt((x[1]-cen[1])^2 + (x[2]-cen[2])^2))
tau <- 0.6
w0 <- 0.9
PrMark <- w0*exp(-dc^2/(2*tau^2))
hist(PrMark)
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

sum(rowSums(yM)>0) # 16 of m guys actually detected
rowSums(y>0)

# Observed data
#y <- yM[rowSums(yM)>0,,]
#dim(y)
y<-yM[w==1,,]

rowSums(y>0)
apply(y>0, c(1,2), sum)

(nind <- nrow(y))

nU <- apply(yU, c(2,3), sum) # Trap counts, excluding marked guys
nU




# JAGS


# Augment data

nz <- 80

#yz <- array(0, c(nind+nz, J, K))
#yz[1:nind,,] <- y

rowSums(yz)



paste("yu[", (nind+1):(nind+nz), ",j,k]",
      sep="", collapse=",")

library(rjags)

###run in r
source('c:/users/rs/dropbox/scrPIDHn.R')
inits<-function(){list(S=cbind( runif(m+nz,xlimS[1],xlimS[2]), runif(m+nz,ylimS[1],ylimS[2])),
			sigA=0.6, w0=0.9,sigma=0.25, lam0=0.2, psi=0.7 )}
out<-scrPIDHn(n=nU, X=X, y=y, M=m+nz, obsmod ="pois",niters=5000, cB=cen,h=c(rep(1,m), rep(0,nz)),
    xlims=xlimS, ylims=ylimS, inits=inits(), delta=c(0.1, 0.1, 0.5, 0.1, 0.1) ) 

yJ<-apply(y,1:2,sum)
data<-list(m=m, M=nz, J=J, K=K,y=yJ,nU=nU, h=c(rep(1,m), rep(0,nz)), xlim=xlimS, ylim=ylimS,
		cent=cen, X=X )

yui<-array(NA,c(M,J,K))
for (j in 1:J){
for(k in 1:K){
yui[,j,k]<-rmultinom(1,nU[j,k], rep(1/M,M))
}}

inits<-function(){list(s=cbind( runif(m+nz,xlimS[1],xlimS[2]), runif(m+nz,ylimS[1],ylimS[2])),
			sigA=0.6, w0=0.9,sigma=0.25, lam0=0.2, psi=0.7,z=rep(1,nz) )}
params=c("N", "psi","sigma", "lam0", "sigA", "w0")
mod<-jags.model('c:/users/rs/dropbox/mknown2Hn.jag', data, inits, n.chains=1, n.adapt=800)
out<-coda.samples(mod, params, n.iter=5000)

dat1 <- list(y=yz, nU=nU, X=X, M=nrow(yz), J=J, K=K,
             h=c(rep(1, nind), rep(NA, nz)),
             cen=cen,
             xlimS=c(xlimS[1], xlimS[2]), ylimS=c(xlimS[1],xlimS[2]))

yui <- array(0, c(dat1$M, J, K))
for(j in 1:J) {
    for(k in 1:K) {
        yui[sample((m+1):dat1$M, dat1$nU[j,k]),j,k] <- 1
    }
}
yui[1:nind,,] <- 0

init1 <- function() list(h=c(rep(NA, nind), rep(2, nz)),
                         yu=yui)




str(dat1)
str(init1())


pars1 <- c("N", "sigma", "lam0", "D", "ED", "EN",
           "tau", "w0", "m")

jm1 <- jags.model("munknown3.jag", dat1, init1, n.chains=1,
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
w0
tau
N
m
N / 16










