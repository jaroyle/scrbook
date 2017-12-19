





# spatial covariate (with mean 0)
#elev.fn.v <- function(x) x[,1]+x[,2]-1
#elev.fn <- function(x) x[1]+x[2]-1


# 2-dimensional integration over unit square
#int2d <- function(alpha, delta=0.02) {
#  z <- seq(delta/2, 1-delta/2, delta)
#  len <- length(z)
#  cell.area <- delta*delta
#  S <- cbind(rep(z, each=len), rep(z, times=len))
#  sum(exp(alpha*elev.fn(S)) * cell.area)
#  }



# Create spatial covariate

spcov <- function(B=1, pix=0.05, cor=2) {
    cell <- seq(0+pix/2, B-pix/2, pix)
    v <- length(cell)
    R <- data.frame(x=rep(cell, each=v),
                    y=rep(cell, times=v))
    D<-e2dist(R,R)
    V<-exp(-D/cor)
    Vi<-solve(V)
    cov1<-t(chol(V))%*%rnorm(nrow(R))
    R$cov1<-cov1-mean(cov1)
    cov1.fn<-function(newpt,cov1,cov1.coords=cbind(R$x,R$y),Vi){
        newpt<-matrix(newpt,ncol=2)
        k<- exp(-e2dist(newpt,cov1.coords)/cor)
        pred<-k%*%Vi%*%cov1
        as.numeric(pred)
    }
    return(list(R=R, Vi=Vi, cov1.fn=cov1.fn))
}













# MCMC. SCR model with inhomogenous point process
scrIPP <- function(y, X, M, niters, xlims, ylims, space.cov,
                   init=list(beta0=-5, beta1=0, sigma=5, lam0=1,
                             s=cbind(runif(M, xlims[1], xlims[2]),
                                     runif(M, ylims[1], ylims[2]))),
                   tune=rep(0.1, 5))
{
    if(!require(R2Cuba))
        stop("Requires the R2Cuba package")

    ydims <- dim(y)
    if(length(ydims) != 3)
        stop("y should be a 3D array of capture data")
    n <- ydims[1]
    if(n > M)
        stop("nrow(y)<M is not allowed")
    J <- ydims[2]
    K <- ydims[3]

    yold <- y
    y <- array(0L, c(M, J, K)) # created augmented dataset
    y[1:n,,] <- yold

    if((nrow(X) != J) | (ncol(X) != 2))
        stop("X should be a 2-column matrix of trap coordinates")

    # initial values
    for(i in 1:length(init))
        assign(names(init)[i], init[[i]])
    pn <- c("beta0", "beta1", "sigma", "lam0", "s")
    noi <- !(pn %in% ls())
    if(any(noi)) {
        stop("Need initial values for ", pn[noi])
    }

    D <- e2dist(s, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))

    mu <- function(x, beta0, beta1) exp(beta0 + beta1*space.cov(s=x))
    EN <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                upper=c(xlims[2], ylims[2]),
                beta0=beta0, beta1=beta1,
                flags=list(verbose=0))$value
    psi <- EN / M
    if(psi > 1)
        stop("Bad initial values for beta0 or beta1. Or M is too low")

    z <- rbinom(M, 1, psi)
    z[rowSums(y)>0] <- 1
#    z[] <- 1

    ll.y <- sum(dpois(y, lam*z, log=TRUE))

    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=7)
    colnames(out) <- c("sigma", "lam0", "beta0", "beta1", "N", "EN", "deviance")

    cat("\ninitial values =",
        c(sigma, lam0, beta0, beta1, sum(z), EN, -2*ll.y), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 ==0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("current =", out[iter-1,], "\n")
            cat("  Acceptance rates\n")
            cat("    s =", sups/M, "\n")
            cat("    z =", zUps/M, "\n")
        }


        # update sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            ll.y.cand <- sum(dpois(y, lam.cand*z, log=TRUE) )
            if(runif(1) < exp( ll.y.cand  - ll.y ) ){
                ll.y <- ll.y.cand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }

        # update lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            ll.y.cand<- sum(dpois(y, lam.cand*z, log=TRUE) )
            if(runif(1) < exp( ll.y.cand - ll.y ) ) {
                lam0<-lam0.cand
                lam<-lam.cand
                ll.y <- ll.y.cand
            }
        }

        # update z
        zUps <- 0
        seen <- apply(y>0, 1, any)
        for(i in 1:M) {
            if(seen[i])
                next
            zcand <- z
            if(z[i]==0) {
                zcand[i] <- 1
                ll.y <- 0
                ll.y.cand <- sum(dpois(y[i,,], lam[i,]*zcand[i], log=TRUE))
            } else {
                zcand[i] <- 0
                ll.y <- sum(dpois(y[i,,], lam[i,]*z[i], log=TRUE))
                ll.y.cand <- 0
            }
            ll.z <- dbinom(z[i], 1, psi, log=TRUE)
            ll.z.cand <- dbinom(zcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp((ll.y.cand+ll.z.cand) - (ll.y+ll.z))) {
                z <- zcand
                zUps <- zUps+1
            }
        }

        # update beta0
        beta0.cand <- rnorm(1, beta0, tune[3])
        EN.cand <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                         upper=c(xlims[2], ylims[2]),
                         beta0 = beta0.cand, beta1=beta1,
                         flags=list(verbose=0))$value
        psi.cand <- EN.cand/M
##        ll.beta <- sum(((beta0 + beta1*space.cov(s)) - log(EN))*z) +
##            dbinom(sum(z), M, EN/M, log=TRUE)
        ll.s <- beta0 + beta1*space.cov(s) - log(EN)
        ll.z <- dbinom(z, 1, psi, log=TRUE)
        if(EN.cand < M) {
##            ll.beta.cand <- sum(((beta0.cand + beta1*space.cov(s)) -
##                log(EN.cand))*z) + dbinom(sum(z), M, EN.cand/M, log=TRUE)
            ll.s.cand <- beta0.cand + beta1*space.cov(s) - log(EN.cand)
            ll.z.cand <- dbinom(z, 1, psi.cand, log=TRUE)
##            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
            if(runif(1) < exp((sum(ll.s.cand)+sum(ll.z.cand)) -
                              (sum(ll.s)+sum(ll.z)) ))  {
                beta0 <- beta0.cand
                EN <- EN.cand
                psi <- psi.cand
##                ll.beta <- ll.beta.cand
                ll.s <- ll.s.cand
                ll.z <- ll.z.cand
            }
        }

        # update beta1
        beta1.cand <- rnorm(1, beta1, tune[4])
        EN.cand <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                         upper=c(xlims[2], ylims[2]),
                         beta0 = beta0, beta1=beta1.cand,
                         flags=list(verbose=0))$value
        psi.cand <- EN.cand/M
        if(EN.cand < M) {
##            ll.beta.cand <- sum(((beta0 + beta1.cand*space.cov(s)) -
##                log(EN.cand))*z) + dbinom(sum(z), M, EN.cand/M, log=TRUE)
            ll.s.cand <- beta0 + beta1.cand*space.cov(s) - log(EN.cand)
            ll.z.cand <- dbinom(z, 1, psi.cand, log=TRUE)
##            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
            if(runif(1) < exp((sum(ll.s.cand)+sum(ll.z.cand)) -
                              (sum(ll.s)+sum(ll.z)) ))  {
                beta1 <- beta1.cand
                EN <- EN.cand
                psi <- psi.cand
##                ll.beta <- ll.beta.cand
                ll.s <- ll.s.cand
                ll.z <- ll.z.cand
            }
        }

##        # update psi
##        psi <- EN / M

        # update s
        sups <- 0
        for(i in 1:M) {
            scand <- c(rnorm(1, s[i,1], tune[5]),
                       rnorm(1, s[i,2], tune[5]))
            inbox <- scand[1]>=xlims[1] & scand[1]<=xlims[2] &
                     scand[2]>=ylims[1] & scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt( (scand[1] - X[,1])^2 + (scand[2] - X[,2])^2 )
            lam.cand <- lam
            lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
            if(z[i]==0) {
##                ll.s <- ll.s.cand <- 0
                ll.y <- ll.y.cand <- 0
            } else {
##                ll.s <- sum(dpois(y[i,,], lam[i,], log=TRUE) )
##                ll.s.cand <- sum(dpois(y[i,,], lam.cand[i,], log=TRUE) )
                ll.y <- sum(dpois(y[i,,], lam[i,], log=TRUE) )
                ll.y.cand <- sum(dpois(y[i,,], lam.cand[i,], log=TRUE) )
            }
            #ln(prior), denominator is constant
##            prior.s <- beta0 + beta1*space.cov(s[i,])
##            prior.s.cand <- beta0 + beta1*space.cov(scand)
            ll.s <- beta0 + beta1*space.cov(s[i,])
            ll.s.cand <- beta0 + beta1*space.cov(scand)

##           if(runif(1)< exp((ll.s.cand+prior.s.cand) - (ll.s+prior.s))) {
            if(runif(1) < exp((ll.y.cand+ll.s.cand) -
                              (ll.y+ll.s))) {
                s[i,] <- scand
                lam <- lam.cand
                D[i,] <- dtmp
                ##psi[i] <- psi.cand
                sups <- sups+1
            }
        }
        ll.y <- sum(dpois(y, lam*z, log=TRUE))
        out[iter,] <- c(sigma, lam0, beta0, beta1, sum(z), EN, -2*ll.y)
    }
    last <- list(s=s, lam=lam, z=z)
    list(out=out, last=last)
}










if(1==2) {









# Spatial covariate (with mean 0)
elev.fn <- function(s) {
    s <- matrix(s, ncol=2)        # Force s to be a matrix
    (s[,1] + s[,2] - 100) / 40.8  # Returns (standardized) "elevation"
}

mu <- function(s, beta0, beta1) exp(beta0 + beta1*elev.fn(s=s))


library(R2Cuba)
xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(100,100), beta0=0, beta1=2)

xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta0=0, beta1=2,
            flags=list("verbose"=0))


# Simulate PP using rejection sampling
##set.seed(31025)
beta0 <- -6 # intercept of intensity function
beta1 <- -1  # effect of elevation on intensity
# Next line computes integral, which is expected value of N
EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
            lower=c(0,0), upper=c(100,100))$value
EN
(N <- rpois(1, EN)) # Realized N
s <- matrix(NA, N, 2) # This matrix will hold the coordinates
elev.min <- elev.fn(c(0,0))
elev.max <- elev.fn(c(100, 100))
Q <- max(c(exp(beta0 + beta1*elev.min),
           exp(beta0 + beta1*elev.max)))
counter <- 1
while(counter <= N) {
  x.c <- runif(1, 0, 100); y.c <- runif(1, 0, 100)
  s.cand <- c(x.c,y.c)
  pr <- mu(s.cand, beta0, beta1) #/ EN
  if(runif(1) < pr/Q) {
    s[counter,] <- s.cand
    counter <- counter+1
    }
  }

plot(s)



xsp <- seq(20, 80, by=10); len <- length(xsp)
X <- cbind(rep(xsp, each=len), rep(xsp, times=len)) # traps
ntraps <- nrow(X); noccasions <- 5
y <- array(NA, c(N, ntraps, noccasions)) # capture data
sigma <- 5  # scale parameter
lam0 <- 1   # basal encounter rate
lam <- matrix(NA, N, ntraps)
set.seed(5588)
for(i in 1:N) {
    for(j in 1:ntraps) {
        # The object "s" was simulated in previous section
        distSq <- (s[i,1]-X[j,1])^2 + (s[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(noccasions, lam[i,j])
    }
}
# data augmentation
nz <- 80
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, noccasions))
yz[1:nrow(y),,] <- y # Fill data augmentation array


sum(rowSums(y)>0)



# Fit the model using MCMC
# Sample the parameters: "sigma", "lam0", "beta0", "beta1", "N", "EN"

## Not run:
##set.seed(3434)
system.time({
fm1 <- scrIPP(yz, X, M, 5000, xlims=c(0,100), ylims=c(0,100),
              space.cov=elev.fn,
              tune=c(0.45, 0.25, 0.4, 0.3, 20))
}) # 328s


library(coda)

plot(mc1 <- mcmc(fm1$out), ask=TRUE)

summary(mc1)


rejectionRate(mc1)


}
