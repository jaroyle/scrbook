





# spatial covariate (with mean 0)
elev.fn.v <- function(x) x[,1]+x[,2]-1
elev.fn <- function(x) x[1]+x[2]-1


# 2-dimensional integration over unit square
#int2d <- function(alpha, delta=0.02) {
#  z <- seq(delta/2, 1-delta/2, delta)
#  len <- length(z)
#  cell.area <- delta*delta
#  S <- cbind(rep(z, each=len), rep(z, times=len))
#  sum(exp(alpha*elev.fn(S)) * cell.area)
#  }



# Create spatial covariate

spcov <- function(B=1, pix=0.05) {
    cell <- seq(0+pix/2, B-pix/2, pix)
    v <- length(cell)
    R <- data.frame(x=rep(cell, each=v),
                    y=rep(cell, times=v))
    D<-e2dist(R,R)
    V<-exp(-D/2)
    Vi<-solve(V)
    cov1<-t(chol(V))%*%rnorm(nrow(R))
    R$cov1<-cov1-mean(cov1)
    cov1.fn<-function(newpt,cov1,cov1.coords=cbind(R$x,R$y),Vi){
        newpt<-matrix(newpt,ncol=2)
        k<- exp(-e2dist(newpt,cov1.coords)/2)
        pred<-k%*%Vi%*%cov1
        as.numeric(pred)
    }
    return(list(R=R, Vi=Vi, cov1.fn=cov1.fn))
}













# MCMC. SCR model with inhomogenous point process
scrIPP <- function(Z, X, M, niters, xlims, ylims, tune=rep(0.1, 4))
{

    if(!require(R2Cuba))
        stop("Requires the R2Cuba package")

    Zdims <- dim(Z)
    R <- Zdims[2]
    T <- Zdims[3]

    # initial values
    S <- cbind(runif(M,xlims[1],xlims[2]),runif(M,ylims[1],ylims[2]))
    D <- e2dist(S, X)
    sigma <-runif(1, .3, .6)
    lam0 <- runif(1, 4, 6)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))

    psi <- runif(1, 0.4, 0.6)
    beta1 <- rnorm(1, 0)

    w <- rbinom(M, 1, psi)
    w[rowSums(Z)>0] <- 1

    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=5)
    colnames(out) <- c("sigma", "lam0", "psi", "beta1", "N")

    mu <- function(x, beta) exp(beta*elev.fn(x=x))

    cat("\ninitial values =", c(sigma, lam0, psi, beta1, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 ==0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("current =", out[iter-1,], "\n")
            cat("  Acceptance rates\n")
            cat("    S =", Sups/M, "\n")
            cat("    w =", wUps/M, "\n")
        }

        ll<- sum(dpois(Z, lam*w, log=TRUE))

        # update sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE) )
            if(runif(1)<exp( llcand  - ll ) ){
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }

        # update lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE) )
            if(runif(1) < exp( llcand - ll ) ) {
                lam0<-lam0.cand
                lam<-lam.cand
                ll <- llcand
            }
        }

        # update w
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i])
                next
            wcand<-w
            if(w[i]==0) {
                wcand[i] <- 1
                ll.w <- 0
                ll.w.cand <- sum(dpois(Z[i,,], lam[i,]*wcand[i], log=TRUE))
            } else {
                wcand[i] <- 0
                ll.w <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                ll.w.cand <- 0
            }
            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp((ll.w.cand+prior.cand) - (ll.w+prior))) {
                w <- wcand
                wUps <- wUps+1
            }
        }

        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update beta1
#        D1 <- int2d(beta1, delta=.05)
        sink(file="NUL")
        D1 <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                    upper=c(xlims[2], ylims[2]), beta=beta1,
                    flags=list(verbose=0))$value
        beta1.cand <- rnorm(1, beta1, tune[3])
#        D1.cand <- int2d(beta1.cand, delta=0.05)
        D1.cand <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                    upper=c(xlims[2], ylims[2]), beta=beta1.cand,
                    flags=list(verbose=0))$value
        sink()
        ll.beta1 <- sum(  beta1*elev.fn.v(S) - log(D1) )
        ll.beta1.cand <- sum( beta1.cand*elev.fn.v(S) - log(D1.cand) )
        if(runif(1) < exp(ll.beta1.cand - ll.beta1) )  {
          beta1<-beta1.cand
          }


        # update S
        Sups <- 0
        for(i in 1:M) {
#            Scand <- matrix(c(rnorm(1, S[i,1], tune[4]),
#                              rnorm(1, S[i,2], tune[4])), nrow=1)
            Scand <- c(rnorm(1, S[i,1], tune[4]),
                              rnorm(1, S[i,2], tune[4]))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
                     Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
            lam.cand <- lam
            lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
            if(w[i]==0)
                ll.S <- ll.S.cand <- 0
            else {
                ll.S <- sum(dpois(Z[i,,], lam[i,], log=TRUE) )
                ll.S.cand <- sum(dpois(Z[i,,], lam.cand[i,], log=TRUE) )
            }
            #ln(prior), denominator is constant
            prior.S <- beta1*elev.fn(S[i,]) # - log(D1)
            prior.S.cand <- beta1*elev.fn(Scand) # - log(D1)

           if(runif(1)< exp((ll.S.cand+prior.S.cand) - (ll.S+prior.S))) {
                S[i,] <- Scand
                lam <- lam.cand
                D[i,] <- dtmp
                ##psi[i] <- psi.cand
                Sups <- Sups+1
            }
        }


        out[iter,] <- c(sigma,lam0,psi,beta1,sum(w) )
    }
    last <- list(S=S, lam=lam, w=w)
    list(out=out, last=last)
}









