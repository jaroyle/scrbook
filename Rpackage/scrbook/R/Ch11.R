





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
scrIPP <- function(Z, X, M, niters, xlims, ylims, tune=rep(0.1, 5),
                   beta.init=c(5,2))
{
    if(!require(R2Cuba))
        stop("Requires the R2Cuba package")

    Zdims <- dim(Z)
    R <- Zdims[2]
    T <- Zdims[3]

    elev.fn <- function(x) x[1]+x[2]-200
    elev.fn.v <- function(x) x[,1]+x[,2]-200

    # initial values
    S <- cbind(runif(M,xlims[1],xlims[2]),runif(M,ylims[1],ylims[2]))
    D <- e2dist(S, X)
    sigma <-runif(1, 20, 30) #.3, .6)
    lam0 <- runif(1, 4, 6)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))

    beta0 <- beta.init[1]
    beta1 <- beta.init[2]

    mu <- function(x, beta0, beta1) exp(beta0 + beta1*elev.fn(x=x))
    EN <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                upper=c(xlims[2], ylims[2]),
                beta0=beta0, beta1=beta1,
                flags=list(verbose=0))$value
    psi <- EN / M
    if(psi > 1)
        stop("Bad initial values for beta0 or beta1. Or M is too low")

    w <- rbinom(M, 1, psi)
    w[rowSums(Z)>0] <- 1
#    w[] <- 1

    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=5)
    colnames(out) <- c("sigma", "lam0", "beta0", "beta1", "N")

    cat("\ninitial values =", c(sigma, lam0, beta0, beta1, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 ==0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("current =", out[iter-1,], "\n")
            cat("  Acceptance rates\n")
            cat("    S =", Sups/M, "\n")
            cat("    w =", wUps/M, "\n")
            cat("    EN =", EN, "\n")
        }

        ll <- sum(dpois(Z, lam*w, log=TRUE))

        # update sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE) )
            if(runif(1) < exp( llcand  - ll ) ){
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
            wcand <- w
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

        # update beta0
        beta0.cand <- rnorm(1, beta0, tune[3])
        EN.cand <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                         upper=c(xlims[2], ylims[2]),
                         beta0 = beta0.cand, beta1=beta1,
                         flags=list(verbose=0))$value
        ll.beta <- sum((beta0 + beta1*elev.fn.v(S))*w) - EN
#        ll.beta <- sum(beta0 + beta1*elev.fn.v(S)) - EN
        if(EN.cand < M) {
#            warning("uh-oh")
#        }
            ll.beta.cand <- sum((beta0.cand + beta1*elev.fn.v(S))*w) -
                EN.cand
#            ll.beta.cand <- sum(beta0.cand + beta1*elev.fn.v(S)) -
#                EN.cand
            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
                beta0 <- beta0.cand
                EN <- EN.cand
                ll.beta <- ll.beta.cand
            }
        }

        # update beta1
        beta1.cand <- rnorm(1, beta1, tune[4])
        EN.cand <- cuhre(2, 1, mu, lower=c(xlims[1], ylims[1]),
                         upper=c(xlims[2], ylims[2]),
                         beta0 = beta0, beta1=beta1.cand,
                         flags=list(verbose=0))$value
        if(EN.cand < M) {
            ll.beta.cand <- sum((beta0 + beta1.cand*elev.fn.v(S))*w) -
                EN.cand
#            ll.beta.cand <- sum(beta0 + beta1.cand*elev.fn.v(S)) -
#                EN.cand
            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
                beta1 <- beta1.cand
                EN <- EN.cand
                ll.beta <- ll.beta.cand
            }
        }

        # update psi
        psi <- EN / M

        # update S
        Sups <- 0
        for(i in 1:M) {
#            Scand <- matrix(c(rnorm(1, S[i,1], tune[4]),
#                              rnorm(1, S[i,2], tune[4])), nrow=1)
            Scand <- c(rnorm(1, S[i,1], tune[5]),
                       rnorm(1, S[i,2], tune[5]))
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
            prior.S <- beta0 + beta1*elev.fn(S[i,])  - log(EN)
            prior.S.cand <- beta0 + beta1*elev.fn(Scand)  - log(EN)

           if(runif(1)< exp((ll.S.cand+prior.S.cand) - (ll.S+prior.S))) {
                S[i,] <- Scand
                lam <- lam.cand
                D[i,] <- dtmp
                ##psi[i] <- psi.cand
                Sups <- Sups+1
            }
        }
        out[iter,] <- c(sigma,lam0,beta0,beta1,sum(w) )
    }
    last <- list(S=S, lam=lam, w=w)
    list(out=out, last=last)
}









