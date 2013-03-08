
# SCR model for unmarked populations
scrUN <- function(n, X, M, #obsmod=c("pois", "bern"),
                  updateY=TRUE,
                  niters, xlims, ylims,
                  inits=NULL, priors=NULL,
                  tune=c(0.2, 0.1, 2)) {

#    obsmod <- match.arg(obsmod)
    obsmod <- "pois"
    J <- nrow(n)
    K <- ncol(n)
    if(!is.null(inits)) {
        for(i in 1:length(inits))
            assign(names(inits)[i], inits[[i]])
        pn <- c("sigma", "lam0", "s", "z")
        if(updateY)
            pn <- c(pn, "y")
        noi <- !(pn %in% ls())
        if(any(noi)) {
            cat("\nWarning: Generating initial values for:", pn[noi], "\n")
        }
    }
    if(!("s" %in% ls()))
        s <- cbind(runif(M, xlims[1], xlims[2]),
                   runif(M, ylims[1], ylims[2]))
    if(!("sigma" %in% ls()))
        sigma <- runif(1, min(abs(diff(X[,1]))), max(abs(diff(X[,1]))))
    if(!("lam0" %in% ls()))
        lam0 <- runif(1, 0.1, 0.5) # This is p0 when obsmod="bern"
    if(!("psi" %in% ls()))
       psi <- runif(1, 0.9, 0.99)

    dist <- e2dist(s, X)
    lam <- lam0*exp(-(dist*dist)/(2*sigma*sigma))

    if(!("y" %in% ls())) {
        y <- array(0L, c(M, J, K))
        for(j in 1:J) {
            for(k in 1:K) {
                y[sample(1:M, n[j,k], replace=TRUE,
                         lam[,j]/sum(lam[,j])), j,k] <- 1
            }
        }
    }
    if(!("z" %in% ls()))
       z <- ifelse(rowSums(y)>0, 1L, 0L)

    if(!all(dim(s) == c(M, 2)))
        stop("The dimensions of 2 should be ", c(M, 2))
    if(!all(dim(y) == c(M, J, K)) & updateY)
        stop("The dimensions of y should be ", c(M, J, K))
    if(length(z) != M)
        stop("length(z) should be ", M)


    out <- matrix(NA, nrow=niters, ncol=4)
    nought <- ifelse(obsmod=="pois", "lam0", "p0")
    colnames(out) <- c("sigma", nought, "psi", "N")

    cat("\nstarting values =", c(sigma, lam0, psi, sum(z)), "\n\n")

    if(updateY) {
        for(iter in 1:niters) {

            if(iter %% 100 == 0) {
                cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
                cat("   current =", out[iter-1,], "\n")
            }

            ll <- sum(dpois(y, lam*z, log=TRUE))

            # update sigma
            sigma.cand <- rnorm(1, sigma, tune[1])
            if(sigma.cand > 0) {
                prior <- prior.cand <- 0
                lam.cand <- lam0*exp(-(dist*dist) /
                                     (2*sigma.cand*sigma.cand))
                llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
                if(runif(1) < exp((llcand + prior.cand) - (ll + prior))) {
                    ll <- llcand
                    lam <- lam.cand
                    sigma <- sigma.cand
                }
            }

            # update lam0
            lam0.cand <- rnorm(1, lam0, tune[2])
            if(lam0.cand >= 0) {
                lam.cand <- lam0.cand*exp(-(dist*dist)/(2*sigma*sigma))
                llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
                if(runif(1) < exp((llcand) - (ll))) {
                    ll <- llcand
                    lam0 <- lam0.cand
                    lam <- lam.cand
                }
            }

            # update z
            zUps <- 0
            seen <- rowSums(y) > 0 #(y>0, 1, any)
            for(i in 1:M) {
                if(seen[i])
                    next
                zcand <- z
                zcand[i] <- ifelse(z[i]==0, 1, 0)
                ll <- sum(dpois(y[i,,], lam[i,]*z[i], log=TRUE))
                llcand <- sum(dpois(y[i,,], lam[i,]*zcand[i], log=TRUE))
                prior <- dbinom(z[i], 1, psi, log=TRUE)
                prior.cand <- dbinom(zcand[i], 1, psi, log=TRUE)
                if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                    z[i] <- zcand[i]
                    zUps <- zUps+1
                }
            }

            # update y
            for(j in 1:J) {
                probs <- lam[,j]*z
                for(k in 1:K) {
                    if(n[j,k] == 0) {
                        y[,j,k] <- 0
                        next
                    }
                    probs <- probs/sum(probs)
                    y[,j,k] <- rmultinom(1, n[j,k], probs)
                }
            }

            # update psi
            psi <- rbeta(1, 1+sum(z), 1+M-sum(z))

            # update s
            sups <- 0
            for(i in 1:M) {   # note this is "M" in general
                scand <- c(rnorm(1, s[i,1], tune[3]),
                           rnorm(1, s[i,2], tune[3]))
                inbox <- (scand[1] >= xlims[1] & scand[1] <= xlims[2]) &
                         (scand[2] >= ylims[1] & scand[2] <= ylims[2])
                if(!inbox)
                    next
                dtmp <- sqrt((scand[1] - X[,1])^2 + (scand[2] - X[,2])^2)
                lam.cand <- lam
                lam.cand[i,] <- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )
                ll <- sum(dpois(y[i,,], lam[i,]*z[i], log=TRUE))
                llcand <- sum(dpois(y[i,,], lam.cand[i,]*z[i], log=TRUE))
                if(runif(1) < exp(llcand - ll)) {
                    ll <- llcand
                    s[i,] <- scand
                    lam[i,] <- lam.cand[i,]
                    dist[i,] <- dtmp
                    sups <- sups+1
                }
            }
            if(iter %% 100 == 0) {
                cat("   Acceptance rates\n")
                cat("     z =", zUps/M, "\n")
                cat("     s =", sups/M, "\n")
            }
            out[iter,] <- c(sigma,lam0,psi,sum(z) )
        }

    } else if(!updateY) {

        for(iter in 1:niters) {

            if(iter %% 100 == 0) {
                cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
                cat("   current =", out[iter-1,], "\n")
            }

            ll <- sum(dpois(n, colSums(lam*z), log=TRUE))

            # update sigma
            sigma.cand <- rnorm(1, sigma, tune[1])
            if(sigma.cand > 0) {
                prior <- prior.cand <- 0
                lam.cand <- lam0*exp(-(dist*dist) /
                                     (2*sigma.cand*sigma.cand))
                llcand <- sum(dpois(n, colSums(lam.cand*z), log=TRUE))
                if(runif(1) < exp((llcand + prior.cand) - (ll + prior))) {
                    ll <- llcand
                    lam <- lam.cand
                    sigma <- sigma.cand
                }
            }

            # update lam0
            lam0.cand <- rnorm(1, lam0, tune[2])
            if(lam0.cand >= 0) {
                lam.cand <- lam0.cand*exp(-(dist*dist)/(2*sigma*sigma))
                llcand <- sum(dpois(n, colSums(lam.cand*z), log=TRUE))
                if(runif(1) < exp((llcand) - (ll))) {
                    ll <- llcand
                    lam0 <- lam0.cand
                    lam <- lam.cand
                }
            }

            # update z
            zUps <- 0
            for(i in 1:M) {
                zcand <- z
                zcand[i] <- ifelse(z[i]==0, 1, 0)
#                ll <- sum(dpois(n, lam[i,]*z[i], log=TRUE))
                llcand <- sum(dpois(n, colSums(lam*zcand), log=TRUE))
                prior <- dbinom(z[i], 1, psi, log=TRUE)
                prior.cand <- dbinom(zcand[i], 1, psi, log=TRUE)
                if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                    z[i] <- zcand[i]
                    zUps <- zUps+1
                    ll <- llcand
                }
            }

            # update psi
            psi <- rbeta(1, 1+sum(z), 1+M-sum(z))

            # update s
            sups <- 0
            for(i in 1:M) {   # note this is "M" in general
                scand <- c(rnorm(1, s[i,1], tune[3]),
                           rnorm(1, s[i,2], tune[3]))
                inbox <- (scand[1] >= xlims[1] & scand[1] <= xlims[2]) &
                         (scand[2] >= ylims[1] & scand[2] <= ylims[2])
                if(!inbox)
                    next
                dtmp <- sqrt((scand[1] - X[,1])^2 + (scand[2] - X[,2])^2)
                lam.cand <- lam
                lam.cand[i,] <- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )
                llcand <- sum(dpois(n, colSums(lam.cand*z), log=TRUE))
                if(runif(1) < exp(llcand - ll)) {
                    ll <- llcand
                    s[i,] <- scand
                    lam[i,] <- lam.cand[i,]
                    dist[i,] <- dtmp
                    sups <- sups+1
                }
            }
            if(iter %% 100 == 0) {
                cat("   Acceptance rates\n")
                cat("     z =", zUps/M, "\n")
                cat("     s =", sups/M, "\n")
            }
            out[iter,] <- c(sigma, lam0, psi, sum(z) )
        }
    }
    last <- list(sigma=sigma, lam0=lam0, psi=psi, z=z, s=s)
    if(updateY)
        last$y <- y
    ret <- list(sims=out, last=last)
    return(ret)
#    return(out)
}


















