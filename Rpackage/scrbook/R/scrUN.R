



# MCMC with Z partially observed
# Need to allow for different observation models
scrUN <- function(y, X, M, obsmod=c("pois", "bern"),
                  niters, xlims, ylims, a, b, tune=c(0.2, 0.1, 2)) {

    obsmod <- match.arg(obsmod)
    J <- nrow(y)
    K <- ncol(y)
    S <- cbind(runif(M, xlims[1], xlims[2]),
               runif(M, ylims[1], ylims[2]))
    D <- e2dist(S, X)
    sigma <- runif(1, 0.5, 1.5)
    lam0 <- runif(1, 0.7, 0.9) # This is p0 when obsmod="bern"
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    psi <- runif(1, 0.7, 1)
    w <- rbinom(M, 1, psi)

    Z <- array(0, c(M,J,K))
    up <- 0
    for(j in 1:J) {
        for(k in 1:K) {
            if(y[j,k]==0) {
                Z[,j,k] <- 0
                up <- up+1
                next
            }
            probs <- lam[,j]*w
            if(identical(obsmod, "pois")) {
                probs <- probs/sum(probs)
                Z[,j,k] <- rmultinom(1, y[j,k], probs)
            }
            else if(identical(obsmod, "bern")) {
#                Z[,j,k] <- 0
#                guys <- sample(1:M, y[j,k], prob=probs)
#                Z[guys,j,k] <- 1
                while(sum(Z[,j,k]) != y[j,k]) {
                    Z[,j,k] <- rbinom(M, 1, probs)
                    }
                up <- up+1
                cat("  z init", up, "of", J*K, "\n")
            }
        }
    }

    out <- matrix(NA,nrow=niters,ncol=4)
    colnames(out) <- c("sigma", "lam0", "psi", "N")

    cat("\nstarting values =", c(sigma, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }


        if(identical(obsmod, "pois")) {
            ll <- sum(dpois(Z, lam*w, log=TRUE))
        } else if(identical(obsmod, "bern")) {

            ll <- sum(dbinom(Z, 1, lam*w, log=TRUE))
        }

        # update sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            if(!missing(a) && !missing(b)) {
                prior <- dgamma(sigma, a, b, log=TRUE)
                prior.cand <- dgamma(sigma.cand, a, b, log=TRUE)
            } else {
                prior <- prior.cand <- 0
            }
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            # w is recycled over lam, R times
            # lam*w is recycled over Z, T times
            if(identical(obsmod, "pois")) {
                llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            } else if(identical(obsmod, "bern")) {
                llcand <- sum(dbinom(Z, 1, lam.cand*w, log=TRUE))
            }
            if(runif(1) < exp((llcand + prior.cand) -
                              (ll + prior))) {
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
       }

        # update lam0, this is p0 if obsmod="bern"
        lam0.cand <- rnorm(1, lam0, tune[2])
        test2 <- TRUE
        if(identical(obsmod, "bern"))
            test2 <- lam0.cand <= 1
        if(lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            if(identical(obsmod, "pois"))
                llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            else if(identical(obsmod, "bern"))
                llcand <- sum(dbinom(Z, 1, lam.cand*w, log=TRUE))
            if(runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

        # update w
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i])
                next
            wcand <- ifelse(w[i]==0, 1, 0)
            if(identical(obsmod, "pois")) {
                ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(Z[i,,], lam[i,]*wcand, log=TRUE))
            } else if(identical(obsmod, "bern")) {
                ll <- sum(dbinom(Z[i,,], 1, lam[i,]*w[i], log=TRUE))
                llcand <- sum(dbinom(Z[i,,], 1, lam[i,]*wcand, log=TRUE))
            }

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand
                wUps <- wUps+1
            }
        }

        # update Z
        Zups <- 0
        for(j in 1:J) {
            probs <- lam[,j]*w
            for(k in 1:K) {
                if(y[j,k]==0) {
                    Z[,j,k] <- 0
                    next
                }
#                probs <- probs/sum(probs)
                if(identical(obsmod, "pois")) {
                    probs <- probs/sum(probs)
                    Z[,j,k] <- rmultinom(1, y[j,k], probs)
                    Zups <- Zups+1
                }
                else if(identical(obsmod, "bern")) {
#                    Z[,j,k] <- 0
#                    guy <- sample(1:M, y[j,k], prob=probs)
#                    Z[guy,j,k] <- 1
                    zcand <- rbinom(M, 1, probs)
                    if(sum(zcand) != y[j,k])
                        next
                    Z[,j,k] <- zcand
                    Zups <- Zups+1
                }
            }
        }

        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update S
        Sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            Scand <- c(rnorm(1, S[i,1], tune[3]),
                       rnorm(1, S[i,2], tune[3]))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
            Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
            lam.cand <- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )

            # recycle lam over Z
            if(identical(obsmod, "pois")) {
                ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(Z[i,,], lam.cand*w[i], log=TRUE))
            } else if(identical(obsmod, "bern")) {
                ll <- sum(dbinom(Z[i,,], 1, lam[i,]*w[i], log=TRUE))
                llcand <- sum(dbinom(Z[i,,],1,lam.cand*w[i], log=TRUE))
            }

            if(runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i,] <- Scand
                lam[i,] <- lam.cand
                D[i,] <- dtmp
                Sups <- Sups+1
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     Z =", round(Zups/(J*K), 2), "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(sigma,lam0,psi,sum(w) )

    }
    return(out)
}


















