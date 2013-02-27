



# MCMC with y partially observed
# Need to allow for different observation models
scrUN <- function(n, X, M, obsmod=c("pois", "bern"),
                  niters, xlims, ylims, a, b, tune=c(0.2, 0.1, 2),
                  zGibbs=FALSE) {

    obsmod <- match.arg(obsmod)
    J <- nrow(n)
    K <- ncol(n)
    s <- cbind(runif(M, xlims[1], xlims[2]),
               runif(M, ylims[1], ylims[2]))
    D <- e2dist(s, X)
    sigma <- runif(1, 0, 0.2)
    lam0 <- runif(1, 0.1, 0.8) # This is p0 when obsmod="bern"
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    psi <- runif(1, 0.9, 0.99)
    w <- rbinom(M, 1, psi)
    N0 <- sum(w)
    nmx <- max(n, na.rm=TRUE)
    if(N0 < nmx) {
        is0 <- which(w==0)
        w[sample(is0, nmx-N0)] <- 1
    }

    y <- array(0, c(M,J,K))
    up <- 0
    for(j in 1:J) {
        for(k in 1:K) {
            if(n[j,k]==0) {
                y[,j,k] <- 0
                up <- up+1
                next
            }
            probs <- lam[,j]*w
            if(identical(obsmod, "pois")) {
                probs <- probs/sum(probs)
                y[,j,k] <- rmultinom(1, n[j,k], probs)
            }
            else if(identical(obsmod, "bern")) {
                y[,j,k] <- 0
                guys <- sample(1:M, n[j,k], prob=probs/sum(probs))
                y[guys,j,k] <- 1
            }
        }
    }

    out <- matrix(NA,nrow=niters,ncol=4)
    nought <- ifelse(obsmod=="pois", "lam0", "p0")
    colnames(out) <- c("sigma", nought, "psi", "N")

    cat("\nstarting values =", c(sigma, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }


        if(identical(obsmod, "pois")) {
            ll <- sum(dpois(y, lam*w, log=TRUE))
        } else if(identical(obsmod, "bern")) {
            ll <- sum(dbinom(y, 1, lam*w, log=TRUE))
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
            # lam*w is recycled over y, T times
            if(identical(obsmod, "pois")) {
                llcand <- sum(dpois(y, lam.cand*w, log=TRUE))
            } else if(identical(obsmod, "bern")) {
                llcand <- sum(dbinom(y, 1, lam.cand*w, log=TRUE))
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
                llcand <- sum(dpois(y, lam.cand*w, log=TRUE))
            else if(identical(obsmod, "bern"))
                llcand <- sum(dbinom(y, 1, lam.cand*w, log=TRUE))
            if(runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

        # update w
        wUps <- 0
        seen <- apply(y>0, 1, any)
        for(i in 1:M) {
            if(seen[i])
                next
            wcand <- ifelse(w[i]==0, 1, 0)
            if(identical(obsmod, "pois")) {
                ll <- sum(dpois(y[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(y[i,,], lam[i,]*wcand, log=TRUE))
            } else if(identical(obsmod, "bern")) {
                ll <- sum(dbinom(y[i,,], 1, lam[i,]*w[i], log=TRUE))
                llcand <- sum(dbinom(y[i,,], 1, lam[i,]*wcand, log=TRUE))
            }

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand
                wUps <- wUps+1
            }
        }

        # update y
        yups <- 0
        for(j in 1:J) {
            probs <- lam[,j]*w
            for(k in 1:K) {
                if(n[j,k]==0) {
                    y[,j,k] <- 0
#                    yups <- yups+1
                    next
                }
#                probs <- probs/sum(probs)
                if(identical(obsmod, "pois")) {
                    probs <- probs/sum(probs)
                    y[,j,k] <- rmultinom(1, n[j,k], probs)
                    yups <- yups+1
                }
                else if(identical(obsmod, "bern")) {
                    if(zGibbs) {
                        y[,j,k] <- 0
                        probs <- probs/sum(probs)
                        guy <- sample(1:M, n[j,k], prob=probs)
                        y[guy,j,k] <- 1
                        yups <- 1
                    } else {
#                    zcand <- rbinom(M, 1, probs)
#                    if(sum(zcand) != n[j,k])
#                        next
#                    y[,j,k] <- zcand
#                    yups <- yups+1

#                    z.cand <- sample(y[,j,k]) # random permutation
                        # alternative: just move one 1
                        z.cand <- y[,j,k]
                        w1 <- w==1
                        z1w1 <- which(z.cand==1 & w1)
                        to0 <- sample(z1w1, 1)
                        z.cand[to0] <- 0
                        z0w1 <- which(z.cand==0 & w1)
                        to1 <- sample(z0w1, 1)
                        if(identical(to0, to1))
                            next
                        z.cand[to1] <- 1
                        prior <- sum(dbinom(y[,j,k], 1, probs, log=TRUE))
                        prior.cand <- sum(dbinom(z.cand, 1, probs,
                                                 log=TRUE))
                        # no likelihood contribution
                        if(runif(1) < exp(prior.cand - prior)) {
                            y[,j,k] <- z.cand
                            yups <- yups+1
                        }
                    }
                }
            }
        }

        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update s
        sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            scand <- c(rnorm(1, s[i,1], tune[3]),
                       rnorm(1, s[i,2], tune[3]))
            inbox <- scand[1]>=xlims[1] & scand[1]<=xlims[2] &
            scand[2]>=ylims[1] & scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt((scand[1] - X[,1])^2 + (scand[2] - X[,2])^2)
            lam.cand <- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )

            # recycle lam over y
            if(identical(obsmod, "pois")) {
                ll <- sum(dpois(y[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(y[i,,], lam.cand*w[i], log=TRUE))
            } else if(identical(obsmod, "bern")) {
                ll <- sum(dbinom(y[i,,], 1, lam[i,]*w[i], log=TRUE))
                llcand <- sum(dbinom(y[i,,],1,lam.cand*w[i], log=TRUE))
            }

            if(runif(1) < exp(llcand - ll)) {
                ll <- llcand
                s[i,] <- scand
                lam[i,] <- lam.cand
                D[i,] <- dtmp
                sups <- sups+1
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     y =", round(yups/(J*K), 2), "\n")
            cat("     s =", sups/M, "\n")
        }

        out[iter,] <- c(sigma,lam0,psi,sum(w) )

    }
    return(out)
}


















