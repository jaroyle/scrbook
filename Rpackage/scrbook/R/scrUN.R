



# MCMC with Z partially observed
# Need to allow for different observation models
scrUN <- function(y, X, Zknown, M, niters, xlims, ylims, a, b) {

    R <- nrow(y)
    T <- ncol(y)
    S <- cbind(runif(M, xlims[1], xlims[2]), runif(M, ylims[1], ylims[2]))
    D <- e2dist(S, X)
    theta <- runif(1, 1, 3)
    lam0 <- runif(1, .1, 0.6)
    lam <- lam0*exp(-(D*D)/(2*theta*theta))
    psi <- runif(1,.1,.6)
    w <- rbinom(M,1,psi)

    Z <- array(NA, c(M,R,T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if(!missing(Zknown)) {# || !is.null(Zknown)) {
        nMarked <- nrow(Zknown)
        marked[1:nMarked] <- TRUE
        Z[1:nMarked,,] <- Zknown
    }
    w[marked] <- 1
    Zdata <- !is.na(Z)
    for(r in 1:R) {
        for(t in 1:T) {
            if(y[r,t]==0) {
                Z[,r,t] <- 0
                next
            }
            unmarked <- !Zdata[,r,t]
            nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
            if(nUnknown < 0)
                browser()
            probs <- lam[,r]*w
            probs <- probs[unmarked]
            Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs)
        }
    }

    out <- matrix(NA,nrow=niters,ncol=4)
    colnames(out) <- c("sigma", "lam0", "psi", "N")

    cat("\nstarting values =", c(theta, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

        # update theta
#        theta.cand <- rnorm(1, theta, 0.1)
#        prop <- prop.back <- 0
#        if(theta.cand>0){
        tune.theta <- 0.1
        theta.cand <- rlnorm(1, log(theta), tune.theta)
        prop <- dlnorm(theta.cand, log(theta), tune.theta, log=TRUE)
        prop.back <- dlnorm(theta, log(theta.cand), tune.theta, log=TRUE)

        if(!missing(a) && !missing(b)) {
            prior <- dgamma(theta, a, b, log=TRUE)
            prior.cand <- dgamma(theta.cand, a, b, log=TRUE)
        } else {
            prior <- prior.cand <- 0
        }
            lam.cand <- lam0*exp(-(D*D)/(2*theta.cand*theta.cand))
            # w is recycled over lam, R times
            # lam*w is recycled over Z, T times
            ll <- sum(dpois(Z, lam*w, log=TRUE))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp((llcand+prior.cand+prop.back) -
                              (ll+prior+prop))) {
                ll<-llcand
                lam<-lam.cand
                theta<-theta.cand
            }
#       }

        # update lam0
#        lam0.cand <- rnorm(1, lam0, .03)
#        prop <- prop.back <- 0
#        if(lam0.cand>0) {
        tune.lam0 <- 0.25
        lam0.cand <- rlnorm(1, log(lam0), tune.lam0)
        prop <- dlnorm(lam0.cand, log(lam0), tune.lam0, log=TRUE)
        prop.back <- dlnorm(lam0, log(lam0.cand), tune.lam0, log=TRUE)
            lam.cand <- lam0.cand*exp(-(D*D)/(2*theta*theta))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp((llcand+prop.back) - (ll+prop))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
#        }

        ### update "w" here
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i] | marked[i])
                next
            wcand <- ifelse(w[i]==0, 1, 0)

            ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
            llcand <- sum(dpois(Z[i,,], lam[i,]*wcand, log=TRUE))

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand
                wUps <- wUps+1
            }
        }

        # update Z
        for(r in 1:R) {
            zip <- lam[,r]*w
            probs <- zip/sum(zip)
            for(t in 1:T) {
                if(y[r,t]==0) {
                    Z[,r,t] <- 0
                    next
                }
                unmarked <- !Zdata[,r,t]
                nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown == 0)
                    next
                Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs[unmarked])
            }
        }


        # update psi
        psi<-rbeta(1,1+sum(w),1+M-sum(w))

        # update S
        Sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            Scand <- c(rnorm(1, S[i,1], 2),
                       rnorm(1, S[i,2], 2))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
                     Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(inbox) {
                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
                lam.cand <- lam0*exp(-(dtmp*dtmp)/(2*theta*theta) )

                # recycle lam over Z
                ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(Z[i,,], lam.cand*w[i], log=TRUE))

                if(runif(1) < exp(llcand - ll)) {
                    ll <- llcand
                    S[i,] <- Scand
                    lam[i,] <- lam.cand
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(theta,lam0,psi,sum(w) )

    }

    return(out)
}


















