scrPID<-function (n, X, y, M, obsmod = c("pois", "bern"),niters, npics,
    xlims, ylims, a, b, inits, delta ) 
{
    obsmod <- match.arg(obsmod)

    J <- nrow(n)
    K <- ncol(n)
    S <- inits$S
    D <- e2dist(S, X)
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    psi <- inits$psi
    z <- rbinom(M, 1, psi)
    Y <- array(NA, c(M, J, K))
    nMarked <- nrow(y)
    marked <- rep(FALSE, M)
        marked[1:nMarked] <- TRUE
        Y[1:nMarked, , ] <- y
    z[marked] <- 1
    Ydata <- !is.na(Y)
    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            unmarked <- !Ydata[, j, k]
            nUnknown <- n[j, k] 
            probs <- lam[, j] * z
            probs <- probs[unmarked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois")) 
                Y[unmarked, j, k] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[unmarked, j, k] <- 0
                guys <- sample(which(unmarked), nUnknown, prob = probs)
                Y[guys, j, k] <- 1
            }
        }
    }


	cr<-rep(1,M) 
	if(missing(npics)){
	crat<-1} else{
	crat<-npics[1]/npics[2]}
        cr[marked]<-crat

    out <- matrix(NA, nrow = niters, ncol = 5)
    colnames(out) <- c("sigma", "lam0","c", "psi", "N")
    cat("\nstarting values =", c(sigma, lam0, crat, psi, sum(z)), "\n\n")
    for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), 
                "\n")
            cat("   current =", out[iter - 1, ], "\n")
        }


        if(identical(obsmod, "pois")) {
            ll <- sum(dpois(Y, lam*cr * z, log = TRUE))
       } else if(identical(obsmod, "bern")) {
           ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) 
       }

	if(!missing(npics)){
	crat<-rbeta(1, 1+npics[1], 1+npics[2]-npics[1]) #npics[1]=identified, npics[2]=all marked
	cr[marked]<-crat
 	}

        sigma.cand <- rnorm(1, sigma, delta[1])
        if (sigma.cand > 0) {
            if (!missing(a) && !missing(b)) {
                prior <- dgamma(sigma, a, b, log = TRUE)
                prior.cand <- dgamma(sigma.cand, a, b, log = TRUE)
            }
            else {
                prior <- prior.cand <- 0
            }
            lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand * sigma.cand))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam*cr * z, log = TRUE))
                llcand <- sum(dpois(Y, lam.cand*cr * z, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) 
                llcand <- sum(dbinom(Y, 1, lam.cand*cr * z, log = TRUE)) 
            }
            if (runif(1) < exp((llcand + prior.cand) - (ll + prior))) {
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }
        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern")) 
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma * sigma))
            if (identical(obsmod, "pois")) 
                llcand <- sum(dpois(Y, lam.cand*cr * z, log = TRUE))
            else if (identical(obsmod, "bern")) 
                llcand <- sum(dbinom(Y, 1, lam.cand*cr * z, log = TRUE)) 
            if (runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }
        zUps <- 0
        seen <- apply(Y > 0, 1, any)
        for (i in 1:M) {
            if (seen[i] | marked[i]) 
                next
            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand, 
                  log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i], 
                  log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] * 
                  zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
            if (runif(1) < exp((llcand + prior.cand) - (ll + 
                prior))) {
                z[i] <- zcand
                zUps <- zUps + 1
            }
        }
        for (j in 1:J) {
            zip <- lam[, j] * z
            for (k in 1:K) {
                if (n[j, k] == 0) {
                  Y[!marked, j, k] <- 0
                  next
                }
                unmarked <- !Ydata[, j, k]
                nUnknown <- n[j, k]
                if (nUnknown == 0) 
                  next
                probs <- zip[unmarked]
                probs <- probs/sum(probs)
                if (identical(obsmod, "pois")) 
                  Y[unmarked, j, k] <- rmultinom(1, nUnknown, probs)
                else if (identical(obsmod, "bern")) {
                  Y[unmarked, j, k] <- 0
                  guy <- sample(which(unmarked), nUnknown, prob = probs)
                  Y[guy, j, k] <- 1
                }
            }
        }
               

        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + (M-sum(marked)) - sum(z[!marked]))


        Sups <- 0
        for (i in 1:M) {
            Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1, S[i, 2], delta[3]))
            inbox <- Scand[1] >= xlims[1] & Scand[1] <= xlims[2] & 
                Scand[2] >= ylims[1] & Scand[2] <= ylims[2]
            if (!inbox) 
                next
            dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2*sigma*sigma))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] *cr[i] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam.cand *cr[i]* z[i], 
                  log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] *cr[i]* z[i], 
                  log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam.cand *cr[i]* z[i], log = TRUE))
            }
            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i, ] <- Scand
                lam[i, ] <- lam.cand
                D[i, ] <- dtmp
                Sups <- Sups + 1
            }
        }
        if (iter%%100 == 0) {
            cat("   Acceptance rates\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, sum(z))
    }
    return(out)
}
