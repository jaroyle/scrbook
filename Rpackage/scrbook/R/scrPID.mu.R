scrPIDu<-function (n, X, y, M, mmax, obsmod = c("pois", "bern"),niters, npics,
    xlims, ylims, inits, delta ) 
{
    obsmod <- match.arg(obsmod)

    J <- nrow(n)
    K <- ncol(n)
    S <- inits$S
    D <- e2dist(S, X)
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    nObs <- nrow(y)
    Y<-array(0,c(M+mmax, J,K))
        Y[1:nObs, , ] <- y  #augmented SCR matrix
    marked<-rep(FALSE, M+mmax)
    marked[1:mmax]<-TRUE
    psi <- inits$psi
    psim<-inits$psim
    z <- rbinom(M+mmax, 1, psi)
    z[1:nObs]<-1


    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            nUnknown <- n[j, k] 
            probs <- lam[!marked, j] * z[!marked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois")) 
                Y[!marked, j, k] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[!marked, j, k] <- 0
                guys <- sample((mmax+1):(M+mmax), nUnknown, prob = probs)
                Y[guys, j, k] <- 1
           }
        }
    }


	cr<-rep(1,M+mmax) 
	if(missing(npics)){
	crat<-1} else{
	crat<-npics[1]/npics[2]}
        cr[marked]<-crat

    out <- matrix(NA, nrow = niters, ncol = 7)
    colnames(out) <- c("sigma", "lam0","c", "psi","psim", "m", "N")
    cat("\nstarting values =", c(sigma, lam0, crat, psi,psim, sum(z[marked]), sum(z)), "\n\n")
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
 
            lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand * sigma.cand))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam*cr * z, log = TRUE))
                llcand <- sum(dpois(Y, lam.cand*cr * z, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) 
                llcand <- sum(dbinom(Y, 1, lam.cand*cr * z, log = TRUE)) 
            }
            if (runif(1) < exp(llcand  - ll )) {
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

##update z for marked
        zUpsm<-zUps <- 0

        for (i in (nObs+1):mmax) {

            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                llz <- sum(dpois(Y[i, , ], lam[i, ] * cr[i]* z[i], log = TRUE))
                llcandz <- sum(dpois(Y[i, , ], lam[i, ] *cr[i]* zcand, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                llz <- sum(dbinom(Y[i, , ], 1, lam[i, ] *cr[i]* z[i], log = TRUE))
                llcandz <- sum(dbinom(Y[i, , ], 1, lam[i, ] *cr[i]* zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psim, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psim, log = TRUE)
            if (runif(1) < exp((llcandz + prior.cand) - (llz + prior))) {
                z[i] <- zcand
                zUpsm <- zUpsm + 1
            }
        }

	seen<-apply(Y>0,1,any)
        for (i in (mmax+1):(M+mmax)) {
		if(seen[i]) next

            zcand <- ifelse(z[i] == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand, log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i], log = TRUE))
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] * zcand, log = TRUE))
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log = TRUE)
		rat<-(llcand + prior.cand) - (ll + prior)
		#if(is.na(rat)) print(c(i,ll,llcand, prior, prior.cand))
            if (runif(1) < exp(rat)) {
                z[i] <- zcand
                zUps <- zUps + 1
            }
        }

        for (j in 1:J) {
            zip <- lam[!marked, j] * z[!marked]
            for (k in 1:K) {
                if (n[j, k] == 0) {
                  Y[!marked, j, k] <- 0
                  next
                }
                nUnknown <- n[j, k] 

                probs <- zip/sum(zip)
                if (identical(obsmod, "pois")) 
                  Y[!marked, j, k] <- rmultinom(1, nUnknown, probs)
                else if (identical(obsmod, "bern")) {
                  Y[!marked, j, k] <- 0
                  guy <- sample((mmax+1):(M+mmax), nUnknown, prob = probs)
                  Y[guy, j, k] <- 1
                }
            }
        }
               

        psim <- rbeta(1, 1 + sum(z[marked]), 1 + mmax - sum(z[marked]))

        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + M - sum(z[!marked]))


        Sups <- 0
        for (i in 1:(M+mmax)) {
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
            cat("     zm =", zUpsm/mmax, "\n")
            cat("     S =", Sups/(M+mmax), "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, psim, sum(z[marked]), sum(z))
    }
    return(out)
}
