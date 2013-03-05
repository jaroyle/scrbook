scrPID.tel<-function (n, X, y, M, locs,telID, obsmod = c("pois", "bern"),npics, niters, 
    xlims, ylims, a, b, inits, delta ) 
{
library(mvtnorm)
    obsmod <- match.arg(obsmod)

    R <- nrow(n)
    T <- ncol(n)
    S <- inits$S
    Sin<-t(sapply(locs, colMeans))
    S[telID,]<-Sin
    D <- e2dist(S, X)
    ntot<-length(locs)
    tel<-rep(FALSE, M)
    tel[telID]<- TRUE #does ind. have telemetry tag?
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))
    lam[lam==0] <- 1e-300
    psi <- inits$psi
    z <- rbinom(M, 1, psi)
    Y <- array(NA, c(M, R, T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if (!missing(y)) {
        nMarked <- nrow(y)
        marked[1:nMarked] <- TRUE
        Y[1:nMarked, , ] <- y
    }
    z[marked] <- 1
    Ydata <- !is.na(Y)
    for (r in 1:R) {
        for (t in 1:T) {
            if (n[r, t] == 0) {
                Y[!marked, r, t] <- 0
                next
            }
            unmarked <- !Ydata[, r, t]
            nUnknown <- n[r, t] 
            if (nUnknown < 0) 
                browser()
            probs <- lam[, r] * z
            probs <- probs[unmarked]
            probs <- probs/sum(probs)
            if (identical(obsmod, "pois")) 
                Y[unmarked, r, t] <- rmultinom(1, nUnknown, probs)
            else if (identical(obsmod, "bern")) {
                Y[unmarked, r, t] <- 0
                guys <- sample(which(unmarked), nUnknown, prob = probs)
                Y[guys, r, t] <- 1
            }
        }
    }

	cr<-rep(1,M) 

	if(missing(npics)){
	crat<-1} else{
	crat<-npics[1]/npics[2]}

	cr[marked]<-crat

    out <- matrix(NA, nrow = niters, ncol = 5)
    colnames(out) <- c("sigma", "lam0", "c", "psi", "N")
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


########################################################################################################################
######## telemetry model to estimate sigma #############################################################################


### update sigma

sigma.cand <- rnorm(1, sigma, delta[1])
if (sigma.cand > 0) {

lls<-lls.cand<-rep(NA, ntot) 

for (x in 1:ntot) {
	lls[x]<-sum(dmvnorm(x=locs[[x]],mean=c(S[telID[x],1],S[telID[x],2]), sigma=cbind(c(sigma^2,0), c(0,sigma^2)), log=T))   
	lls.cand[x]<-sum(dmvnorm(x=locs[[x]],mean=c(S[telID[x],1],S[telID[x],2]), sigma=cbind(c(sigma.cand^2,0), c(0,sigma.cand^2)), log=T))   
	}
   if(runif(1) < exp( sum(lls.cand)  - sum(lls) ) ){
    sigma<-sigma.cand
    lam <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
					}

}

########################################################################################################################
########################################################################################################################

	if(!missing(npics)){
	crat<-rbeta(1, 1+npics[1], 1+npics[2]-npics[1]) #npics[1]=identified, npics[2]=all marked
	cr[marked]<-crat
 	}

        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern")) 
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma * sigma))
#            lam.cand[lam.cand==0] <- 1e-300 
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

        for (r in 1:R) {
            zip <- lam[, r] * z
            for (t in 1:T) {
                if (n[r, t] == 0) {
                  Y[!marked, r, t] <- 0
                  next
                }
                unmarked <- !Ydata[, r, t]
                nUnknown <- n[r, t] 
                if (nUnknown == 0) 
                  next
                probs <- zip[unmarked]
                probs <- probs/sum(probs)
                if (identical(obsmod, "pois")) 
                  Y[unmarked, r, t] <- rmultinom(1, nUnknown, probs)
                else if (identical(obsmod, "bern")) {
                  Y[unmarked, r, t] <- 0
                  guy <- sample(which(unmarked), nUnknown, prob = probs)
                  Y[guy, r, t] <- 1
                }
            }
        }
               
        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + (M-sum(marked)) - sum(z[!marked])

        Sups <-Skups<- 0
        for (i in 1:M) {

	if (tel[i]) {
            Scand <- c(rnorm(1, S[i, 1], delta[3]), rnorm(1, S[i, 2], delta[3]))
		    } else {
            Scand <- c(rnorm(1, S[i, 1], delta[4]), rnorm(1, S[i, 2], delta[4]))
                           }

            inbox <- Scand[1] >= xlims[1] & Scand[1] <= xlims[2] & 
                Scand[2] >= ylims[1] & Scand[2] <= ylims[2]

            if (!inbox) 
                next
            dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2*sigma*sigma))


############# telemetry model ####################################################
	if (tel[i]) {
# still need to match locs with i vector
	ll<-sum(dmvnorm(x=locs[[sum(tel[1:i])]],mean=c(S[i,1],S[i,2]), sigma=cbind(c(sigma^2,0), c(0,sigma^2)), log=T))   
	llcand<-sum(dmvnorm(x=locs[[sum(tel[1:i])]],mean=c(Scand[1],Scand[2]), sigma=cbind(c(sigma^2,0), c(0,sigma^2)), log=T))   

            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i, ] <- Scand
                lam[i, ] <- lam.cand
                D[i, ] <- dtmp
                Skups <- Skups + 1
            }

} else { #end of if tel statement, begin non tel update

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
          } # end of no tel part
        }
        if (iter%%100 == 0) {
            cat("   Acceptance rates\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/(M-length(locs)), "\n")
            cat("     Sk =", Skups/length(locs), "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, sum(z))
    }
    return(out)
}
