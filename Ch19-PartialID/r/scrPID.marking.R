multidens<-function(x,p){
n=nrow(x)
llp<-NULL
for (i in 1:n){
llp[i]<-dmultinom(x[i,],1,p[i,], log=TRUE)
}
ll<-sum(llp)
return(ll)

}


scrPIDm<-function (n, X, y, ymark, capsites, M, obsmod = c("pois", "bern"),nmarked=c("known", "unknown"), niters, npics,
    xlims, ylims, inits, delta ) 
{
    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)

    if(identical(nmarked, "unknown") & !missing(npics)) 
	stop ("Need to know number of marked individuals if individual identification of marks is imperfect")

    J <- nrow(n)
    K <- ncol(n)
    S <- inits$S

for(i in 1:dim(ymark)[1]){
    S[i,]<-X[capsites[which(ymark[i,]==1)],]
}
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
                Y[, j, k] <- 0
                next
            }
            unmarked <- !Ydata[, j, k]
            nUnknown <- n[j, k] - sum(Y[!unmarked, j,k])
            if (nUnknown < 0) 
                browser()
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

### add in piece about marking event ###
Icap<-matrix(0, nrow=M, ncol=length(capsites)+1)
Icap[,length(capsites)+1]<-1
Icap[1:nMarked,]<-ymark
#Icap<-ymark

p0<-inits$p0
lam.cap<-exp(p0 +(- (D[,capsites]*D[,capsites]))/(2*sigma*sigma)) 
ptot<-rowSums(lam.cap)
lam.cap<-(lam.cap / (1 + ptot)) * z
lam.cap<-cbind(lam.cap, 1-rowSums(lam.cap[, 1:length(capsites)]))

#lam.cap<-exp(- (D[marked,capsites]*D[marked,capsites])/(2*sigma*sigma)) 
#ptot<-rowSums(lam.cap)
#lam.cap<-lam.cap/ptot

    out <- matrix(NA, nrow = niters, ncol = 6)
    colnames(out) <- c("sigma", "lam0","p0","c", "psi", "N") #"p0",
    cat("\nstarting values =", c(sigma, lam0, p0,crat, psi, sum(z)), "\n\n")
    
for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), 
                "\n")
            cat("   current =", out[iter - 1, ], "\n")
        }


     #   if(identical(obsmod, "pois")) {
     #       ll <- sum(dpois(Y, lam*cr * z, log = TRUE)) + multidens(Icap, lam.cap)
     #  } else if(identical(obsmod, "bern")) {
     #      ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) + multidens(Icap, lam.cap)
     #  }

	if(!missing(npics)){
	crat<-rbeta(1, 1+npics[1], 1+npics[2]-npics[1]) #npics[1]=identified, npics[2]=all marked
	cr[marked]<-crat
 	}

        sigma.cand <- rnorm(1, sigma, delta[1])
        if (sigma.cand > 0) {

            lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand * sigma.cand))

	lam.cap.cand<-exp(p0 + (- (D[,capsites]*D[,capsites]))/(2*sigma.cand*sigma.cand)) 
	ptot<-rowSums(lam.cap.cand)
	lam.cap.cand<-(lam.cap.cand / (1 + ptot)) * z
	lam.cap.cand<-cbind(lam.cap.cand, 1-rowSums(lam.cap.cand[, 1:length(capsites)]))

            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam*cr * z, log = TRUE)) + multidens(Icap, lam.cap)
                llcand <- sum(dpois(Y, lam.cand*cr * z, log = TRUE)) + multidens(Icap, lam.cap.cand)
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) + multidens(Icap, lam.cap)
                llcand <- sum(dbinom(Y, 1, lam.cand*cr * z, log = TRUE)) + multidens(Icap, lam.cap.cand)
            }
            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                lam <- lam.cand
                lam.cap <- lam.cap.cand
                sigma <- sigma.cand
            }
        }

####lam0 only conditional on resighting data
        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern")) 
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
            lam.cand <- lam0.cand * exp(-(D * D)/(2 * sigma * sigma))
            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y, lam*cr * z, log = TRUE)) 
                llcand <- sum(dpois(Y, lam.cand*cr * z, log = TRUE)) }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y, 1, lam*cr * z, log = TRUE)) 
                llcand <- sum(dbinom(Y, 1, lam.cand*cr * z, log = TRUE)) }
            if (runif(1) < exp((llcand) - (ll))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

###p0 only conditional on initial marking data
        p0.cand <- rnorm(1, p0, delta[3])

            lam.cap.cand <-exp(p0.cand +(- (D[,capsites]*D[,capsites]))/(2*sigma*sigma)) 
	ptot<-rowSums(lam.cap.cand)
	lam.cap.cand<-(lam.cap.cand / (1 + ptot)) * z
	lam.cap.cand<-cbind(lam.cap.cand, 1-rowSums(lam.cap.cand[, 1:length(capsites)]))

	ll<-multidens(Icap, lam.cap)
	llcand<-multidens(Icap, lam.cap.cand)

            if (runif(1) < exp((llcand) - (ll))) {
                p0 <- p0.cand
                lam.cap <- lam.cap.cand
            }
        


        zUps <- 0
        seen <- apply(Y > 0, 1, any)

        zcand <- ifelse(z == 0, 1, 0)
        lam.cap.cand <-exp(p0 +(- (D[,capsites]*D[,capsites]))/(2*sigma*sigma)) 
	ptot<-rowSums(lam.cap.cand)
	lam.cap.cand<-(lam.cap.cand / (1 + ptot)) * zcand
	lam.cap.cand<-cbind(lam.cap.cand, 1-rowSums(lam.cap.cand[, 1:length(capsites)]))

        for (i in 1:M) {
            if (seen[i] | marked[i]) 
                next

            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] * z[i], log = TRUE)) #+ dmultinom(Icap[i,],1, lam.cap[i,], log=TRUE)
                llcand <- sum(dpois(Y[i, , ], lam[i, ] * zcand[i], log = TRUE)) #+ dmultinom(Icap[i,],1, lam.cap.cand[i,], log=TRUE)
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] * z[i], log = TRUE)) #+ dmultinom(Icap[i,],1, lam.cap[i,], log=TRUE)
                llcand <- sum(dbinom(Y[i, , ], 1, lam[i, ] * zcand[i], log = TRUE)) #+ dmultinom(Icap[i,],1, lam.cap.cand[i,], log=TRUE)
            }
            prior <- dbinom(z[i], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand[i], 1, psi, log = TRUE)
            if (runif(1) < exp((llcand + prior.cand) - (ll + prior))) {
                z[i] <- zcand[i]
		#lam.cap[i,]<-lam.cap.cand[i,]
                zUps <- zUps + 1
            }
        }


        for (j in 1:J) {
            zip <- lam[, j] * z
            for (k in 1:K) {
                if (n[j, k] == 0) {
                  Y[, j, k] <- 0
                  next
                }
                unmarked <- !Ydata[, j, k]
                nUnknown <- n[j, k] - sum(Y[!unmarked, j, k])
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
               
            if (identical(nmarked, "known")) {
        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + (M-sum(marked)) - sum(z[!marked]))
	} else if (identical(nmarked, "unknown")) {
        psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))
	}

        Sups <- 0

            Scand <- cbind(rnorm(M, S[, 1], delta[4]), rnorm(M, S[, 2], delta[4]))
            inbox <- Scand[,1] >= xlims[1] & Scand[,1] <= xlims[2] & 
                Scand[,2] >= ylims[1] & Scand[,2] <= ylims[2]

            dtmp <- e2dist(Scand, X)
            lam.cand <- lam0 * exp(-(dtmp * dtmp)/(2*sigma*sigma))

            lam.cap.cand <-exp(p0 + (- (dtmp[,capsites]*dtmp[,capsites]))/(2*sigma*sigma)) 
	ptot<-rowSums(lam.cap.cand)
	lam.cap.cand<-(lam.cap.cand / (1 + ptot)) * z
	lam.cap.cand<-cbind(lam.cap.cand, 1-rowSums(lam.cap.cand[, 1:length(capsites)]))


        for (i in 1:M) {

            if (!inbox[i]) 
                next


	   # if(marked[i]){  #marked

            if (identical(obsmod, "pois")) {
                ll <- sum(dpois(Y[i, , ], lam[i, ] *cr[i] * z[i], log = TRUE)) + dmultinom(Icap[i,],1, lam.cap[i,], log=TRUE)
                llcand <- sum(dpois(Y[i, , ], lam.cand[i,] *cr[i]* z[i], log = TRUE)) + dmultinom(Icap[i,],1, lam.cap.cand[i,], log=TRUE)
            }
            else if (identical(obsmod, "bern")) {
                ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] *cr[i]* z[i],log = TRUE)) + dmultinom(Icap[i,],1, lam.cap[i,], log=TRUE)
                llcand <- sum(dbinom(Y[i, , ], 1, lam.cand[i,] *cr[i]* z[i], log = TRUE)) + dmultinom(Icap[i,],1, lam.cap.cand[i,], log=TRUE)
            }

	   # } else { #unmarked

           # if (identical(obsmod, "pois")) {
           #     ll <- sum(dpois(Y[i, , ], lam[i, ] *cr[i] * z[i], log = TRUE))
           #     llcand <- sum(dpois(Y[i, , ], lam.cand[i,] *cr[i]* z[i], log = TRUE)) 
           # }
           # else if (identical(obsmod, "bern")) {
           #     ll <- sum(dbinom(Y[i, , ], 1, lam[i, ] *cr[i]* z[i],log = TRUE)) 
           #     llcand <- sum(dbinom(Y[i, , ], 1, lam.cand[i,] *cr[i]* z[i], log = TRUE)) 
           # }

	   # }

            if (runif(1) < exp(llcand - ll)) {
                ll <- llcand
                S[i, ] <- Scand[i, ]
                lam[i, ] <- lam.cand[i, ]

	    #if(marked[i]){ 
                lam.cap[i, ] <- lam.cap.cand[i, ]
	    #	}

                D[i, ] <- dtmp[i, ]
                Sups <- Sups + 1
            }
        }

        if (iter%%100 == 0) {
            cat("   Acceptance rates\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }
        out[iter, ] <- c(sigma, lam0,p0,  crat, psi, sum(z)) #p0,
    }
    return(out)
}
