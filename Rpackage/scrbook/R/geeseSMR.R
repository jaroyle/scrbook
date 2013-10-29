

# MCMC with Z partially observed

geeseSMR <- function (n, X, y, M, mmax, obsmod = c("pois", "bern"),niters, npics, Sex,EffAr,
    inits, delta, SSp ) 
{
    obsmod <- match.arg(obsmod)


Eff<-apply(EffAr, 1:2, sum)

    SEX<-c(Sex, rep(NA, (M+mmax)-nrow(y)))
 knownS<-!is.na(SEX)

    J <- nrow(n)
    K <- ncol(n)
    S <- cbind(inits$S$x, inits$S$y)

iv<-(1:nrow(y))[apply(y,1,sum)>0]
lc<-list()
for (i in 1:length(iv)){
lc[[i]]<-(1:J)[apply(y[iv[i],,],1,sum)>0]
SxK<-mean(X[lc[[i]],1])
SyK<-mean(X[lc[[i]],2])
S[iv[i],]<-c(SxK, SyK)
}

    D <- e2dist(S, X)
    sigma <- inits$sigma
    lam0 <- inits$lam0
    lam<-array(NA, c(M+mmax, J,2))

    lam[,,1] <- lam0 * exp(-(D * D)/(2 * sigma[1] * sigma[1]))
    lam[,,2] <- lam0 * exp(-(D * D)/(2 * sigma[2] * sigma[2]))

    nObs <- nrow(y)
    Y<-array(0,c(M+mmax, J,K))
        Y[1:nObs, , ] <- y  #augmented SCR matrix
    marked<-rep(FALSE, M+mmax)
    marked[1:mmax]<-TRUE
    psi <- inits$psi
    psim<-inits$psim
    z <- rbinom(M+mmax, 1, psi)
    z[1:nObs]<-1
    phi<-inits$phi
    sex<-rbinom(M+mmax, 1, phi)+1
    sex[knownS]<-SEX[knownS]

    probs<-matrix(NA, nrow=M,ncol=J)
    for (i in 1:M){
    probs[i,]<-lam[(mmax+i),,sex[(mmax+i)]]*z[(mmax+i)]
    }

    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            nUnknown <- n[j, k] 
	    pt<-probs[,j]*EffAr[!marked,j,k] 
            pt <- pt/sum(pt)
            if (identical(obsmod, "pois")) 
                Y[!marked, j, k] <- rmultinom(1, nUnknown, pt)
            else if (identical(obsmod, "bern")) {
                Y[!marked, j, k] <- 0
                guys <- sample((mmax+1):(M+mmax), nUnknown, prob = pt)
                Y[guys, j, k] <- 1
            }
        }
    }

Ymat<-apply(Y,1:2,sum)

	cr<-rep(1,M+mmax) 
	if(missing(npics)){
	crat<-1} else{
	crat<-npics[1]/npics[2]}
        cr[marked]<-crat

    out <- matrix(NA, nrow = niters, ncol = 8)
    colnames(out) <- c("sigmaF","sigmaM", "lam0","c", "psi", "phi","m", "N")
    cat("\nstarting values =", c(sigma, lam0, crat, psi, phi,sum(z[marked]), sum(z)), "\n\n")
    for (iter in 1:niters) {
        if (iter%%100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), 
                "\n")
            cat("   current =", out[iter - 1, ], "\n")
        }
	ll<-NULL

	for (sx in 1:2){
          if(identical(obsmod, "pois")) {
            ll[sx] <- sum(dpois(Ymat[sex==sx,], lam[sex==sx,,sx]*cr[sex==sx] *Eff[sex==sx,]* z[sex==sx], log = TRUE))
          } else if(identical(obsmod, "bern")) {
           ll[sx] <- sum(dbinom(Ymat[sex==sx,], Eff[sex==sx,], lam[sex==sx,,sx]*cr[sex==sx] * z[sex==sx], log = TRUE))
          }
	}

	if(!missing(npics)){
	crat<-rbeta(1, 1+npics[1], 1+npics[2]-npics[1]) #npics[1]=identified, npics[2]=all marked
	cr[marked]<-crat
 	}


	for (sx in 1:2){

        sigma.cand <- rnorm(1, sigma[sx], delta[1])
        if (sigma.cand > 0) {
 
	 lam.cand <- lam0 * exp(-(D * D)/(2 * sigma.cand * sigma.cand))
            if (identical(obsmod, "pois")) {
                ll[sx] <- sum(dpois(Ymat[sex==sx,], lam[sex==sx,,sx]*cr[sex==sx] *Eff[sex==sx,]* z[sex==sx], log = TRUE))
                llcand <- sum(dpois(Ymat[sex==sx,], lam.cand[sex==sx,]*cr[sex==sx]*Eff[sex==sx,] * z[sex==sx], log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll[sx] <- sum(dbinom(Ymat[sex==sx,], Eff[sex==sx,], lam[sex==sx,,sx]*cr[sex==sx] * z[sex==sx], log = TRUE)) 
                llcand <- sum(dbinom(Ymat[sex==sx,], Eff[sex==sx,], lam.cand[sex==sx,]*cr[sex==sx] * z[sex==sx], log = TRUE)) 
            }
            if (runif(1) < exp(llcand  - ll[sx] )) {
                ll[sx] <- llcand
                lam[,,sx] <- lam.cand
                sigma[sx] <- sigma.cand
            }
        }
	}

        lam0.cand <- rnorm(1, lam0, delta[2])
        test2 <- TRUE
        if (identical(obsmod, "bern")) 
            test2 <- lam0.cand <= 1
        if (lam0.cand >= 0 & test2) {
	    lam.cand<-array(NA, c(M+mmax, J,2))
            lam.cand[,,1] <- lam0.cand * exp(-(D * D)/(2 * sigma[1] * sigma[1]))
            lam.cand[,,2] <- lam0.cand * exp(-(D * D)/(2 * sigma[2] * sigma[2]))
            if (identical(obsmod, "pois")) 
	     llcand <- (sum(dpois(Ymat[sex==1,], Eff[sex==1,]*lam.cand[sex==1,,1]*z[sex==1]*cr[sex==1], log=TRUE)) 
			+ sum(dpois(Ymat[sex==2,], Eff[sex==2,]*lam.cand[sex==2,,2]*z[sex==2]*cr[sex==2], log=TRUE))   )

            else if (identical(obsmod, "bern")) 

	     llcand <- (sum(dbinom(Ymat[sex==1,], Eff[sex==1,],lam.cand[sex==1,,1]*z[sex==1]*cr[sex==1], log=TRUE)) 
			+ sum(dbinom(Ymat[sex==2,], Eff[sex==2,],lam.cand[sex==2,,2]*z[sex==2]*cr[sex==2], log=TRUE))   )

            if (runif(1) < exp((llcand) - (sum(ll)))) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

##update z, different priors for marked and unmarked
	probs<-matrix(nrow=M+mmax, ncol=J)
   	for (i in 1:(M+mmax)){
	probs[i,]<-lam[i,,sex[i]]
	}

zcand<-ifelse(z == 0, 1, 0)
            if (identical(obsmod, "pois")) {
                llz <- rowSums(dpois(Ymat, probs * cr* z * Eff, log = TRUE))
                llcandz <- rowSums(dpois(Ymat, probs * cr* zcand * Eff, log = TRUE))
            } else if (identical(obsmod, "bern")) {
                llz <- rowSums(dbinom(Ymat,Eff, probs * cr* z , log = TRUE))
                llcandz <- rowSums(dbinom(Ymat,Eff, probs * cr* zcand , log = TRUE))
            }

            prior <- dbinom(z, 1, psim, log = TRUE)
	    prior[(mmax+1):(M+mmax)]<-dbinom(z[(mmax+1):(M+mmax)], 1, psi, log = TRUE)
            prior.cand <- dbinom(zcand, 1, psim, log = TRUE)
            prior.cand[(mmax+1):(M+mmax)] <- dbinom(zcand[(mmax+1):(M+mmax)], 1, psi, log = TRUE)

		rat<-(llcandz + prior.cand) - (llz + prior)

            kp<-runif(M+mmax,0,1) < exp(rat) 
                z[kp] <- zcand[kp]
                zUps <- sum(kp)
            


 probs<-matrix(NA, nrow=M,ncol=J)
    for (i in 1:M){
    probs[i,]<-lam[(mmax+i),,sex[(mmax+i)]]*z[(mmax+i)]
    }

    for (j in 1:J) {
        for (k in 1:K) {
            if (n[j, k] == 0) {
                Y[!marked, j, k] <- 0
                next
            }
            nUnknown <- n[j, k] 
	    pt<-probs[,j]*EffAr[!marked,j,k] 
            pt <- pt/sum(pt)
            if (identical(obsmod, "pois")) 
                Y[!marked, j, k] <- rmultinom(1, nUnknown, pt)
            else if (identical(obsmod, "bern")) {
                Y[!marked, j, k] <- 0
                guys <- sample((mmax+1):(M+mmax), nUnknown, prob = pt)
                Y[guys, j, k] <- 1
            }
        }
    }

       Ymat<-apply(Y,1:2,sum)        

        psim <- rbeta(1, 1 + sum(z[marked]), 1 + mmax - sum(z[marked]))

        psi <- rbeta(1, 1 + sum(z[!marked]), 1 + M - sum(z[!marked]))

#####update sex
            sex.cand <- ifelse(sex==1, 2, 1)

	probs<-probs.cand<-matrix(nrow=M+mmax, ncol=J)
   	for (i in 1:(M+mmax)){
	probs[i,]<-lam[i,,sex[i]]
	probs.cand[i,]<-lam[i,,sex.cand[i]]
	}

if (identical(obsmod, "pois")) {
                lls <- rowSums(dpois(Ymat[!knownS,], probs[!knownS,] * cr[!knownS]* z[!knownS] * Eff[!knownS,], log = TRUE))
                llcand <- rowSums(dpois(Ymat[!knownS,], probs.cand[!knownS,] * cr[!knownS]* z[!knownS] * Eff[!knownS,], log = TRUE))
            }else if (identical(obsmod, "bern")) {
                lls <- rowSums(dbinom(Ymat[!knownS,], Eff[!knownS,], probs[!knownS,] * cr[!knownS]* z[!knownS] , log = TRUE))
                llcand <- rowSums(dbinom(Ymat[!knownS,], Eff[!knownS,], probs.cand[!knownS,] * cr[!knownS]* z[!knownS] , log = TRUE))
            }
            prior <- dbinom(sex[!knownS]-1, 1, phi, log=TRUE)
            prior.cand <- dbinom(sex.cand[!knownS]-1, 1, phi, log=TRUE)
		kks<-runif(sum(!knownS),0,1) < exp( (llcand+prior.cand) - (lls+prior) )
                sex[!knownS][kks] <- sex.cand[!knownS][kks]
                sUps <- sum(kks)

	#update phi
	allmale=sum(sex-1)
        phi<-rbeta(1,1+allmale,1+(M+mmax)-allmale)  #


### update activity centers
     #    Sups <- 0

            Scand <- cbind(rnorm(M+mmax, S[,1], delta[3]),
                       rnorm(M+mmax, S[,2], delta[3]))

	Scoord<-SpatialPoints(Scand)
	SinPoly<-over(Scoord,SSp)

            dtmp <- e2dist(Scand, X)
	IN<-!is.na(SinPoly)

	  lam.cand=array(NA, c(M+mmax,J,2))
                lam.cand[,,1] <- lam0*exp(-(dtmp*dtmp)/(2*sigma[1]*sigma[1]) )
                lam.cand[,,2] <- lam0*exp(-(dtmp*dtmp)/(2*sigma[2]*sigma[2]) ) 

	probs<-probs.cand<-matrix(nrow=M+mmax, ncol=J)
   	for (i in 1:(M+mmax)){
	probs[i,]<-lam[i,,sex[i]]
	probs.cand[i,]<-lam.cand[i,,sex[i]]
	}

            if (identical(obsmod, "pois")) {
                ll <- rowSums(dpois(Ymat[IN, ], probs[IN,] *cr[IN] * z[IN] * Eff[IN,], log = TRUE))
                llcand <- rowSums(dpois(Ymat[IN, ], probs.cand[IN,] *cr[IN] * z[IN] * Eff[IN,], log = TRUE))
            }
            else if (identical(obsmod, "bern")) {
                ll <- rowSums(dbinom(Ymat[IN, ],Eff[IN,], probs[IN,] *cr[IN] * z[IN], log = TRUE))
                llcand <- rowSums(dbinom(Ymat[IN, ],Eff[IN,], probs.cand[IN,] *cr[IN] * z[IN], log = TRUE))
            }
		kkS<-runif(sum(IN),0,1) < exp(llcand - ll)
		S[IN,1][kkS]<-Scand[IN,1][kkS]
		S[IN,2][kkS]<-Scand[IN,2][kkS]
	      Sups<-sum(kkS)


        if (iter%%100 == 0) {
            cat("   Acceptance rates\n")
            cat("     z =", zUps/(M+mmax), "\n")
            cat("     S =", Sups/(M+mmax), "\n")
            cat("     sex =", sUps/((M+mmax)-sum(knownS)), "\n")
        }
        out[iter, ] <- c(sigma, lam0, crat, psi, phi,sum(z[marked]), sum(z))
    }
    return(out)
}

