

# MCMC with Z partially observed

geeseSMR <- function(y, X, Zknown, M,EffAr, niters, inits,Sex, tune, SSp) {

Eff<-apply(EffAr, 1:2, sum)

    SEX<-c(Sex, rep(NA, M-nrow(Zknown)))
 knownS<-!is.na(SEX)
    R <- nrow(y)
    T <- ncol(y)
    S <- cbind(inits$S$x, inits$S$y)

iv<-(1:nrow(Zknown))[apply(Zknown,1,sum)>0]
lc<-list()
for (i in 1:length(iv)){
lc[[i]]<-(1:R)[apply(Zknown[iv[i],,],1,sum)>0]
SxK<-mean(X[lc[[i]],1])
SyK<-mean(X[lc[[i]],2])
S[iv[i],]<-c(SxK, SyK)
}

    D <- e2dist(S, X)
    theta <- inits$theta
    lam0 <- inits$lam0
lam<-array(NA, dim=c(M,R, 2))
    lam[,,1] <- lam0*exp(-(D*D)/(2*theta[1]*theta[1]))
    lam[,,2] <- lam0*exp(-(D*D)/(2*theta[2]*theta[2]))
    lam[lam<1e-30] <- 1e-30
    w <- inits$w
    psi <- inits$psi
    phi<-inits$phi
    sex<-rbinom(M,1,phi)+1

Z <- array(NA, c(M,R,T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if( !is.null(Zknown)) {   #!missing(Zknown) ||
        nMarked <- nrow(Zknown)
        marked[1:nMarked] <- TRUE
        Z[1:nMarked,,] <- Zknown
    }
    w[marked] <- 1
    sex[knownS]<-SEX[knownS]
    Zdata <- !is.na(Z)

probs<-matrix(NA, nrow=M,ncol=R)
for (i in 1:M){
probs[i,]<-lam[i,,sex[i]]*w[i]
}

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

	pt<-probs[,r]*EffAr[,r,t]  #pt is an M x R array 

	guyIDs <- sample((nMarked+1):M, nUnknown, replace=FALSE, prob=pt[!marked]) 		
	noGuys<-(1:M)[-c(1:nMarked,guyIDs)]
            Z[guyIDs,r,t] <- 1
		Z[noGuys,r,t] <- 0
        }
    }

Zmat<-apply(Z, 1:2, sum)

    out <- matrix(NA,nrow=niters,ncol=6)
    colnames(out) <- c("sigmaF", "sigmaM","lam0", "psi","phi", "N")

    cat("\nstarting values =", c(theta, lam0, psi, phi, sum(w)), "\n\n")


    for(iter in 1:niters) {

        if(iter %% 10 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }


ll<-c(NA, NA)
for (se in 1:2){
ll[se]<-sum(dbinom(Zmat[sex==se,],Eff, lam[sex==se,,se]*w[sex==se], log=TRUE)) 
}

        # update theta

	    #llcand<-c(NA, NA)
	for (se in 1:2){
        theta.cand <- rnorm(1, theta[se], tune$sig[se])
       
        if(theta.cand>0){

            lam.cand <- lam0*exp(-(D*D)/(2*theta.cand*theta.cand))

	    ll[se]<-sum(dbinom(Zmat[sex==se,],Eff[sex==se,], lam[sex==se,,se]*w[sex==se], log=TRUE)) 
            llcand <- sum(dbinom(Zmat[sex==se,],Eff[sex==se,], lam.cand[sex==se,]*w[sex==se], log=TRUE))

            if(runif(1) < exp( llcand - ll[se] ) ){
                ll[se]<-llcand
                lam[,,se]<-lam.cand
                theta[se]<-theta.cand
            }
        }
	} #end sex loop

        # update lam0
        lam0.cand <- rnorm(1, lam0, tune$lam0)
        if(lam0.cand>0 & lam0.cand<1) {
            
            lam.cand1 <- lam0.cand*exp(-(D*D)/(2*theta[1]*theta[1]))
            lam.cand2 <- lam0.cand*exp(-(D*D)/(2*theta[2]*theta[2]))

            llcand <- sum(dbinom(Zmat[sex==1,], Eff[sex==1,],lam.cand1[sex==1,]*w[sex==1], log=TRUE)) + sum(dbinom(Zmat[sex==2,], Eff[sex==2,],lam.cand2[sex==2,]*w[sex==2], log=TRUE)) 

            if(runif(1) < exp(llcand - sum(ll))) {
                lam0 <- lam0.cand
                lam[,,1] <- lam.cand1
                lam[,,2] <- lam.cand2
            }
        }



      ### update "w" here
   	for (i in 1:M){
	probs[i,]<-lam[i,,sex[i]]
	}

       seen <- apply(Z>0, 1, any)
 	    
           wcand <- ifelse(w==0, 1, 0)

         ll <- rowSums(dbinom(Zmat,Eff, probs*w, log=TRUE))
           llcand <- rowSums(dbinom(Zmat,Eff, probs*wcand, log=TRUE))

           prior <- dbinom(w, 1, psi, log=TRUE)
          prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
           kp<-runif(M, 0, 1) < exp( (llcand+prior.cand) - (ll+prior) ) & !marked
               w[kp] <- wcand[kp]
               wUps <- sum(kp)


       ### update "sex" here
        sUps <- 0
        for(i in 1:M) {
            if(knownS[i]) #no need to update marked guys with known sex
                next

	    if (w[i]==1) {

            sex.cand <- ifelse(sex[i]==1, 2, 1)

            lls <- sum(dbinom(Zmat[i,],Eff[i,], lam[i,,sex[i]], log=TRUE))
            llcand <- sum(dbinom(Zmat[i,],Eff[i,], lam[i,,sex.cand], log=TRUE))

            prior <- dbinom(sex[i]-1, 1, phi, log=TRUE)
            prior.cand <- dbinom(sex.cand-1, 1, phi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (lls+prior) )) {
                sex[i] <- sex.cand
                sUps <- sUps+1
            } 								# end second if statement
        } else { 							#end first if statement; add what happens if W1[i]=0
	sex[i]<-rbinom(1,1,phi)+1
		}

	} #end M loop

	#update phi
	allmale=sum(sex-1)
        phi<-rbeta(1,1+allmale,1+M-allmale)  #


        # update Z

	Z[(nMarked+1):M,,]<-0  #give'em all 0's 
	for (i in 1:M){
	probs[i,]<-lam[i,,sex[i]]*w[i]
	}

        for(r in 1:R) {
            for(t in 1:T) {
		pt<-probs[,r]*EffAr[,r,t] 

                unmarked <- !Zdata[,r,t]
                nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown < 0)
                    browser()
                if(nUnknown == 0)
   		next

	guyIDs <- sample((nMarked+1):M, nUnknown, replace=FALSE, prob=pt[!marked]) 		
	noGuys<-(1:M)[-c(1:nMarked,guyIDs)]
            Z[guyIDs,r,t] <- 1
		Z[noGuys,r,t] <- 0

            }
        }

	Zmat<-apply(Z, 1:2, sum)

        # update psi
        psi<-rbeta(1,1+sum(w),1+M-sum(w))

        # update S
        Sups <- 0

            Scand <- cbind(rnorm(M, S[,1], tune$S),
                       rnorm(M, S[,2], tune$S))

	Scoord<-SpatialPoints(Scand)
	SinPoly<-over(Scoord,SSp)


        for(i in 1:M) {   # note this is "M" in general

            if(!is.na(SinPoly[i])) {
                dtmp <-  sqrt((Scand[i,1] - X[,1])^2 + (Scand[i,2] - X[,2])^2)

	  lam.cand=matrix(NA, nrow=R, ncol=2)
                lam.cand[,1] <- lam0*exp(-(dtmp*dtmp)/(2*theta[1]*theta[1]) )
                lam.cand[,2] <- lam0*exp(-(dtmp*dtmp)/(2*theta[2]*theta[2]) ) 


                ll <- sum(dbinom(Zmat[i,],Eff[i,], lam[i,,sex[i]]*w[i], log=TRUE))
                llcand <- sum(dbinom(Zmat[i,],Eff[i,], lam.cand[,sex[i]]*w[i], log=TRUE))

                if(runif(1) < exp(llcand - ll)) {
                    S[i,] <- Scand[i,]
                    lam[i,,1] <- lam.cand[,1]
                    lam[i,,2] <- lam.cand[,2]
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
        }

        if(iter %% 10 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/sum(!marked), "\n")
            cat("     sex =", sUps/sum(!knownS), "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(theta,lam0,psi,phi,sum(w) )

    }

    return(out)
}


















