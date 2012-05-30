

# MCMC with Z partially observed

pIDgeese <- function(y, X, Zknown, M,Eff, niter, xl,xu,yl,yu, inits,Sex, delta) {
M<-M
    SEX<-c(Sex, rep(NA, M-nrow(Zknown)))
 knownS<-!is.na(SEX)
    R <- nrow(y)
    T <- ncol(y)
    S <- inits$S

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

	    probs<-rep(NA, M)

		for (m in 1:M){
            probs[m] <- lam[m,r,sex[m]]*w[m]*Eff[m,r,t] 
			}

	guyIDs <- sample((nMarked+1):M, nUnknown, replace=FALSE, prob=probs[!marked]) 		
	noGuys<-(1:M)[-c(1:nMarked,guyIDs)]
            Z[guyIDs,r,t] <- 1
		Z[noGuys,r,t] <- 0
        }
    }

    out <- matrix(NA,nrow=niters,ncol=6)
    colnames(out) <- c("sigmaF", "sigmaM","lam0", "psi","phi", "N")

    cat("\nstarting values =", c(theta, lam0, psi, sum(w)), "\n\n")


    for(iter in 1:niters) {

        if(iter %% 10 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

ll<-c(NA, NA)
for (se in 1:2){
lam.arr<- array(lam[,,se],c(M,R,T))  * (array(w,c(M,R,T)) )* Eff
ll[se]<-sum(dbinom(Z[sex==se,,],1, lam.arr[sex==se,,], log=TRUE)) 
}
        # update theta

	for (se in 1:2){
        theta.cand <- rnorm(1, theta[se], tune$sig[se])
       
        if(theta.cand>0){

            lam.cand <- lam0*exp(-(D*D)/(2*theta.cand*theta.cand))

    	    lam.arr<-array(lam[,,se], c(M,R,T)) * array(w, c(M,R,T)) * Eff
  	    lam.cand.arr<-array(lam.cand, c(M,R,T)) * array(w, c(M,R,T)) * Eff

            ll[se] <- sum(dbinom(Z[sex==se,,], 1,lam.arr[sex==se,,], log=TRUE))
            llcand <- sum(dbinom(Z[sex==se,,],1, lam.cand.arr[sex==se,,], log=TRUE))

            if(runif(1) < exp( llcand  - ll[se] ) ){
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
  	    lam.cand.arr1<-array(lam.cand1, c(M,R,T)) * array(w, c(M,R,T)) * Eff
            lam.cand2 <- lam0.cand*exp(-(D*D)/(2*theta[2]*theta[2]))
  	    lam.cand.arr2<-array(lam.cand2, c(M,R,T)) * array(w, c(M,R,T)) * Eff

            llcand <- sum(dbinom(Z[sex==1,,], 1,lam.cand.arr1[sex==1,,], log=TRUE)) + sum(dbinom(Z[sex==2,,], 1,lam.cand.arr2[sex==2,,], log=TRUE)) 

            if(runif(1) < exp(llcand - sum(ll))) {
                lam0 <- lam0.cand
                lam[,,1] <- lam.cand1
                lam[,,2] <- lam.cand2
            }
        }


        ### update "w" here
        wUps <- 0
       seen <- apply(Z>0, 1, any)
      for(i in 1:M) {
           if(seen[i] | marked[i])
               next
           wcand <- ifelse(w[i]==0, 1, 0)

	lam.mat<-matrix(lam[i,,sex[i]], nrow=R, ncol=T) * matrix(w[i], nrow=R, ncol=T) *Eff[i,,]
	lam.cand.mat<-matrix(lam[i,,sex[i]], nrow=R, ncol=T) * matrix(wcand, nrow=R, ncol=T) *Eff[i,,]

         ll <- sum(dbinom(Z[i,,],1, lam.mat, log=TRUE))
           llcand <- sum(dbinom(Z[i,,],1, lam.cand.mat, log=TRUE))

           prior <- dbinom(w[i], 1, psi, log=TRUE)
          prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
           if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
               w[i] <- wcand
               wUps <- wUps+1
           }
       }


       ### update "sex" here
        sUps <- 0
        for(i in 1:M) {
            if(knownS[i]) #no need to update marked guys with known sex
                next

	    if (w[i]==1) {

            sex.cand <- ifelse(sex[i]==1, 2, 1)

lam.mat<- matrix(lam[i,,sex[i]], nrow=R, ncol=T) * Eff[i,,]
lam.mat.cand<- matrix(lam[i,,sex.cand], nrow=R, ncol=T) * Eff[i,,]

            lls <- sum(dpois(Z[i,,],lam.mat, log=TRUE))
            llcand <- sum(dpois(Z[i,,], lam.mat.cand, log=TRUE))

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

        for(r in 1:R) {
            for(t in 1:T) {
		for (m in 1:M){
            probs[m] <- lam[m,r,sex[m]]*w[m]*Eff[m,r,t] 
			}
            probs <- probs[!marked]

                unmarked <- !Zdata[,r,t]
                nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown < 0)
                    browser()
                if(nUnknown == 0)
   		next

	guyIDs <- sample((nMarked+1):M, nUnknown, replace=FALSE, prob=probs) 		
	noGuys<-(1:M)[-c(1:nMarked,guyIDs)]
            Z[guyIDs,r,t] <- 1
		Z[noGuys,r,t] <- 0

            }
        }


        # update psi
        psi<-rbeta(1,1+sum(w[!marked]),1+sum(!marked)-sum(w[!marked]))

        # update S
        Sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            Scand <- c(rnorm(1, S[i,1], tune$S),
                       rnorm(1, S[i,2], tune$S))
            inbox <- Scand[1]>=xl & Scand[1]<=xu &
                     Scand[2]>=yl & Scand[2]<=yu
            if(inbox) {
                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)

	  lam.cand=matrix(NA, nrow=R, ncol=2)
                lam.cand[,1] <- lam0*exp(-(dtmp*dtmp)/(2*theta[1]*theta[1]) )
                lam.cand[,2] <- lam0*exp(-(dtmp*dtmp)/(2*theta[2]*theta[2]) ) 

	lam.mat<-matrix(lam[i,, sex[i]], nrow=R, ncol=T) * matrix(w[i], nrow=R, ncol=T) *Eff[i,,]
	lam.cand.mat<-matrix(lam.cand[,sex[i]], nrow=R, ncol=T) * matrix(w[i], nrow=R, ncol=T) *Eff[i,,]

                ll <- sum(dbinom(Z[i,,],1, lam.mat, log=TRUE))
                llcand <- sum(dbinom(Z[i,,],1, lam.cand.mat, log=TRUE))

                if(runif(1) < exp(llcand - ll)) {
                    S[i,] <- Scand
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


















