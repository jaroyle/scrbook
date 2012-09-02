
# MCMC with Z partially observed for raccoon data set

raccoon.MCMC <- function(n,X,Eff,y,M, niters,SSp,collar, locs,mi,mall, inits) {  #mi = all identified marked pics per occasion; mall=all marked pics per occasion

    ntot<-length(locs)  #number of known overall
    Sk<-   t(sapply(locs,colMeans))  #HR centers of known
    theta<-inits$theta 
    lam0 <- inits$lam0  #vector of 6
    prat<-mi/mall
    Dk<-e2dist(Sk, X)
  
    R <- nrow(n)
    T <- ncol(n)

lamk<-array(NA, c(ntot,R,T))  
for (t in 1:T){
    lamk [,,t]<- lam0[t]*exp(-(Dk*Dk)/(2*theta*theta))
}

    S <- cbind(inits$xz$x/1000,inits$xz$y/1000) 
    S[collar==TRUE,]<-Sk

#start of known photographed uncollared at the camera they were photographed at
iv<-(1:104)[apply(y,1,sum)>0]
lc<-NULL
for (i in 1:length(iv)){
lc[i]<-(1:R)[apply(y[iv[i],,],1,sum)>0]
S[iv[i],]<-X[lc[i],]
}


    D <- e2dist(S, X)
  
lam<-array(NA, c(M,R,T))  
for (t in 1:T){
    lam[,,t] <- lam0[t]*exp(-(D*D)/(2*theta*theta))
}

    w <- inits$w 
    psi <- inits$psi 
    SSp<-SSp

### yr 1
    Z <- array(NA, c(M,R,T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if(!missing(y) || !is.null(y)) {
        nMarked <- nrow(y)
        marked[1:nMarked] <- TRUE
        Z[1:nMarked,,] <- y
    }
    w[marked] <- 1
    Zdata <- !is.na(Z)
    for(r in 1:R) {
        for(t in 1:T) {
            if(n[r,t]==0) {
                Z[,r,t]  <- 0
                next
            }
            unmarked <- !Zdata[,r,t]
            nUnknown <- n[r,t] - sum(Z[!unmarked,r,t])
            if(nUnknown < 0)
                browser()
            probs <- lam[,r,t]*w*Eff[,r,t]
            probs <- probs[unmarked]
            Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs) 
        }
    }


    out <- matrix(NA,nrow=niters,ncol=15)
    colnames(out) <- c("sigma", "lam0(1)","lam0(2)","lam0(3)","lam0(4)","lam0(5)","lam0(6)" ,"prat(1)","prat(2)","prat(3)","prat(4)","prat(5)","prat(6)" ,"psi","N") #

    cat("\nstarting values =", c(theta, lam0,prat, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }


########################################################################################################################
######## telemetry model to estimate sigma #############################################################################


### update sigma

theta.cand<-rnorm(1, theta, 0.03)

if(theta.cand>0) {

llsig<-llsig.cand<-rep(NA, ntot)

for (ind in 1:ntot) {
llsig[ind]<-sum(dmvnorm(x=locs[[ind]],mean=c(Sk[ind,1],Sk[ind,2]), sigma=cbind(c(theta^2,0), c(0,theta^2)), log=T))   
llsig.cand[ind]<-sum(dmvnorm(x=locs[[ind]],mean=c(Sk[ind,1],Sk[ind,2]), sigma=cbind(c(theta.cand^2,0), c(0,theta.cand^2)), log=T))   
}

   if(runif(1) < exp( sum(llsig.cand)  - sum(llsig) ) ){
theta<-theta.cand

for (t in 1:T){
    lam[,,t] <- lam0[t]*exp(-(D*D)/(2*theta.cand*theta.cand))
    lamk[,,t] <- lam0[t]*exp(-(Dk*Dk)/(2*theta.cand*theta.cand))
}
}
}


############################################################################################################################
########## camera trap model ###############################################################################################

	#update rat directly from full conditional using Gibbs sampling

for (t in 1:T) {

  prat[t]<-rbeta(1,1+mi[t],1+mall[t]-mi[t]) 

}

#make matrix
rat<-matrix(rep(prat,R), nrow=R, ncol=T, byrow=T)

        # update lam0
for ( t in 1:T) {

        lam0.cand <- rnorm(1, lam0[t], .004)

        if(lam0.cand>0) {
           
	lam.cand.arr <- matrix(w,nrow=M,ncol=R) * ( lam0.cand*exp(-(D*D)/(2*theta*theta)) ) * Eff[,,t]  #occ. 1 

	lam.arr<- matrix(w,nrow=M,ncol=R) * lam[,,t] * Eff[,,t]


	ll <- sum(dpois(Z[marked,,t], lam.arr[marked,]*rat[1,t], log=TRUE)) + sum(dpois(Z[!marked,,t], lam.arr[!marked,], log=TRUE))

	llcand <- sum(dpois(Z[marked,,t], lam.cand.arr[marked,]*rat[1,t], log=TRUE)) + sum(dpois(Z[!marked,,t], lam.cand.arr[!marked,], log=TRUE))


            if(runif(1) < exp(llcand - ll) ) { 
                lam0[t] <- lam0.cand
                lam [,,t]<- lam.cand.arr
		lamk[,,t] <- lam0.cand*exp(-(Dk*Dk)/(2*theta*theta))
            }
        }
}

        ### update "w" here
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i] | marked[i])
                next
            wcand <- ifelse(w[i]==0, 1, 0)

	lam.mat<- matrix(w[i], nrow=R,ncol=T) * lam[i,,] * Eff[i,,]
	lam.mat.cand<- matrix(wcand, nrow=R,ncol=T) * lam[i,,]* Eff[i,,]

            ll <- sum(dpois(Z[i,,],lam.mat, log=TRUE))
            llcand <- sum(dpois(Z[i,,], lam.mat.cand, log=TRUE))

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand
                wUps <- wUps+1
            }
        }


        # update Z
        for(r in 1:R) {
            for(t in 1:T) {
		zip <- lam[,r,t]*w *Eff[,r,t]

	#zipE<-zip*Eff[,r,t]
            #probs <- zipE/sum(zipE)

                if(n[r,t]==0) {
                    Z[,r,t] <- 0
                    next
                }
                unmarked <- !Zdata[,r,t]
                nUnknown <- n[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown < 0)
                    browser()
		if(nUnknown == 0)
   next
                Z[unmarked,r,t] <- rmultinom(1, nUnknown, zip[unmarked])
            }
        }

#dput(Z, file = paste('Z_', iter, '.R', sep=''))

        # update psi
        psi<-rbeta(1,1+sum(w[!marked]),1+sum(!marked)-sum(w[!marked]))  #update psi according to new values of w



### update Sk from telemetry data
	Skups<-0
        Skcand <- as.matrix(cbind(rnorm(ntot, Sk[,1], 0.5),
                   rnorm(ntot, Sk[,2], 0.05)))

	#SkcandX<-rnorm(ntot, Sk[,1], 0.5)
	#ak<-3338.222+SkcandX*1.373+0.6
	#bk<-3338.222+SkcandX*1.373 -0.6
        #Skcand <- as.matrix(cbind(SkcandX,
         #         runif(ntot,bk,ak)))

	Scoord<-SpatialPoints(Skcand*1000)
	SinPoly<-over(Scoord,SSp)

        for(i in 1:ntot) {   # have to be attributed to the right S's
	
	if(is.na(SinPoly[i,1])==FALSE) {

	llsk<-sum(dmvnorm(x=locs[[i]],mean=c(Sk[i,1],Sk[i,2]), sigma=cbind(c(theta^2,0), c(0,theta^2)), log=T))   
	llsk.cand<-sum(dmvnorm(x=locs[[i]],mean=c(Skcand[i,1],Skcand[i,2]), sigma=cbind(c(theta^2,0), c(0,theta^2)), log=T))   

                if(runif(1) < exp(llsk.cand - llsk)) {
		Sk[i,]<-Skcand[i,]
		Dk[i,]<-sqrt((Sk[i,1] - X[,1])^2 + (Sk[i,2] - X[,2])^2)
		for (t in 1:T){
		lamk[i,,t]<-lam0[t]*exp(-(Dk[i,]*Dk[i,])/(2*theta*theta) )
				}
		Skups<-Skups+1
		}
			}
		}


        # update S for unknown
        Sups <- 0
	ScandX<-rnorm(M, S[,1], 1.5)
	a<-3338.222+ScandX*1.373+0.6
	b<-3338.222+ScandX*1.373 -0.6
        Scand <- as.matrix(cbind(ScandX,
                  runif(M,b,a)))


	Scoord<-SpatialPoints(Scand*1000)
	SinPoly<-over(Scoord,SSp)

        for(i in 1:M) {   # note this is "M" in general

	if (collar[i]==TRUE)   #skip collared

		next

	if(is.na(SinPoly[i,1])==FALSE) {

                dtmp <- sqrt((Scand[i,1] - X[,1])^2 + (Scand[i,2] - X[,2])^2)

		lam.cand<-matrix(NA, nrow=R, ncol=T)
		for (t in 1:T){
                lam.cand[,t] <- lam0[t]*exp(-(dtmp*dtmp)/(2*theta*theta) )
			}

                #create lam matrix

		lam.mat<- matrix(w[i], nrow=R,ncol=T) * lam[i,,]* Eff[i,,]
		lam.mat.cand<- matrix(w[i], nrow=R,ncol=T) * lam.cand* Eff[i,,]

			if (marked[i]){
			
                ll <- sum(dpois(Z[i,,], lam.mat*rat, log=TRUE)) 
                llcand <- sum(dpois(Z[i,,], lam.mat.cand*rat, log=TRUE)) 
					} else {
			
                ll <- sum(dpois(Z[i,,], lam.mat, log=TRUE)) 
                llcand <- sum(dpois(Z[i,,], lam.mat.cand, log=TRUE)) 
					} 

                if(runif(1) < exp(llcand - ll)) {
                    ll <- llcand
                    S[i,] <- Scand[i,]
                    lam[i,,] <- lam.cand
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
	}

### substitute S for known guys with current activity center from telemetry model, adjust D and lam

	S[collar==TRUE,]<-Sk 
	D[collar==TRUE,]<-Dk
	lam[collar==TRUE,,]<-lamk



        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/(M-nMarked), "\n")
            cat("     S =", Sups/(M-ntot), "\n")
            cat("     Sk =", Skups/ntot, "\n")
        }

        out[iter,] <- c(theta,lam0,prat, psi,sum(w))

    }


   # return(list(out=out, last=list(theta=theta,lam0=lam0,psi1=psi1, psi2=psi2, Z1=Z1,Z2=Z2,w1=w1,w2=w2, S1=S1, S2=S2)))
return(out)
}


















