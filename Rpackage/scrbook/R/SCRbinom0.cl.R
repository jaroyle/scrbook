SCR0binom.cl <-
function(y,X,M, xl,xu,yl,yu,K,delta, niter) {

#create initial values
S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
sigma<-runif(1,5,10)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
seen <- apply(y>0, 1, any)
z[seen]<-1		#set seen individuals' z=1
K=K
delta=delta

#initiate distance matrix and lamij matrix
D <- e2dist(S, X)
lam<-lam0*exp(-(D*D)/(2*sigma*sigma))
pmat<- 1-exp(-lam)
pmat[pmat<1e-30] <- 1e-30 #ensures >0 detection probabilities in first iteration


#set up matrix to hold results
out<-matrix(nrow=niter, ncol=4)
colnames(out)<-c('sigma', 'lam0', 'psi', 'N')

#have R print starting values
cat("\nstarting values =", c(sigma, lam0, psi, sum(z),"\n\n"))

#start iterations of the chain

for (iter in 1:niter) {

#have R output the time and parameter estimates at every 100th iteration
        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

#update sigma
        sig.cand <- rnorm(1, sigma, delta[1])
        if(sig.cand>0){   #automatically reject sig.cand that are <0
            lam.cand <- lam0*exp(-(D*D)/(2*sig.cand*sig.cand))
	    p.cand<-1-exp(-lam.cand)

		#if (iter==1) p.cand[p.cand<1e-30] <- 1e-30 

            ll<- sum(dbinom(y, K,pmat*z, log=TRUE))
            llcand <- sum(dbinom(y,K, p.cand*z, log=TRUE))
            if(runif(1) < exp( llcand  - ll) ){
                ll<-llcand
                pmat<-p.cand
                sigma<-sig.cand
            }
        }

#update lam0
        lam0.cand <- rnorm(1, lam0,delta[2])
        if(lam0.cand>0){   #automatically reject lam0.cand that are <0
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
	    p.cand<-1-exp(-lam.cand)
            llcand <- sum(dbinom(y,K,p.cand*z, log=TRUE))
            if(runif(1) < exp( llcand  - ll) ){
                pmat<-p.cand
                lam0<-lam0.cand
            }
        }

#update z
        zUps <- 0		#set counter to monitor acceptance rate
        for(i in 1:M) {
            if(seen[i])	#no need to update seen individuals, since their z =1
                next
            zcand <- ifelse(z[i]==0, 1, 0)
            llz <- sum(dbinom(y[i,],K,pmat[i,]*z[i], log=TRUE))
            llcand <- sum(dbinom(y[i,],K, pmat[i,]*zcand, log=TRUE))

            prior <- dbinom(z[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (llz+prior) )) {
                z[i] <- zcand
                zUps <- zUps+1
            }
        }	#end M loop

#update psi
	psi<-rbeta(1, 1+sum(z), 1 + M-sum(z)) #is that correct? or for this to be correct, do we have to update all z's and not set them =1 for the seen individuals??

#update s
       Sups <- 0
        for(i in 1:M) {   
        Scand <- c(rnorm(1, S[i,1], delta[3]), rnorm(1, S[i,2], delta[3])) 
		inbox<-Scand[1]<xu & Scand[1]>xl & Scand[2]<yu & Scand[2]>yl
		if(inbox){
                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
                lam.cand<- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )
		p.cand<-1-exp(-lam.cand)

                llS <- sum(dbinom(y[i,],K, pmat[i,]*z[i], log=TRUE)) 
                llcand <- sum(dbinom(y[i,],K, p.cand*z[i], log=TRUE)) 
                if(runif(1) < exp(llcand - llS)) {
                    S[i,] <- Scand
                    pmat[i,] <- p.cand
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
}	#end M loop

#prompt R to output acceptance rates of z and S
        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     z =", zUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(sigma,lam0,psi,sum(z))

}  #end of iteration loop

    return(out)

}  #end of function call

