SCR0poisSSp <-
function(y,X,M, SSp, delta, niter) { #SSp is the polygon of the state-space

#create initial values
nSSp<-as.owin(SSp)
S<-cbind( (runifpoint(n=M, win=nSSp))[[3]]/1000, (runifpoint(n=M, win=nSSp))[[4]]/1000 )
sigma<-runif(1,0.5,5)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
seen <- apply(y>0, 1, any)
z[seen]<-1#set seen individuals' z=1

#initiate distance matrix and lamij matrix
d <- e2dist1(S, X)
lam<-lam0*exp(-(d*d)/(2*sigma*sigma))

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
            lam.cand <- lam0*exp(-(d*d)/(2*sig.cand*sig.cand))
            ll<- sum(dpois(y, lam*z, log=TRUE))
            llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
            if(runif(1) < exp( llcand  - ll) ){
                ll<-llcand
                lam<-lam.cand
                sigma<-sig.cand
            }
        }

#update lam0
        lam0.cand <- rnorm(1, lam0, delta[2])
        if(lam0.cand>0){   #automatically reject lam0.cand that are <0
            lam.cand <- lam0.cand*exp(-(d*d)/(2*sigma*sigma))
            ll<- sum(dpois(y, lam*z, log=TRUE))
            llcand <- sum(dpois(y, lam.cand*z, log=TRUE))
            if(runif(1) < exp( llcand  - ll) ){
                ll<-llcand
                lam<-lam.cand
                lam0<-lam0.cand
            }
        }

#update z
        zUps <- 0#set counter to monitor acceptance rate
        for(i in 1:M) {
            if(seen[i])#no need to update seen individuals, since their z =1
                next
            zcand <- ifelse(z[i]==0, 1, 0)
            llz <- sum(dpois(y[i,],lam[i,]*z[i], log=TRUE))
            llcand <- sum(dpois(y[i,], lam[i,]*zcand, log=TRUE))

            prior <- dbinom(z[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(zcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (llz+prior) )) {
                z[i] <- zcand
                zUps <- zUps+1
            }
        }#end M loop

#update psi
psi<-rbeta(1, 1+sum(z), 1 + M-sum(z)) 

#update s
       Sups <- 0
        Scand <- as.matrix(cbind(rnorm(M, S[,1], delta[3]), rnorm(M, S[,2], delta[3])) )
Scoord<-SpatialPoints(Scand*1000)       #convert to spatial points on UTM (m) scale
SinPoly<-over(Scoord,SSp)# check if scand is within the polygon

        for(i in 1:M) {   

if(is.na(SinPoly[i])==FALSE) {

                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
                lam.cand<- lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma) )

                llS <- sum(dpois(y[i,], lam[i,]*z[i], log=TRUE)) 
                llcand <- sum(dpois(y[i,], lam.cand*z[i], log=TRUE)) 
                if(runif(1) < exp(llcand - llS)) {
                    S[i,] <- Scand
                    lam[i,] <- lam.cand
                    d[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
}#end M loop

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

