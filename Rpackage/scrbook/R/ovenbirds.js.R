
ovenbirds.js<-function(model=c('NSJS', 'SJS', 'SMS'), n.chains, n.adapt, n.iter) {

mod<-match.arg(model)

#install required packages
library(rjags)
library(scrbook)
library(secr)

##data setup and organization
data(ovenbird)

##set up the traps and state space
X<-traps<-traps(ovenCH)
xlim<-c(min(X[[1]][,1])-150,max(X[[1]][,1])+150)
ylim<-c(min(X[[1]][,2])-150,max(X[[1]][,2])+150)
ntraps<- nrow(traps[[1]])
Y<-ovenCH
K<-10
# data augmentation to all years
M<-200 
Sst<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
Sst<-array(Sst,dim=c(M,2,5))


#Link individuals across years
hold<- unique(c(unlist(dimnames(Y[[1]])[1]), unlist(dimnames(Y[[2]])[1]), unlist(dimnames(Y[[3]])[1]), unlist(dimnames(Y[[4]])[1]), unlist(dimnames(Y[[5]])[1])))

Yarr<-array(ntraps+1,dim=c(M,K,5))
for(i in 1:5){
tmp<-Y[[i]]
tmp[tmp<0]<-tmp[tmp<0]*(-1) ## one guy died, we ignore that here and treat it as a  normal event
tmp[tmp==0]<-ntraps+1
nind<-nrow(tmp)
nrep<-ncol(tmp)
tmp2<-matrix(ntraps+1,nrow=M,ncol=10)  # pad last col with NA for year 1
tmp2[pmatch(unlist(dimnames(Y[[i]])[1]), hold),1:nrep]<-tmp

#set inital values for activity centers
Stmp<-Sst[,,i]
Stmp[pmatch(unlist(dimnames(Y[[i]])[1]), hold),1:2]<-spiderplot(tmp2[pmatch(unlist(dimnames(Y[[i]])[1]), hold),1:nrep],as.matrix(X[[i]]))$avg.s
Sst[,,i]<-Stmp
Yarr[,,i]<-tmp2
}

#######################################################################################################
###### Spatial Jolly Seber  SJS #######################################################################
if(mod=='SJS') {

cat("
model {
#PRIORS

psi ~ dunif(0,1)
phi ~ dunif(0,1)
alpha0 ~ dnorm(0,10)
sigma ~dunif(0,200)
alpha1<- 1/(2*sigma*sigma)

A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))

for(t in 1:5){
N[t] <- sum(z[1:M,t])
D[t] <- N[t]/A
gamma[t] ~ dunif(0,1)

}

for(i in 1:M){
  z[i,1] ~ dbern(psi)
    
  #to estimate the number of recruits, we need a few derivations 
  R[i,1]<- z[i,1]
  R[i,2]<-(1-z[i,1])*z[i,2]
  R[i,3]<- (1-z[i,1])*(1-z[i,2])*z[i,3]
  R[i,4] <-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*z[i,4])
  R[i,5] <-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*z[i,5]



  for(t in 1:5){

  S[i,1,t] ~ dunif(xlim[1],xlim[2])
  S[i,2,t] ~ dunif(ylim[1],ylim[2])

  for(j in 1:ntraps){
    d[i,j,t] <- pow(pow(S[i,1,t]-X[j,1],2) + pow(S[i,2,t]-X[j,2],2),1)
     }

  for(k in 1:K){
    for(j in 1:ntraps){
      lp[i,k,j,t] <- exp(alpha0 - alpha1*d[i,j,t])*z[i,t]           
      cp[i,k,j,t] <- lp[i,k,j,t]/(1+sum(lp[i,k,,t]))
    }
    cp[i,k,ntraps+1,t] <- 1-sum(cp[i,k,1:ntraps,t])  # last cell = not captured
    Ycat[i,k,t] ~ dcat(cp[i,k,,t])
  } 
}  

a[i,1]<-(1-z[i,1])

for(t in 2:T){                           #ensure that individuals are not recruited into the population more than once
      a1[i,t] <- sum(z[i, 1:t])          #have you ever been alive (0 = no, >1 = yes)
       a[i,t] <- 1-step(a1[i,t] - 1)     #use the step function to make a1 binary 

       mu[i,t]<- (phi*z[i,t-1]) + (gamma[t]*a[i,t-1])
        z[i,t]~dbern(mu[i,t])
          }
        }

R1<-sum(R[1:M,1])
R2<-sum(R[1:M,2])
R3<-sum(R[1:M,3])
R4<-sum(R[1:M,4])
R5<-sum(R[1:M,5])        
}

",file="modelJS.txt")
###
###
###
###

 

zst<-c(rep(1,M/2),rep(0,M/2))
zst<-cbind(zst,zst,zst,zst,zst)

inits <- function(){list (z=zst,sigma=runif(1,25,100), gamma=runif(5,0,1) ,S=Sst,alpha0=runif(1,-2,-1) ) }             

parameters <- c("psi","alpha0","alpha1","sigma","N","D", "phi", "gamma", "R2", "R3", "R4", "R5")
                                                                  
data <- list (X=as.matrix(X[[1]]),K=10,T=5, Ycat=Yarr,M=M,ntraps=ntraps,ylim=ylim,xlim=xlim)        

modelFile= "modelJS.txt"

} #end if

######################################################################################################
###### Non-Spatial Jolly Seber #####################################################################


if(mod=='NSJS') {

cat("
model {
#PRIORS

psi ~ dunif(0,1)
phi ~ dunif(0,1)
p.mean ~ dunif(0,1)


for(t in 1:T){
N[t] <- sum(z[1:M,t])
gamma[t] ~ dunif(0,1)
}

for(i in 1:M){
  z[i,1] ~ dbern(psi)            #Alive state for the first year
 cp[i,1] <- z[i,1]*p.mean
  Y[i,1] ~ dbinom(cp[i,1], K)    #Y are the number of encounters for the secondary occasions
  A[i,1]<-(1-z[i,1])

  #to estimate the number of recruits, we need a few derivations 
  R[i,1]<- z[i,1]
  R[i,2]<-(1-z[i,1])*z[i,2]
  R[i,3]<- (1-z[i,1])*(1-z[i,2])*z[i,3]
  R[i,4] <-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*z[i,4]
  R[i,5] <-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*z[i,5]


for(t in 2:T){    #for loop over years 2 to T
      a1[i,t] <- sum(z[i, 1:t])   #have you ever been alive (0 = no, >1 = yes)
       A[i,t] <- 1-step(a1[i,t] - 1)     #use the step function to make a1 binary 
  #A is the indicator if an individual is available to be recruited
      mu[i,t]<- (phi*z[i,t-1]) + (gamma[t]*A[i,t-1])
       
       z[i,t]~dbern(mu[i,t])
      cp[i,t] <- z[i,t]*p.mean
       Y[i,t] ~ dbinom(cp[i,t], K)                  
          }
        }

R1<-sum(R[1:M,1])
R2<-sum(R[1:M,2])
R3<-sum(R[1:M,3])
R4<-sum(R[1:M,4])
R5<-sum(R[1:M,5])        
}

",file="modelNSJS.txt")
###
###

#set up data for non-spatial model
Yarr[Yarr < 45] <- 1
Yarr[Yarr == 45] <- 0
Ybin=matrix(NA, M, 5)
for(t in 1:5){
Ybin[,t] <- rowSums(Yarr[,,t])
}

zst<-c(rep(1,M/2),rep(0,M/2))
zst<-cbind(zst,zst,zst,zst,zst)

inits <- function(){list (z=zst, gamma=runif(5,0,1), phi=runif(1,0,1), p.mean=runif(1,0,1) ) }             

parameters <- c("psi","p.mean","N", "phi", "gamma", "R2", "R3", "R4", "R5")
                                                                  
data <- list (Y=Ybin, M=M, T=5, K=10)

modelFile= "modelNSJS.txt"

} #end if

######################################################################################################
###### Spatial Multi-Session SMS #####################################################################


if(mod=='SMS') {

cat("
model {
#PRIORS
for(t in 1:5){
psi[t] ~ dunif(0,1)
}
alpha0 ~ dnorm(0,10)
sigma ~dunif(0,200)
alpha1<- 1/(2*sigma*sigma)

A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
for(t in 1:5){
N[t] <- sum(z[1:M,t])
D[t] <- N[t]/A

for(i in 1:M){
  z[i,t] ~ dbern(psi[t])
  S[i,1,t] ~ dunif(xlim[1],xlim[2])
  S[i,2,t] ~ dunif(ylim[1],ylim[2])
  for(j in 1:ntraps){
    #distance from capture to the center of the home range
    d[i,j,t] <- pow(pow(S[i,1,t]-X[j,1],2) + pow(S[i,2,t]-X[j,2],2),1)
  }
  for(k in 1:K){
    for(j in 1:ntraps){
      lp[i,k,j,t] <- exp(alpha0 - alpha1*d[i,j,t])*z[i,t]           
      cp[i,k,j,t] <- lp[i,k,j,t]/(1+sum(lp[i,k,,t]))
    }
    cp[i,k,ntraps+1,t] <- 1-sum(cp[i,k,1:ntraps,t])  # last cell = not captured
    Ycat[i,k,t] ~ dcat(cp[i,k,,t])
  } 
}  
}

}",file="modelSMS.txt")
###
###



zst<-c(rep(1,M/2),rep(0,M/2))
zst<-cbind(zst,zst,zst,zst,zst)

inits <- function(){list (z=zst,sigma=runif(1,25,100),S=Sst,alpha0=runif(1,-2,-1) ) }             
                                                           
data <- list (X=as.matrix(X[[1]]),K=10,Ycat=Yarr,M=M,ntraps=ntraps,ylim=ylim,xlim=xlim)        

parameters <- c("psi","alpha0","alpha1","sigma","N","D")

modelFile= "modelSMS.txt"

} #end if



######################################################################################################

#run model
mod.out <- jags.model(modelFile, data, inits, n.chains=n.chains, n.adapt=n.adapt)
out <- coda.samples(mod.out,  parameters, n.iter=n.iter)


return(out)

} #end function





