

library("scrbook")
library("secr")
data(ovenbird)
set.seed(2013)

png("Ovenbird_traps.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,3))
 plot(ovenCH[["2005"]])
 plot(ovenCH[["2007"]])
 plot(ovenCH[["2009"]])
dev.off()

## extract the trap locations and create a state-space by adding 150 m
X<-traps<-traps(ovenCH)
xlim<-c(min(X[[1]][,1])-200,max(X[[1]][,1])+200)
ylim<-c(min(X[[1]][,2])-200,max(X[[1]][,2])+200)
ntraps<- nrow(traps[[1]])

## Y are the encounter history data
Y<-ovenCH
K<-10  # number of samples in each year
M<-100 # do constant data augmentation to all years

## starting values for each individual's activity centers
Sst0<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
Sst<-NULL
Ymat<-NULL
died<-NULL
# make the data into a 3-d array

for(i in 1:5){
tmp<-Y[[i]]
nind<-nrow(tmp)
nrep<-ncol(tmp)
D<-matrix(0,nrow=M,ncol=10)
D[1:nind,1:nrep]<-tmp
D[D>0]<- 0
D[D<0]<- 1
died<-rbind(died,D)

tmp[tmp<0]<- 0 ## dead guy set to 0 b/c "died" created above
tmp[tmp==0]<-ntraps+1
tmp2<-matrix(NA,nrow=M,ncol=10)  # pad last col with NA for year 1
tmp2[,1:nrep]<-ntraps+1
tmp2[1:nind,1:nrep]<-tmp
Ymat<- rbind(Ymat, tmp2)

sout<-spiderplot(tmp2[1:nind,1:nrep],as.matrix(X[[i]]))$avg.s
Stmp<-Sst0
Stmp[1:nind,1:2]<-sout
Sst<-rbind(Sst,Stmp)
}
for(i in 1:nrow(died)){
xx<-died[i,]
if(sum(xx)>0){
first<-(1:length(xx))[xx==1]
died[i,first:ncol(died)]<-1
died[i,first]<-0
}
}

## 
## This bit of code dumps out the BUGS model file
##
cat("
model {
 # year-specific N parameterized in DA parameter 

alpha0 ~ dnorm(0,.01)
sigma ~ dunif(0,200)
alpha1<- 1/(2*sigma*sigma)

A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
for(t in 1:5){
N[t] <- inprod(z[1:bigM],yrdummy[,t])
D[t] <- (N[t]/A)*10000  # put in units of per ha
psi[t] ~ dunif(0,1)
}

for(i in 1:bigM){

  z[i] ~ dbern(psi[year[i]])
  S[i,1] ~ dunif(xlim[1],xlim[2])
  S[i,2] ~ dunif(ylim[1],ylim[2])

 for(j in 1:ntraps){
    d2[i,j] <- pow(pow(S[i,1]-X[j,1],2) + pow(S[i,2]-X[j,2],2),1)
  }
for(k in 1:K){
 Ycat[i,k] ~ dcat(cp[i,k,])

  for(j in 1:ntraps){
    lp[i,k,j] <- exp(alpha0 - alpha1*d2[i,j])*z[i]*(1-died[i,k])        
    cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,1:ntraps]))
  }
  cp[i,k,ntraps+1] <- 1-sum(cp[i,k,1:ntraps])  # last cell = not captured

}
}
   
}
",file="model.txt")
###
###
###

nind<-c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]))
yrid<-sort(rep(1:5,M))
yrdummy<-as.numeric(c(yrid==1,yrid==2,yrid==3,yrid==4,yrid==5))
yrdummy<-matrix(yrdummy,ncol=5,byrow=FALSE)

bigM<- 5*M
## starting values for data augmentation parameters
zst<-NULL
for(i in 1:5){
zst<-c(zst,rep(1,nind[i]),rep(0,M-nind[i]))
}

inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst,alpha0=runif(1,-2,-1) 
,alpha2=-2,alpha3=-2) }              
## parameters to monitor
parameters <- c("psi","alpha0","alpha1","sigma","N","D")
## data used in BUGS model                                                                  
data <- list (died=died,yrdummy=yrdummy,year=yrid,X=as.matrix(X[[1]]),K=10,Ycat=Ymat,bigM=bigM,ntraps=ntraps,ylim=ylim,xlim=xlim)         

##
### This takes ~1 hour or so to run
##
library("R2jags")
out.200 <- jags(data, inits, parameters, "model.txt", n.thin=1,n.chains=3, 
n.burnin=1000,n.iter=5000,DIC=FALSE)




###
###
### Now fit some models in secr
###
###
###

## fit constant-density model
ovenbird.model.1 <- secr.fit(ovenCH,buffer=300)
## fit net avoidance model
ovenbird.model.1b <- secr.fit(ovenCH, model =   list(g0~b),buffer=300)
## fit model with time trend in detection
ovenbird.model.1T <- secr.fit(ovenCH, model =    list(g0 ~ T),buffer=300)
## fit model with 2-class mixture for g0
ovenbird.model.h2 <- secr.fit(ovenCH, model =    list(g0~h2),buffer=300)
## fit a model with session (year)-specific Density
ovenbird.model.D <- secr.fit(ovenCH, model =    list(D~session),buffer=300)



## compare & average pre-fitted models
AIC (ovenbird.model.1, ovenbird.model.1b, ovenbird.model.1T,
    ovenbird.model.h2,ovenbird.model.D)



