wolvSCR0Dssv2 <-function(y3d,traps,nb=1000,ni=2000,M=200,
Sgrid=wolverine$grid8,engine="bugs",area=8){


ssarea<- (nrow(Sgrid)*area)/100

# trapping grid scaled appropriately
traplocs<-as.matrix(traps[,2:3])
mingridx<-min(traplocs[,1])
mingridy<-min(traplocs[,2])
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 ###units of 10 km
ntraps<- nrow(traplocs)

## Standardize the state-space grid in the same manner
Sgrid[,1]<-Sgrid[,1]-mingridx
Sgrid[,2]<-Sgrid[,2]-mingridy
Sgrid<-Sgrid/10000 # units of 10 km
plot(Sgrid,pch=20)
points(traplocs,pch=20,cex=5,col="red")
Dmat<-e2dist(traplocs,Sgrid)
dist2<- t(Dmat)^2

### ARRAY having dimensions individual x trap x rep
## MASK is trap x rep
Y<-y3d
nz<-M-dim(Y)[1]
MASK<-traps[,4:ncol(traps)]
Dmat<-as.matrix(dist(traplocs))
nind<-dim(Y)[1]
K<-dim(Y)[3]

newy<-array(0,dim=c(nind+nz,ntraps, K))
for(j in 1:nind){
newy[j,1:ntraps,1:K]<-Y[j,1:ntraps,1:K]
}
Y<-newy
ndays<-apply(MASK,1,sum)
ncaps<- apply(Y,c(1,2),sum)

start.time = Sys.time()



sink("modelfile.txt")
cat("
model {

for(g in 1:ngrid){
 probs[g]<- 1/ngrid
}

alpha1 ~ dnorm(0,.1)

loglam0~dnorm(0,.01)
lam0<-exp(loglam0)
p0<-exp(loglam0)/(1+exp(loglam0))
sigma<- sqrt(1/(2*alpha1))
psi ~ dunif(0,1)

for(i in 1:M){
 w[i]~dbern(psi)
 s[i]~dcat(probs[1:ngrid])
for(j in 1:ntraps){
  mu[i,j]<-w[i]*p[i,j]
 ncaps[i,j]~ dbin(mu[i,j],ndays[j])
p[i,j] <- p0*exp(-alpha1*dist2[s[i],j] )
##dist2[i,j]<-  pow(Sgrid[s[i],1] - traplocs[j,1],2)   + pow(Sgrid[s[i],2] - traplocs[j,2],2)
}
}


N<-sum(w[1:M])
D<-N/ssarea
}
",fill=TRUE)
sink()

wst<-c(rep(1,M))
ngrid<-nrow(Sgrid)

data <- list ("ncaps","traplocs","M","Sgrid","ntraps","ndays","ngrid","ssarea")
sst<-rep(NA,M)
for(i in 1:nind){
    tmp0<-ncaps[i,]>0
    tmp<- dist2[,tmp0]
    if(is.matrix(tmp))
      mn<-apply(tmp,1,min)
      else
      mn<- tmp
    sst[i]<- (1:ngrid)[mn==min(mn)]

}
sst[(nind+1):M]<-sample(1:ngrid,M-nind,replace=TRUE)

inits <- function(){
  list (alpha1=runif(1,1.2,1.6),loglam0=runif(1,-3,-2),w=wst,s=sst,psi=runif(1))
}
nc<-3
nthin<-1

if(engine=="bugs"){
library("R2WinBUGS")
data <- list ("dist2","ncaps","M","ntraps","ndays","ngrid","ssarea")
parameters <- c("psi","sigma","lam0","p0","N","D","alpha1")

out <- bugs (data, inits, parameters, "modelfile.txt", n.thin=nthin,n.chains=nc,
n.burnin=nb,n.iter=ni,debug=FALSE,DIC=FALSE)
}
if(engine=="jags"){
library("rjags")
data <- list (dist2=dist2,ncaps=ncaps,M=M,ntraps=ntraps,ndays=ndays,ngrid=ngrid,ssarea=ssarea)
parameters <- c("psi","sigma","lam0","p0","N","D","alpha1") ##,"s","w")

jm<- jags.model("modelfile.txt", data=data, inits=inits, n.chains=nc,
                 n.adapt=nb)
tm<-unix.time(
jout<- coda.samples(jm, parameters, n.iter=ni-nb, thin=nthin))
out<-list(jout=jout,jm=jm)
}
return(out)
}
