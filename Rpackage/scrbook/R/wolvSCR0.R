wolvSCR0 <-
function(y3d,traps,ni=2000, nb=1000,buffer=2,M=200,keepz=FALSE,coord.scale=10000){

library("R2WinBUGS")

# trapping grid scaled appropriately
traplocs<-as.matrix(traps[,2:3])
mingridx<-min(traplocs[,1])
mingridy<-min(traplocs[,2])
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/coord.scale ###units of 10 km if coord.scale=10000
ntraps<- nrow(traplocs)
## set the state-space
Xl<-min(traplocs[,1] - buffer)
Xu<-max(traplocs[,1] + buffer)
Yl<-min(traplocs[,2] - buffer)
Yu<-max(traplocs[,2] + buffer)
area<- (Xu-Xl)*(Yu-Yl)/10
#####plot(traplocs,pch=20,xlim=c(Xl,Xu),ylim=c(Yl,Yu))
### ARRAY having dimensions individual x rep x trap
## MASK is trap x rep
nz<-M-dim(y3d)[1]
MASK<-traps[,4:ncol(traps)]
Dmat<-as.matrix(dist(traplocs))
nind<-dim(y3d)[1]
K<-dim(y3d)[2]

## Data Augmentation

newy<-array(0,dim=c(nind+nz,K,ntraps))
for(j in 1:nind){
newy[j,1:K,1:ntraps]<-y3d[j,1:K,1:ntraps]
}
y3d<-newy
M<-nind+nz
# compute trap-specific sample size
K<-apply(MASK,1,sum)
y<- apply(y3d,c(1,3),sum)

sink("modelfile.txt")
cat("
model {

sigma~dunif(0,50)
p0~dunif(0,1)
alpha1<- (1/(2*sigma*sigma) )

psi ~ dunif(0,1)
for(i in 1:M){
 z[i]~dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
}
for(i in 1:M){
  for(j in 1:ntraps){
   mu[i,j]<-z[i]*p[i,j]
   y[i,j]~ dbin(mu[i,j],K[j])
   dd[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2)
   p[i,j]  <-  p0*exp( - alpha1*dd[i,j] )
   ncapsnew[i,j]~dbin(mu[i,j],K[j])
   err[i,j]<-  pow(pow(y[i,j],.5) - pow(K[j]*p[i,j],.5),2)
   errnew[i,j]<- pow(pow(ncapsnew[i,j],.5) - pow(K[j]*p[i,j],.5),2)
 }

}
Xobs<-sum(err[,])
Xnew<-sum(errnew[,])

N<-sum(z[1:M])
D<-N/area
}
",fill=TRUE)
sink()

data <- list ("y","traplocs","M","ntraps","K","Xl","Xu","Yl","Yu","area")

sst<-cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))  # starting values for s
for(i in 1:nind){
if(sum(y[i,])==0) next
sst[i,1]<- mean( traplocs[y[i,]>0,1] )
sst[i,2]<- mean( traplocs[y[i,]>0,2] )
}


zst<-c(rep(1,nind),rep(0,M-nind))
inits <- function(){
  list (sigma=runif(1,.4,1),p0=runif(1,.01,.2),z=zst,s=sst)
}
parameters <- c("psi","sigma","p0","N","D","alpha1","s","z") ###,"Xobs","Xnew","s")
if(keepz)
parameters <- c("psi","sigma","p0","N","D","alpha1","s","z") ###,"Xobs","Xnew","s","z")

out <- bugs(data, inits, parameters, "modelfile.txt", n.thin=1,n.chains=3, n.burnin=nb,n.iter=ni,working.dir=getwd(),debug=FALSE)
out
}
