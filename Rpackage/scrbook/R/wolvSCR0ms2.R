wolvSCR0ms2 <-
function(nb=1000,ni=2000,buffer=2,M=200){
library("R2WinBUGS")
library("R2jags")
library("scrbook")

data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
wsex<-wolverine$wsex



# trapping grid scaled appropriately
traplocs<-as.matrix(traps[,2:3])
mingridx<-min(traplocs[,1])
mingridy<-min(traplocs[,2])
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 ###units of 10 km
ntraps<- nrow(traplocs)
## set the state-space
Xl<-min(traplocs[,1] - buffer)
Xu<-max(traplocs[,1] + buffer)
Yl<-min(traplocs[,2] - buffer)
Yu<-max(traplocs[,2] + buffer)
area<- (Xu-Xl)*(Yu-Yl)/10

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


sink("modelfile5.txt")
cat("
model {

for(i in 1:3){
alpha0[i] ~ dnorm(0,.1)
sigma[i] ~ dunif(0,5)
beta[i]<- 1/(2*sigma[i]*sigma[i])
}

psi ~ dunif(0,1)
mod~  dcat(mod.probs[])
mod.probs[1]<- 1/3
mod.probs[2]<- 1/3
mod.probs[3]<- 1/3

for(i in 1:M){
 w[i]~dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)

logit(p0[i,1])<- alpha0[1]
log(p0[i,2])<- alpha0[2]
p0[i,3]<- alpha0[3]
beta.vec[i,1]<- beta[1]
beta.vec[i,2]<- beta[2]
beta.vec[i,3]<- beta[3]

for(j in 1:ntraps){

 mu[i,j]<-w[i]*p[i,j,mod]
 y[i,j]~ dbin(mu[i,j],K[j])

 dist2[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2)
 p[i,j,1]       <-  p0[i,1]*exp( - beta.vec[i,1]*dist2[i,j] )
 p[i,j,2]       <-  1-exp(-p0[i,2]*exp( - beta.vec[i,2]*dist2[i,j] ) )
 logit(p[i,j,3])<- p0[i,3] - beta.vec[i,3]*dist2[i,j]

}

}


N<-sum(w[1:M])
D<-N/area
}
",fill=TRUE)
sink()


sst<-traplocs[sample(1:nrow(traplocs),M,replace=TRUE),]
for(i in 1:nind){
sst[i,]<- matrix(traplocs[y[i,]>0,],ncol=2,byrow=FALSE)[1,]

}

wst<-c(rep(1,nind),rep(0,M-nind))
mod<-1
data <- list ("y","traplocs","M","ntraps","K","Xl","Xu","Yl","Yu","area")

if(engine=="jags")  {
inits <- function(){   list (psi=.5,alpha0=rnorm(3,-2,.5),beta=runif(3,.2,.8),w=wst,s=sst,mod=1) }
parameters <- c("mod","psi","sigma","beta","alpha0","N","D")
out <-  jags(data, inits, parameters, "mf5.txt", n.thin=1,n.chains=3, n.burnin=nb,
n.iter=ni,working.dir=getwd())
}

if(engine=="bugs"){
inits <- function(){   list (psi=.5,alpha0=rnorm(3,-2,.5),beta=runif(3,.2,.8),w=wst,s=sst,mod=2) }
parameters <- c("mod","psi","sigma","beta","alpha0","N","D")
out <-  bugs(data, inits, parameters, "mf5.txt", n.thin=1,n.chains=3, n.burnin=nb,
n.iter=ni,working.dir=getwd(),DIC=FALSE)
}


out



}
