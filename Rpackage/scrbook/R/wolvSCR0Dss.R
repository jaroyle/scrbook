wolvSCR0Dss <-
function(y3d,traps,nb=1000,ni=2000,M=200,
Sgrid=wolverine$grid8,engine="bugs",area=8){

### GET DATA
ssarea<- (nrow(Sgrid)*area)/100
print(ssarea)


# trapping grid scaled appropriately
traplocs<-as.matrix(traps[,1:2])
mingridx<-min(traplocs[,1])
mingridy<-min(traplocs[,2])
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 ###units of 10 km
## set the state-space
ntraps<- nrow(traplocs)

### ARRAY having dimensions individual x rep x trap
## MASK is trap x rep

Y<-y3d
nz<-M-dim(Y)[1]
MASK<-traps[,3:ncol(traps)]
Dmat<-as.matrix(dist(traplocs))
nind<-dim(Y)[1]
K<-dim(Y)[2]

newy<-array(0,dim=c(nind+nz,K,ntraps))
for(j in 1:nind){
newy[j,1:K,1:ntraps]<-Y[j,1:K,1:ntraps]
}
Y<-newy
M<-nind+nz
ndays<-apply(MASK,1,sum)
ncaps<- apply(Y,c(1,3),sum)

start.time = Sys.time()

Sgrid[,1]<-Sgrid[,1]-mingridx
Sgrid[,2]<-Sgrid[,2]-mingridy
Sgrid<-Sgrid/10000 # units of 10 km
plot(Sgrid,pch=20)
points(traplocs,pch=20,cex=5,col="red")
Dmat<-e2dist(traplocs,Sgrid)
probs<-rep(1/nrow(Sgrid),nrow(Sgrid))


sink("modelfile.txt")
cat("
model {

for(g in 1:ngrid){
 probs[g]<- 1/ngrid
}
beta~dnorm(0,.1)
###sigma~dunif(0,50)
loglam0~dnorm(0,.01)
lam0<-exp(loglam0)
p0<-exp(loglam0)/(1+exp(loglam0))
sigma<- sqrt(1/(2*beta))
psi ~ dunif(0,1)
#p0~dunif(0,1)
#lam0<-log(p0/(1-p0))
for(i in 1:M){
 w[i]~dbern(psi)
 s[i]~dcat(probs[1:ngrid])
for(j in 1:ntraps){
  mu[i,j]<-w[i]*p[i,j]
 ncaps[i,j]~ dbin(mu[i,j],ndays[j]) 
p[i,j] <- p0*exp(-beta*dist2[i,j] )
dist2[i,j]<-  pow(Sgrid[s[i],1] - traplocs[j,1],2)   + pow(Sgrid[s[i],2] - traplocs[j,2],2) 
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
 sst[i]<-(1:ntraps)[ncaps[i,]>0][1]
}
sst[(nind+1):M]<-sample(1:ngrid,M-nind,replace=TRUE)
##sst<-sample(1:ngrid,M,replace=TRUE)
inits <- function(){
  list (beta=runif(1,1.2,1.6),loglam0=runif(1,-3,-2),w=wst,s=sst,psi=runif(1))
}
nc<-3
nthin<-1

if(engine=="bugs"){
library("R2WinBUGS")
data <- list ("ncaps","traplocs","M","Sgrid","ntraps","ndays","ngrid","ssarea")
parameters <- c("psi","sigma","lam0","p0","N","D") ###"x0g","y0g","total.exposure","Nin")

out <- bugs (data, inits, parameters, "modelfile.txt", n.thin=nthin,n.chains=nc, 
n.burnin=nb,n.iter=ni,debug=FALSE,DIC=FALSE)
}
if(engine=="jags"){
library("rjags")
data <- list (ncaps=ncaps,traplocs=traplocs,M=M,Sgrid=Sgrid,ntraps=ntraps,ndays=ndays,ngrid=ngrid,ssarea=ssarea)
parameters <- c("psi","sigma","lam0","p0","N","D") ##,"s","w") ###"x0g","y0g","total.exposure","Nin")

jm<- jags.model("modelfile.txt", data=data, inits=inits, n.chains=nc,
                 n.adapt=nb)
tm<-unix.time(
jout<- coda.samples(jm, parameters, n.iter=ni-nb, thin=nthin))
out<-list(jout=jout,jm=jm)
}
return(out)
}
