wolvSCRMb <-
function(nb=1000,ni=2000,buffer=2,M=200){
library("R2jags")
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
wsex<-wolverine$wsex

library("R2WinBUGS")

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
plot(traplocs,pch=20,xlim=c(Xl,Xu),ylim=c(Yl,Yu))
### ARRAY having dimensions individual x rep x trap
## MASK is trap x rep
nz<-M-dim(y3d)[1]
MASK<-as.matrix(traps[,4:ncol(traps)])
Dmat<-as.matrix(dist(traplocs))
nind<-dim(y3d)[1]
K<-dim(y3d)[2]

## Data Augmentation

newy<-array(0,dim=c(nind+nz,K,ntraps))
for(j in 1:nind){
newy[j,1:K,1:ntraps]<-y3d[j,1:K,1:ntraps]
}
y3d<-newy

## create a behavioral response covariate
B<-array(0,dim=dim(y3d))
for(i in 1:dim(y3d)[1]){
for(j in 1:dim(y3d)[3]){
xx<-y3d[i,,j]
if(any(xx>0)){
 fst<- (1:length(xx))[xx==1][1]
 B[i,(fst+1):length(xx),j]<-1
}
}
}

M<-nind+nz

ncaps<- apply(y3d,c(1,3),sum)
y<-y3d

sink("modelfile.txt")
cat("
model {

#beta~dnorm(0,.01)
sigma~dunif(0,50)
p0~dunif(0,1)
logitp0<-log(p0/(1-p0))
beta<- (1/(2*sigma*sigma) )
alpha2 ~ dnorm(0,.1)
psi ~ dunif(0,1)
for(i in 1:M){
 w[i]~dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu) 

ytot[i]<- sum(y2[i,])
for(j in 1:ntraps){
y2[i,j]<- sum(y[i,,j])
}
}
for(i in 1:M){
for(j in 1:ntraps){

dd[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2) 

for(k in 1:K){
logit(pbase[i,k,j])<- logitp0   + bleen[i,k,j]
bleen[i,k,j]<- alpha2*B[i,k,j]
   mu[i,k,j]<-w[i]*p[i,k,j]*MASK[j,k]
 y[i,k,j]~ dbern(mu[i,k,j]) 
  p[i,k,j]  <-  pbase[i,k,j]*exp( - beta*dd[i,j] )
ncapsnew[i,k,j]~dbern(mu[i,k,j])
}
# total over time
ncapsnew2[i,j]<- sum(ncapsnew[i,,j])
expected[i,j]<- sum(mu[i,,j])
err[i,j]<-  pow(pow(y2[i,j],.5) - pow(expected[i,j],.5),2)
errnew[i,j]<- pow(pow(ncapsnew2[i,j],.5) - pow(expected[i,j],.5),2)
}
ytotnew[i]<-sum(ncapsnew2[i,])
ytot.exp[i]<- sum(expected[i,])

err2[i]<- pow(pow(ytot[i],.5) - pow(ytot.exp[i],.5),2)
err2new[i]<- pow(pow(ytotnew[i],.5) - pow(ytot.exp[i],.5),2)
}

for(j in 1:ntraps){
traptotals.obs[j]<-sum(y2[,j])
traptotals.new[j]<-sum(ncapsnew2[,j])
traptotals.exp[j]<-sum(expected[,j])
err3[j]<- pow(pow(traptotals.obs[j],.5) - pow(traptotals.exp[j],.5),2)
err3new[j]<- pow(pow(traptotals.new[j],.5) - pow(traptotals.exp[j],.5),2)
}

X1obs<-sum(err[,])
X1new<-sum(errnew[,])
X2obs<-sum(err2[])
X2new<-sum(err2new[])
X3obs<-sum(err3[])
X3new<-sum(err3new[])


N<-sum(w[1:M])
D<-N/area
}
",fill=TRUE)
sink()

data <- list ("y","traplocs","M","ntraps","K","Xl","Xu","Yl","Yu","area","B","MASK")
sst<-sample(1:nrow(traplocs),M,replace=TRUE)
wst<-c(rep(1,nind),rep(0,M-nind))
inits <- function(){
  list (sigma=runif(1,.4,1),p0=runif(1,.01,.2),w=wst,alpha2=.25)
}
parameters <- c("psi","sigma","p0","N","D","beta","alpha2","X1obs",
"X1new","X2obs","X2new","X3obs","X3new")
out <- bugs(data, inits, parameters, "modelfile.txt", n.thin=1,n.chains=3, n.burnin=nb,
n.iter=ni,working.dir=getwd(),debug=TRUE)
out
}
