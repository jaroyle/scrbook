SCR0bayes <-
function(dataobj,M=200,engine="jags",ni=2000,nb=1000){
if(sum(engine==c("jags","winbugs"))==0) { return("use jags or winbugs!")  }

y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
K<-data$K
J<-nrow(X)
Xl<-data$xlim[1]
Yl<-data$ylim[1]
Xu<-data$xlim[2]
Yu<-data$ylim[2]

## Data augmentation stuff
y<-rbind(y,matrix(0,nrow=M-nind,ncol=ncol(y)))
z<-c(rep(1,nind),rep(0,M-nind))

cat("
model {
alpha0~dnorm(0,.1)
logit(p0)<- alpha0
beta~dnorm(0,.1)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu) 
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dbin(p[i,j],K)
p[i,j]<- z[i]*p0*exp(- beta*d[i,j]*d[i,j])
}
}
N<-sum(z[])
D<- N/64
}
",file = "SCR0a.txt")

sst<-cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))  # starting values for s
for(i in 1:nind){
if(sum(y[i,])==0) next
sst[i,1]<- mean( X[y[i,]>0,1] )
sst[i,2]<- mean( X[y[i,]>0,2] )
}
data <- list (y=y,X=X,K=K,M=M,J=J,Xl=Xl,Yl=Yl,Xu=Xu,Yu=Yu)
inits <- function(){
  list (alpha0=rnorm(1,-4,.4),beta=runif(1,1,2),s=sst,z=z)
}
parameters <- c("alpha0","beta","N","D")

nthin<-1
nc<-3
if(engine=="winbugs"){
library("R2WinBUGS")
out <- bugs (data, inits, parameters, "SCR0a.txt", n.thin=nthin,n.chains=nc,
 n.burnin=nb,n.iter=ni,debug=FALSE,working.dir=getwd())
}
if(engine=="jags"){
library("rjags")
jm<- jags.model("SCR0a.txt", data=data, inits=inits, n.chains=nc,
                 n.adapt=nb)
out<- coda.samples(jm, parameters, n.iter=ni-nb, thin=nthin)
}

return(out)
}
