SCR0bayes <-
function(data,M=200,engine="jags",ni=2000,nb=1000){
if(sum(engine==c("jags","winbugs"))==0) { return("use jags or winbugs!")  }

y<-data$Y
if(length(dim(y))!=2)
     return("Data must be 2-d array, nind x ntraps")

traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
K<-data$K
J<-nrow(X)
xlim<-data$xlim
ylim<-data$ylim

area<- (max(xlim)-min(xlim))*(max(ylim)-min(ylim))

## Data augmentation stuff
y<-rbind(y,matrix(0,nrow=M-nind,ncol=ncol(y)))
z<-c(rep(1,nind),rep(0,M-nind))

cat("
model {
alpha0~dnorm(0,.1)
logit(p0)<- alpha0
alpha1~dnorm(0,.1)
psi~dunif(0,1)
for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(xlim[1],xlim[2])
 s[i,2]~dunif(ylim[1],ylim[2])
 for(j in 1:J){
    d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
    y[i,j] ~ dbin(p[i,j],K)
    p[i,j]<- z[i]*p0*exp(- alpha1*d[i,j]*d[i,j])
 }
}
N<-sum(z[])
D<- N/area
}
",file = "SCR0b.txt")

sst<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))  # starting values for s
# Starting values for activity centers are taken as the mean encounter location
for(i in 1:nind){
if(sum(y[i,])==0) next
sst[i,1]<- mean( X[y[i,]>0,1] )
sst[i,2]<- mean( X[y[i,]>0,2] )
}
# Package up the data for analysis by BUGS or JAGS
data <- list (y=y,X=X,K=K,M=M,J=J,xlim=xlim,ylim=ylim,area=area)
inits <- function(){
  list (alpha0=rnorm(1,-4,.4),alpha1=runif(1,1,2),s=sst,z=z)
}
parameters <- c("alpha0","alpha1","N","D")

nthin<-1
nc<-3
if(engine=="winbugs"){
library("R2WinBUGS")
out <- bugs (data, inits, parameters, "SCR0b.txt", n.thin=nthin,n.chains=nc,
 n.burnin=nb,n.iter=ni,debug=FALSE,working.dir=getwd())
}
if(engine=="jags"){
library("rjags")
jm<- jags.model("SCR0b.txt", data=data, inits=inits, n.chains=nc,
                 n.adapt=nb)
out<- coda.samples(jm, parameters, n.iter=ni-nb, thin=nthin)
}

return(out)
}
