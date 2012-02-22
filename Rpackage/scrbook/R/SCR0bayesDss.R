SCR0bayesDss <-
function(data,ng=10,engine="jags",ni=2000,nb=1000,M=200){
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

dd<- Yu-Yl
eps<- dd/(ng)
gridx<-seq(Xl+eps/2,Xu-eps/2,eps)
gridy<-seq(Yl+eps/2,Yu-eps/2,eps)
S<-cbind(sort(rep(gridx,ng)),rep(gridy,ng))
nG<-nrow(S)

y<-rbind(y,matrix(0,nrow=(M-nrow(y)),ncol=J ) )
sst<-S[sample(1:nrow(S),M,replace=TRUE),]

cat("
model {
alpha~dnorm(0,.1)
beta~dnorm(0,.1)
psi~dunif(0,1)

for(g in 1:nG){
probs[g]<- 1/nG
}

for(i in 1:M){
 z[i]~dbern(psi)
s[i] ~ dcat(probs[1:nG])
for(j in 1:J){
y[i,j] ~ dbin(pnew[i,j],K)
  d2[i,j]<- pow(S[s[i],1]-X[j,1],2) + pow(S[s[i],2]-X[j,2],2)
p[i,j]<- exp(alpha)*exp(-beta*d2[i,j])
pnew[i,j]<-p[i,j]*z[i]
}
}
N<-sum(z[])
}
",file = "SCR0-discrete.txt")

sst<-sample(1:nG,M,replace=TRUE)
zst<-c(rep(1,nind),rep(0,M-nind))
data <- list (y=y,X=X,K=K,M=M,J=J,S=S,nG=nG)
inits <- function(){
  list (alpha=rnorm(1,-2.5,.4),beta=rnorm(1,2,.5),psi=runif(1),z=zst ,s=sst)
}

parameters <- c("alpha","beta","psi","N")
nthin<-1
nc<-3
nb<-nb
ni<-ni
# 6x6 = 1110 seconds    762 secs on work computer
# 8x8 = 1623 seconds    1195 secs on work computer
# 10x10                 1753 seconds
# 12 x 12 = 4390        2654 secs
# 20 x 20 9214 seconds
# 25 x 25 10000 secs on work machine

if(engine=="winbugs"){
library("R2WinBUGS")
tm<-unix.time(
out <- bugs (data, inits, parameters, "SCR0-discrete.txt", n.thin=nthin,n.chains=nc,
 n.burnin=nb,n.iter=ni,debug=FALSE,working.dir=getwd())
)
}
if(engine=="jags"){
library("rjags")
jm<- jags.model("SCR0-discrete.txt", data=data, inits=inits, n.chains=nc,
                 n.adapt=nb)
tm<-unix.time(jm<- coda.samples(jm, parameters, n.iter=ni-nb, thin=nthin))
out<-jm

}

return(list(out=out,tm=tm))
}
