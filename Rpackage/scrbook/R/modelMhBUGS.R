modelMhBUGS <-
function(Yarr,engine="winbugs",M=500,ni=10000,nb=1000){
if( sum(engine == c("jags","winbugs"))==0) return("not a valid MCMC engine: use winbugs or jags")

nind<-dim(Yarr)[1]
K<-dim(Yarr)[3]
ntraps<- dim(Yarr)[2]

nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum) # summarize by ind x rep
y[y>1]<- 1             # toss out duplicate obs
ytot<-apply(y,1,sum)   # total encounters out of K


cat("
model{
p0 ~ dunif(0,1)       # prior distributions
mup<- log(p0/(1-p0))
sigmap ~ dunif(0,10)
taup<- 1/(sigmap*sigmap)
psi~dunif(0,1)

for(i in 1:(nind+nz)){
  z[i]~dbern(psi)     # zero inflation variables
  lp[i] ~ dnorm(mup,taup) # individual effect
  logit(p[i])<-lp[i]
  mu[i]<-z[i]*p[i]
  y[i]~dbin(mu[i],K)  #  observation model
 }

N<-sum(z[1:(nind+nz)])
}
",file="modelMh.txt")


data1<-list(y=ytot, nz=nz, nind=nind,K=K) 
params1= c('p0','sigmap','psi','N')
inits =  function() {list(z=as.numeric(ytot>=1), psi=.6, p0=runif(1),sigmap=runif(1,.7,1.2),lp=rnorm(M,-2)) }

if(engine=="winbugs"){
library("R2WinBUGS")
out = bugs(data1, inits, params1, model.file="modelMh.txt",working.directory=getwd(),    
       debug=FALSE, n.chains=3, n.iter=(ni+nb), n.burnin=nb, n.thin=1)
}


##summary(as.mcmc.list(fit1))
if(engine=="jags"){
library("rjags")
jm<- jags.model("modelMh.txt", data=data1, inits=inits, n.chains=3,
                 n.adapt=ni)
out<- coda.samples(jm, params1, n.iter=ni, thin=1)

}


out
}
