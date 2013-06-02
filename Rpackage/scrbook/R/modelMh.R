modelMh <-
function(ytot,K,nsim=1000){

out<-matrix(NA,nrow=nsim,ncol=4)
dimnames(out)<-list(NULL,c("mu","sigma","psi","N"))
lp<- rnorm(M,-1,1)
p<-plogis(lp)
mu<- -1
p0<-exp(mu)/(1+exp(mu))
sigma<- 1
psi<- .5
z<-rbinom(M,1,psi)
z[ytot>0]<-1

for(i in 1:nsim){

### update the logit(p) parameters
lp.cand<- rnorm(M,lp,1)  # 1 is a tuning parameter
p.cand<-plogis(lp.cand)
ll<-dbinom(ytot,K,z*p, log=T)
prior<-dnorm(lp,mu,sigma, log=T)
llcand<-dbinom(ytot,K,z*p.cand, log=T)
prior.cand<-dnorm(lp.cand,mu,sigma, log=T)

kp<- runif(M) < exp((llcand+prior.cand)-(ll+prior))
p[kp]<-p.cand[kp]
lp[kp]<-lp.cand[kp]

p0.cand<- rnorm(1,p0,.05)
if(p0.cand>0 & p0.cand<1){
mu.cand<-log(p0.cand/(1-p0.cand))
ll<-sum(dnorm(lp,mu,sigma,log=TRUE))
llcand<-sum(dnorm(lp,mu.cand,sigma,log=TRUE))
if(runif(1)<exp(llcand-ll)) {
 mu<-mu.cand
 p0<-p0.cand
}
}

sigma.cand<-rnorm(1,sigma,.5)
if(sigma.cand>0){
ll<-sum(dnorm(lp,mu,sigma,log=TRUE))
llcand<-sum(dnorm(lp,mu,sigma.cand,log=TRUE))
if(runif(1)<exp(llcand-ll))
 sigma<-sigma.cand
}


### update the z[i] variables
z.cand<-  ifelse(z==1,0,1)  # candidate is 0 if current = 1, etc..
ll<- dbinom(ytot,K,z*p, log=TRUE)
prior<-dbinom(z,1,psi, log=TRUE)
llcand<- dbinom(ytot,K,z.cand*p, log=TRUE)
prior.cand<-dbinom(z.cand,1,psi, log=TRUE)
kp<- runif(M) <  exp((llcand+prior.cand)-(ll+prior))
z[kp]<- z.cand[kp]

psi<-rbeta(1, sum(z) + 1, M-sum(z) + 1)

out[i,]<- c(mu,sigma,psi,sum(z))

}

return(out)

}