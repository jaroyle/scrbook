PoisGLMBBS<-function(y, habitat, niter) {
out<-matrix(NA,nrow=niter,ncol=2)   # matrix to store the output
beta0<- -1                         # starting values
beta1 <- -.8

# begin the MCMC loop ; do niter iterations
for(i in 1:niter){

# update the beta0 parameter
lambda<- exp(beta0+beta1*habitat)
lik.curr<- sum(log(dpois(y,lambda)))
prior.curr<- log(dnorm(beta0,0,100))
beta0.cand<-rnorm(1,beta0,.05)         # generate candidate
lambda.cand<- exp(beta0.cand + beta1*habitat)
lik.cand<- sum(log(dpois(y,lambda.cand)))
prior.cand<- log(dnorm(beta0.cand,0,100))
mhratio<- exp(lik.cand +prior.cand - lik.curr-prior.curr)
if(runif(1)< mhratio)
     beta0<-beta0.cand

# update the beta1 parameter
lik.curr<- sum(log(dpois(y,exp(beta0+beta1*habitat))))
prior.curr<- log(dnorm(beta1,0,100))
beta1.cand<-rnorm(1,beta1,.25)
lambda.cand<- exp(beta0+beta1.cand*habitat)
lik.cand<- sum(log(dpois(y,lambda.cand)))
prior.cand<- log(dnorm(beta1.cand,0,100))
mhratio<- exp(lik.cand + prior.cand - lik.curr - prior.curr)
if(runif(1)< mhratio)
     beta1<-beta1.cand

out[i,]<-c(beta0,beta1)             # save the current values
}
return(out)
}