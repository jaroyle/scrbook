library("scrbook")

data(bbsdata)
y<-bbsdata$counts[,"X90"] # pick out 1990
notna<-!is.na(y)
y<-y[notna]
## forest cover already standardized here:
habitat<-bbsdata$counts[notna,"habitat"]
M<-length(y)
library("R2WinBUGS") # load R2WinBUGS
data <- list ( "y","M","habitat") # bundle data for WinBUGS


cat("
model {
for (i in 1:M){
y[i]~dpois(lam[i])
log(lam[i])<- beta0+beta1*habitat[i]
}
beta0~dunif(-5,5)
beta1~dunif(-5,5)
}
",file="PoissonGLM.txt")
inits <- function() list ( beta0=rnorm(1),beta1=rnorm(1))
parameters <- c("beta0","beta1")
out<-bugs (data, inits, parameters, "PoissonGLM.txt", n.thin=2,n.chains=2,
n.burnin=2000,n.iter=6000,debug=TRUE,working.dir=getwd())


set.seed(2013) # so we all get the same result
out<-matrix(NA,nrow=1000,ncol=2) # matrix to store the output
beta0<- -1 # starting values
beta1 <- -.8
# begin the MCMC loop ; do 1000 iterations
for(i in 1:1000){
# update the beta0 parameter
lambda<- exp(beta0+beta1*habitat)
lik.curr<- sum(log(dpois(y,lambda)))
prior.curr<- log(dnorm(beta0,0,100))
beta0.cand<-rnorm(1,beta0,.05) # generate candidate
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
out[i,]<-c(beta0,beta1) # save the current values
}


plot(out[1:300,1],ylim=c(-1.5,3.3),type="l",lwd=2,ylab="parameter value",
xlab="MCMC iteration")
lines(out[,2],lwd=2) 