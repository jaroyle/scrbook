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
n.burnin=500,n.iter=10500,debug=TRUE,working.dir=getwd())

beta0<-out$sims.list$beta0
beta1<-out$sims.list$beta1

par(mfrow=c(2,1))
plot(beta0,type="l",lwd=2,ylab="parameter value",xlab="MCMC iteration")
plot(beta1,type="l",lwd=2,ylab="parameter value",xlab="MCMC iteration")