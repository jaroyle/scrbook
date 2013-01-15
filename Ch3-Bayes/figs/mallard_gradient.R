
library(scrbook)
data(mallard)

mallard<-list(bandings=mallard$bandings,recoveries=mallard$recoveries,locs=mallard$locs)


sink("model.txt")
cat("
model {
 for(t in 1:5){
    for (i in 1:nobs){
       y[i,t] ~ dbin(p[i,t], B[i,t])
       logit(p[i,t]) <- beta0[t] + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,1]*X[i,2]
     }
}
	beta1~dnorm(0,.001)
	beta2~dnorm(0,.001)
	beta3~dnorm(0,.001)
	for(t in 1:5){
 	beta0[t] ~ dnorm(0,.001)  
 }
}
",fill=TRUE)
sink()

library("R2WinBUGS")

data <- list(B=mallard$bandings, y=mallard$recoveries,
X=mallard$locs,nobs=nrow(mallard$locs))
inits <- 	function(){
list(beta0=rnorm(5),beta1=0,beta2=0,beta3=0)
}
parms <- list('beta0','beta1','beta2','beta3')
out <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=TRUE)


### now run to get estimates of p at each grid cell
### and make a nice image plot

parms <- list('beta0','beta1','beta2','beta3','p')
out2 <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=FALSE)

pbar<-apply(out2$sims.list$p,2,mean)
ux<-unique(mallard$locs[,1])
uy<-unique(mallard$locs[,2])
X<-mallard$locs
par(mfrow=c(2,1))
odx<-order(X[,1],X[,2])
I<-matrix(pbar[odx],nrow=length(uy),ncol=length(ux))

png("mallard_gradient.png",width=5,height=3.5, units="in", res=400)
par(mfrow=c(1,1),mar=c(4,4,4,6))
plot(X,pch=" ",ylab="Northing",xlab="Easting")
image(ux,uy,t(I),col=terrain.colors(10),add=TRUE)
image.scale(I,col=terrain.colors(10))
dev.off()