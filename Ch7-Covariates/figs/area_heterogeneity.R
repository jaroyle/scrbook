png("area_heterogeneity.png",width=7,height=3.5, units="in", res=400)

par(mfrow=c(1,2))
 tau<-.4
 mean( a95<- 6*pi*exp(rnorm(100000,-1-.5*(tau*tau),tau)))

sqrt(var(a95) )

hist(a95,nclass=100,xlim=c(0,30),xlab="95% HR area",ylab="Probability",probability=TRUE,main="SD = 2.88")


 tau<-.1
 mean( a95<- 6*pi*exp(rnorm(100000,-1-.5*(tau*tau),tau)))
sqrt(var(a95))
 

hist(a95,nclass=100,xlim=c(0,30),xlab="95% HR area",ylab="Probability",probability=TRUE,main="SD = 0.70")


dev.off()