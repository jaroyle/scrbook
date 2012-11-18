figname<-"det_functions.png"

### not sure of dimensions to match the current version of the fig.
png(figname,width=7,height=4.5, units="in", res=400)

x<- seq(.011,20,,100)
.

k<- exp(-(1/(2*6*6))*x^2)
p1<- .65*k
plot(x,p1,type="l",lwd=2,xlab="distance",ylab="encounter probability")
p2<- 1-exp(-p1)
lines(x,p2,lty=2,lwd=2)


k<- exp(-(1/(2*1*1))*x)
p3<- .65*k
lines(x,p3,lty=3,lwd=2)


lp4<-  0 - (1/(2*6*6))*x*x
p4<-exp(lp4)/(1+exp(lp4))
lines(x,p4,lty=4,lwd=2)

legend(12,.6,legend=c("Gaussian detection","Gaussian hazard","negative exponential","logistic"),lty=1:4,cex=.8)

dev.off()