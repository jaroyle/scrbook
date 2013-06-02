since 2002   274 cites
since 2003   274 cites 0 articles

N<- c(3,2,5,3,8,11,20,46,65,84,27)
yr<- c(2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013)
delta<- c(12,12,12,12,12,12,12,12,12,12,2)
total<-cumsum(c(delta))
elapsed.time<-total

plot(yr,N)
lik<-function(parm){
n0<-N[1]
r<-exp(parm)
pred<- n0*exp(r*(elapsed.time[-1]-12)/12)

ss<- sum( (N[-1] - pred)^2)
return(ss)

}
out<-nlm(lik,c(log(.01)),hessian=TRUE)
r.hat<-exp(out$estimate)

yr.pred<- 2003:2020
delta.pred<- c(12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12)  ## through 2020
total<-cumsum(c(delta.pred))
elapsed.time<-total
N.pred<- N[1]*exp(r.hat*(elapsed.time -12)/12)
png("exp_growth.png",width=7,height=7, units="in", res=400)

plot(yr.pred,N.pred,xlab="Year",ylab="# of Published Articles",cex.axis=1.25,cex.lab=1.25)
points(yr[1:10],N[1:10],pch=20)

dev.off()

since 2004   271 cites 3 articles published in 2003
since 2005   269 cites 2 articles published in 2004  Efford 2004
since 2006   264 cites 5 articles
since 2007   261 cites 3 articles
since 2008   253 cites 8 articles Borchers and Efford and Royle 
since 2009   242 cites 11 articles
since 2010   222 cites 20 articles
since 2011   176 cites 46 articles
since 2012   111 cites 65 articles
since 2013    27 cites 84 articles published in 2012
                       27 so far since March 6
