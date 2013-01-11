figname<-"blahblah"

### not sure of dimensions to match the current version of the fig.

png(figname,width=7,height=7, units="in", res=400)

set.seed(9823)
traplocs<- cbind(sort(rep(1:5,5)),rep(1:5,5))
Dmat<-e2dist(traplocs,traplocs)
ntraps<-nrow(traplocs)

delta<-2
Xl<-min(traplocs[,1] - delta)
Xu<-max(traplocs[,1] + delta)
Yl<-min(traplocs[,2] - delta)
Yu<-max(traplocs[,2] + delta)
#plot(traplocs,pch=20,cex=1.5,xlim=c(Xl,Xu),ylim=c(Yl,Yu),xlab=" ",ylab=" ")
plot(traplocs,pch=20,cex=2.5,xlim=c(Xl,Xu),ylim=c(Yl,Yu),xlab=" ",ylab=" ",cex.axis=1.5)


sx<-runif(20,Xl,Xu)
sy<-runif(20,Yl,Yu)
#points(sx,sy,pch=20,col="red")
points(sx,sy,pch=20,col="red",cex=2)
csx<-cut(sx,breaks=c(Xl:Xu),include.lowest=TRUE)
csy<-cut(sy,breaks=c(Yl:Yu),include.lowest=TRUE)
table(csx,csy)





dev.off()