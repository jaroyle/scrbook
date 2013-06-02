data(scrbook)
data(tortoise)

library("unmarked")
data.mn<- formatDistData(tortoise,"Dist","Transect",dist.breaks=c(0:32))
dimnames(data.mn)<-list(NULL,c(1:32))

png("tortoise.png",width=7,height=7, units="in", res=400)
 plot(apply(data.mn,2,sum))
apply(data.mn,2,sum)->Nv

 barplot(apply(data.mn,2,sum),xlab="distance (m)",ylab="frequency (no. tortoises)",
cex.axis=1.5,cex.lab=1.5)

dev.off()
