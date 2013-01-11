

spatial.plot<-
function (x, y, add = FALSE, cx = 1,col="gray") 
{
    nc <- as.numeric(cut(y, 10))
    if (!add) 
        plot(x, pch = " ")
        if(col=="gray"){
        cc<-seq(3,17,,10)/20
        cc<-gray(cc)
}
else 
 cc<-terrain.colors(10)
  points(x, pch = 20, col = cc[nc], cex = cx)
  image.scale(y, col = cc)
}
 

library("scrbook")
data(bbsdata)
library("maps")

y<-bbsdata$counts[,"X90"]  # pick out 1990
notna<-!is.na(y)
y<-y[notna]
locs<-bbsdata$counts[notna,c("lon","lat")]
sz<- y/max(y)
png("pacounts.png",width=5.5,height=4, units="in", res=400)
par(mar=c(3,3,3,6))
plot(locs,pch=" ",xlim=range(locs[,1])+c(-.3,+.3),ylim=c(range(locs[,2]) + c(-.6,.6)),axes=FALSE,xlab=" ",ylab=" ")
map('state',regions='pennsylvania',add=TRUE,lwd=2)
spatial.plot(bbsdata$counts[notna,2:3],y,cx=1+sz*6,add=TRUE)
dev.off()


library(maps)
habdata<-bbsdata$habitat
map('state',regions="penn",lwd=2)
I<-matrix(NA,nrow=30,ncol=40)
I<- matrix(habdata[,"dfor"],ncol=40,byrow=FALSE)
ux<-unique(habdata[,2])
uy<-sort(unique(habdata[,3]))
png("paforest.png",width=5.5,height=4, units="in", res=400)
par(mar=c(3,3,3,6))
plot(locs,pch=" ",xlim=range(locs[,1])+c(-.3,+.3),ylim=c(range(locs[,2]) + c(-.6,.6)),axes=FALSE,xlab=" ",ylab=" ")

image(ux,uy,rot(I),add=TRUE,col=gray(seq(3,17,,10)/20) )
map('state',regions='pennsylvania',add=TRUE,lwd=3,col="white")
image.scale(I,col=gray(seq(3,17,,10)/20) )
points(locs,pch=20,col="white")
dev.off()

#spatial.plot(habdata[,2:3],habdata[,"dfor"],cx=3)
#map('state',regions="penn",lwd=2,add=TRUE)
