library("scrbook")
data(bbsdata)
library("maps")

y<-bbsdata$counts[,"X90"] # pick out 1990
notna<-!is.na(y)
y<-y[notna]
locs<-bbsdata$counts[notna,c("lon","lat")]




habdata<-bbsdata$habitat
map('state',regions="penn",lwd=2)
I<-matrix(NA,nrow=30,ncol=40)
I<- matrix(habdata[,"dfor"],ncol=40,byrow=FALSE)
ux<-unique(habdata[,2])
uy<-sort(unique(habdata[,3]))

par(mar=c(3,3,3,6))
plot(locs,pch=" ",xlim=range(locs[,1])+c(-.3,+.3),ylim=c(range(locs[,2]) + c(-.6,.6)),axes=FALSE,xlab=" ",ylab=" ")
image(ux,uy,rot(I),add=TRUE,col=gray(seq(3,17,,10)/20) )
map('state',regions='pennsylvania',add=TRUE,lwd=3,col="white")
image.scale(I,col=gray(seq(3,17,,10)/20) )
points(locs,pch=20,col="white")