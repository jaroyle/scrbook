library("scrbook")
data(bbsdata)
library("maps")

y<-bbsdata$counts[,"X90"] # pick out 1990
notna<-!is.na(y)
y<-y[notna]
locs<-bbsdata$counts[notna,c("lon","lat")]


sz<- y/max(y)
par(mar=c(3,3,3,6))
plot(locs,pch=" ",xlim=range(locs[,1])+c(-.3,+.3),ylim=c(range(locs[,2]) + c(-.6,.6)),
axes=FALSE,xlab=" ",ylab=" ")
map('state',regions='pennsylvania',add=TRUE,lwd=2)
spatial.plot(bbsdata$counts[notna,2:3],y,cx=1+sz*6,add=TRUE)