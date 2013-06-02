Tabitha code:


library(raster)
r<-raster(nrows=4,ncols=4)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84" #sets the projection

extent(r)<-c(.5,4.5,.5,4.5) #sets the extent of the raster
costs1<- c(100,100,100,100,1,100,100,100,1,1,100,1,1,1,1,1)
values(r)<-matrix(costs1,4,4,byrow=FALSE) #assign the costs to the raster

library(gdistance)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)

tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
pts<-cbind( sort(rep(1:4,4)),rep(1:4,4))
#the costDistance function calculates the least cost distane among pts using
#the conductances specified in tri1CorrC
costs1<-costDistance(tr1CorrC,pts)
#here we convert costs1 into a matrix
outD<-as.matrix(costs1)
outD[1:5,1:5]

Andy code

r<-raster(nrows=4,ncols=4)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
costs1<- c(100,100,100,100,1,100,100,100,1,1,100,1,1,1,1,1)
values(r)<-matrix(costs1,4,4,byrow=FALSE)

values(r)<-matrix(costs1,4,4,byrow=FALSE)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
pts<-cbind( sort(rep(1:4,4)),rep(4:1,4))
costs1<-costDistance(tr1CorrC,pts)
outD<-as.matrix(costs1)
outD[1:5,1:5]