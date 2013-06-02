library(raster)
r<-raster(nrows=4,ncols=4)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84" #sets the projection
#We use UTM here because distances and directions are correct
#i.e., they are adjusted for the earth's curvature
extent(r)<-c(.5,4.5,.5,4.5) #sets the extent of the raster
costs1<- c(100,100,100,100,1,100,100,100,1,1,100,1,1,1,1,1)
values(r)<-matrix(costs1,4,4,byrow=FALSE) #assign the costs to the raster

#The geoCorrection function corrects the conductances for the diagonal neighbors
#and, if the raster is in latitude longitude, also corrects for 
#curvature of the earth's surface

tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

#here we specify the locations that we'd like to calculate distances among
# 
pts<-cbind( sort(rep(1:4,4)),rep(1:4,4))

#the costDistance function calculates the least cost distane among pts using
#the conductances specified in tri1CorrC
costs1<-costDistance(tr1CorrC,pts)
#here we convert costs1 into a matrix
outD<-as.matrix(costs1)

plot(r)
points(pts[1,1],pts[1,2],col="red")
points(pts[2,1],pts[2,2],col="blue")
lines(c(1,1),c(1,2))



points(pts[4,1],pts[4,2],col="blue")
lines(c(1,1),c(1,4))
#but the least cost distance is to go around the outside of the high cost area.
lines(c(1,2),c(1,1)) #with cost = (100 + 100)/2 = 100
lines(c(2,3),c(1,1)) #with cost = (100 + 1)/2 = 50.5
#To account for the increased distance along the diagonal, we multiply by sqrt(2)
lines(c(3,4),c(1,2)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
lines(c(4,3),c(2,3)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
lines(c(3,2),c(3,4)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
lines(c(2,1),c(4,4)) #with cost = (100 + 1)/2 = 50.5
100+50.5* 2+(1.414214)* 3 # = 205.2426 

