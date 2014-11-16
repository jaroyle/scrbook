make.EDcovariates<-function(){

### We require 3 R libraries
library("gdistance")
library("raster")

###
### following block of code creates a "patchy" looking covariate to use
### in our cost function. It uses a standard method for generating a correlated
### multivariate normal vector of length, in this case, 400 (one value for each
### pixel). One can use any correlation function here but we chose a standard
### exponential model with range parameter 0.5
###png("raster_krige.png",width=5,height=5, units="in", res=400)
par(mfrow=c(1,1))
set.seed(12)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
delta<- (4.5-.5)/20
gx<- seq(.5 + delta/2, 4.5-delta/2,,20)
gy<- rev(gx)
gx<-sort(rep(gx,20))
gy<-rep(gy,20)
grid<-cbind(gx,gy)
Dmat<-as.matrix(dist(grid))
V<-exp(-Dmat/.5)
z<-t(chol(V))%*%rnorm(400)
z<- (z-mean(z))/sqrt(var(as.vector(z)))
#values(r)<-rot(rot(rot(matrix(z,20,20,byrow=TRUE))))
#values(r)<-(matrix(z,20,20,byrow=FALSE))
values(r)<-matrix(z,20,20,byrow=FALSE)
#spatial.plot(grid,z)
par(mfrow=c(2,1))
plot(r)
####dev.off()
covariate.patchy<- z*sqrt(1.68) + .168   #approx. same scale as systematic covariate below
values(r)<-matrix(covariate.patchy,20,20,byrow=FALSE)
covariate.patchy<- r
class(covariate.patchy)

##
## build systematic covariate
## defined here as a trend from NW to SE
##
cost<-matrix(NA,nrow=20,ncol=20)
cost<-row(cost)+col(cost)
covariate.trend<- (cost-20)/10
values(r)<-matrix(covariate.trend,20,20,byrow=FALSE)
covariate.trend<-r
class(covariate.trend)
plot(r)
list(covariate.patchy=covariate.patchy,covariate.trend=covariate.trend)
}
