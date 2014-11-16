\name{intlik3edv2}
\alias{intlik3edv2}
\title{
 Integrated likelihood of the binomial SCR model
}
\description{
  Computes the integrated likelihood for binomial SCR model with
  a prescribed distance matrix, which may be normal Euclidean distance
  or it may be ecological distance. 
}
\usage{
intlik3edv2(start,y,K,X,S,D,inpoly)
}
\arguments{
  \item{start}{
    values of parameters to evaluate likelihood with.
    order is intercept (alpha0) distance coefficent (alpha1) and log(n0)
}
  \item{y}{
    Data....
}
\item{K}{
  number reps
}
\item{X}{
  trap coordinates
}
\item{S}{
 state-space grid
  }
\item{D}{
  distance matrix between traps and every state-space point. e.g., 
  D = e2dist(traplocs,S)
  }
\item{inpoly}{
  binary indicator if state-space point should be used in the marginal
  likelihood calculation
  }
}
\value{
returns -1*loglikelihood
}
\author{
J. Andy Royle
}
\examples{

library("sp")
library("rgeos")
library("scrbook")
library("raster")
data("fakecorridor")

## Step 1
## produce the corridor system:
buffer<- 0.5
par(mfrow=c(1,1))
aa<-gUnion(l1,l2)
plot(gBuffer(aa,width=buffer),xlim=c(0,10),ylim=c(0,10))
pg<-gBuffer(aa,width=buffer)
pg.coords<- pg@polygons[[1]]@Polygons[[1]]@coords
 # note: can you believe this shit?

xg<-seq(0,10,,40)
yg<-seq(10,0,,40)
delta<-mean(diff(xg))
pts<- cbind(sort(rep(xg,40)),rep(yg,40))
points(pts,pch=20,cex=.5)
in.pts<-point.in.polygon(pts[,1],pts[,2],pg.coords[,1],pg.coords[,2])
points(pts[in.pts==1,],pch=20,col="red")

### Step 2
### Assign cost matrix
cost<-rep(NA,nrow(pts))
cost[in.pts==1]<-1   # low cost to move among pixels but not 0
cost[in.pts!=1]<-10000   # high cost 

## Stuff this into a raster
library("raster")
r<-raster(nrows=40,ncols=40)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(0-delta/2,10+delta/2,0-delta/2,10+delta/2)
values(r)<-matrix(cost,40,40,byrow=FALSE)
par(mfrow=c(1,1))
plot(r)
points(pts,pch=20,cex=.4)


library("gdistance")
## use max = doesn't count moving through boundary pixel
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
costs1<-costDistance(tr1CorrC,pts)
outD<-as.matrix(costs1)

plot(pts,pch=".")
points(pts[in.pts==1,],pch=20,col="red")

#
# Note "traps" file in scrbook library based on 40 x 40 raster
#
traplocs<-traps$loc
trap.id<-traps$locid
ntraps<-nrow(traplocs)


set.seed(2013)
###
### Step 3: Simulate some SCR data
###
N<-200
# restrict to points in the corridor system
S.possible<- (1:nrow(pts))[in.pts==1]
S.id<-sample(S.possible,N,replace=TRUE)
S<- pts[S.id,]

D<- outD[S.id,trap.id]
eD<- e2dist(S,traplocs)
# the follow two distance matrices will be used in the likelihood
Dtraps<-outD[trap.id,]
Deuclid<-e2dist(pts[trap.id,],pts)

alpha0<- -1.5
sigma<- 1.5
alpha1<- 1/(2*sigma*sigma)
K<-10

probcap<-plogis(alpha0)*exp(-alpha1*D*D)
probcapE<-plogis(alpha0)*exp(-alpha1*eD*eD)

Y<-matrix(NA,nrow=N,ncol=ntraps)
Ye<-Y
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
 Ye[i,]<-rbinom(ntraps,K,probcapE[i,])
}

Y<-Y[apply(Y,1,sum)>0,]

###
###  Step 4: Run SCR model
###
# first with distance-in-corridor
frog1<-nlm(intlik3edv2,c(-2.5,2,log(4)),hessian=TRUE,y=Y,K=K,X=traplocs,S=pts,D=Dtraps,inpoly=in.pts)
# now run Euclidean distance
frog2<-nlm(intlik3edv2,c(-2.5,2,log(4)),hessian=TRUE,y=Y,K=K,X=traplocs,S=pts,D=Deuclid,inpoly=in.pts)


}
