\name{fakeshapefile}
\alias{fakeshapefile}
\docType{data}
\title{
this is the objec%%   ~~ data name/kind ... ~~
}
\description{
  An object of class "SpatialPolygonsDataFrame" from the package "sp"


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(fakeshapefile)}
\format{
see: str(SSp)
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library("scrbook")
data(fakeshapefile)
library(maptools)
library(sp)

data<-simSCR0.fn(discard0=TRUE,sd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
Xl<-data$xlim[1]
Yl<-data$ylim[1]
Xu<-data$xlim[2]
Yu<-data$ylim[2]

delta<-.3
ssbuffer<-2
Xl<-min(X[,1]) -ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])- ssbuffer
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))

data("fakeshapefile")
#### replaces this:
#####SSp<-readShapeSpatial('Sim_Polygon.shp')
Pcoord<-SpatialPoints(G)
PinPoly<-over(Pcoord,SSp)
Pin<-as.numeric(!is.na(PinPoly[,1]))
### over() returns NA when the point is not within any polygon, or the ID number of the polygon the point is in; 
### so it works for more complex multi-polygons, too. 
G<-G[Pin==1,]
frog<-nlm(intlik4,c(-2.5,2,log(4)),hessian=TRUE,y=y,K=K,delta=.3,X=traplocs,G=G)

## maybe str(fakeshapefile) ; plot(fakeshapefile) ...
}
\keyword{datasets}
