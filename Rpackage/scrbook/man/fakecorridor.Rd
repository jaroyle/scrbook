\name{fakecorridor}
\alias{fakecorridor}
\docType{data}
\title{
Fake buffered river corridor for ecological distance stuff
}
\description{
Fake buffered river corridor for ecological distance stuff
}
\usage{data(fakecorridor)}
\format{
R data format
}
\details{
  produces two files called "l1" and "l2" which are both of type
  "SpatialLines" from the package "sp"
  
  Also has file "traps" which are 5 trap locations

  Requires library "sp"  --
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library("sp")
library("rgeos")
library("scrbook")
data("fakecorridor")
buffer<- 0.5
par(mfrow=c(1,1))
aa<-gUnion(l1,l2)
plot(gBuffer(aa,width=buffer),xlim=c(0,10),ylim=c(0,10))
pg<-gBuffer(aa,width=buffer)
pg.coords<- pg@polygons[[1]]@Polygons[[1]]@coords

xg<-seq(0,10,,30)
yg<-seq(10,0,,30)

delta<-mean(diff(xg))
pts<- cbind(sort(rep(xg,30)),rep(yg,30))
points(pts,pch=20)
points(traps$locs,pch=20,cex=2,col="red")

in.pts<-point.in.polygon(pts[,1],pts[,2],pg.coords[,1],pg.coords[,2])
points(pts[in.pts==1,],pch=20,col="red")
}
\keyword{datasets}
