\name{make.seg}
\alias{make.seg}
\title{
Make a line segment by laying down a few points.
}
\description{
Make a line seg
}
\usage{
make.seg(npts=10)
}
\arguments{
  \item{npts}{
    how many points should your line have?
}
}
\value{
an object of class "SpatialLines" (R package "sp" I think).
}
\author{
J. Andy Royle
}
\examples{
# make two line segments however you wish them to look
library("scrbook")
plot(NULL,xlim=c(0,10),ylim=c(0,10))
l1<-make.seg(9)
plot(l1)
l2<-make.seg(5)
plot(l1)
lines(l2)
}
