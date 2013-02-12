\name{get.traplocs}
\alias{get.traplocs}
\title{
Make a trapping grid by clicking on some points.
}
\description{
Make a trapping grid.
}
\usage{
get.traplocs(ntraps,ssgrid)
}
\arguments{
  \item{ntraps}{    how many traps?}
  \item{ssgrid}{ the state-space grid -- a bunch of points. defaults
  to NULL so that if a plot window is already open you can add points to
   that. }
}
\value{
a list with locs and locid being the integer id of the pixel 
}
\author{
J. Andy Royle
}
\examples{
# make two line segments however you wish them to look
library("scrbook")

xg<-seq(0,10,,40)
yg<-seq(10,0,,40)
pts<- cbind(sort(rep(xg,40)),rep(yg,40))

out<- get.traplocs(ntraps=10,ssgrid=pts)

points(out$loc,pch=20,col="red")


}
