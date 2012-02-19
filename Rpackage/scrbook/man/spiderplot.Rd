\name{spiderplot}
\alias{spiderplot}
\title{
Makes a spider plot of some SCR data.
}
\description{
Makes a spider plot and computes some summaries of individual captures.
}
\usage{
spiderplot(y3d,traplocs)
}
\arguments{
  \item{y3d}{
 3-d encounter history array: nind x ntraps x nreps
}
  \item{traplos}{
   ntraps x 2 matrix of trap coordinates.
}
}
\value{
 Makes a spider plot. Retunrs list with 

  xcent = distance from trap array centroid of each guy.

  avg.s = average capture location

  center = centroid of trap array
}
\author{
J. Andy Royle
}
\examples{
data("scrbook")
data(beardata)
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]
M=500
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-beardata$bearArray
toad<- spiderplot(beardata$bearArray,beardata$trapmat)
# now grab the distance from centroid variable
xcent<-toad$xcent
# see Chapter 3 of the book
}
