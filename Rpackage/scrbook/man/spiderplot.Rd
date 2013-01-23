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

library("scrbook")

# make spiderplot for the beardata the hard way:
data(beardata)
X<-as.matrix(cbind(id=1:nrow(beardata$trapmat),beardata$trapmat))
opps<-matrix(1,nrow=nrow(X),ncol=8)
dimnames(opps)<-list(NULL,1:8)
X<-cbind(X,opps)
a<-SCR23darray(beardata$flat,X)
toad<-spiderplot(a,beardata$trapmat)
# now grab the distance from centroid variable
xcent<-toad$xcent
# see Chapter 3 of the book

# the easy way
spiderplot(beardata$bearArray,beardata$trapmat)

# for the wolverine data:
y3d <- SCR23darray(wolverine$wcaps,wolverine$wtraps)
spiderplot(y3d,wolverine$wtraps[,2:3])





}
