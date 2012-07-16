\name{beardata}
\alias{beardata}
\docType{data}
\title{
The bear data
}
\description{
The data on bears
}
\usage{data(beardata)}
\format{
  The format is:
 chr "beardata"
}
\details{
loads the Fort Drum beardata.
A list with "trapmat" and "bearArray" (3 d array)
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(beardata)
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=500
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum)
y[y>1]<-1
y<-apply(y,1,sum)   # total encounters out of K

set.seed(2013)

out<-modelMh(y,K,nsim=11000)

# MCMC CODE HERE

## maybe str(beardata) ; plot(beardata) ...
}
\keyword{datasets}
