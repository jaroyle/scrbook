\name{modelMh}
\alias{modelMh}
\title{
fits logit-normal model Mh
}
\description{
......more here
}
\usage{
modelMh(y,K,nsim)
}
\arguments{
  \item{y}{
 vector of detections (augmented)
}
  \item{K}{
   number of samples of a closed population
}
\item{nsim}{
  how many mcmc iterations?
  }
}
\value{
 .... a matrix ...
}
\author{
Andy Royle and Rahel Sollmann
}
\examples{
library("scrbook")
data("beardata")
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=500
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum)
y[y>1]<-1
ytot<-apply(y,1,sum)   # total encounters out of K

set.seed(2013)
out<-modelMh(ytot,K,nsim=11000)
}
