\name{intlik3ed}
\alias{intlik3ed}
\title{
Computes the integrated likelihood for binomial SCR model with
ecological or Euclidean distance.
}
\description{
  Computes integrated likelihood for binomial SCR model with ecological
  or Euclidean distance metric.  The model used here is the Bernoulli
  encounter model so that y[i,j,k] = 1 if individual i is encountered in
  trap j during sample occasion k. The function uses the summary
  statstioncs y2d[i,j] = number of encounters out of K. 
}
\usage{
intlik3ed(start,y,K,X,distmet="ecol",covariate,alpha2=NA)
}
\arguments{
  \item{start}{
    values of parameters at which to evaluate likelihood.
    order is alpha0, alpha1, n0, alpha2. If distmet="euclid" then
    dont provide a starting value for alpha2.
}
  \item{y}{
   nind x ntraps matrix of encounter frequencies from a study based on K
   sample occasions. 
}
\item{K}{
  Number of sampling occasions.
}
\item{X}{
  Trap coordinates
}
\item{distmet}{
Either "ecol" or "euclid". "ecol" is the least-cost path model. 
  }
  \item{covariate}{
    The covariate to be used in the cost function.
    cost(x) = exp(alpha2*covariate(x))
  }
  \item{alpha2}{
    You can fix the value of alpha2 using this argument. In that case
    start should be length 3 and alpha2 is not estimated.
    if dist="euclid" then you should have alpha2=NA
  }
}
\details{
  more here
}
\value{
Returns -1*loglikelihood
}
\references{
  SCR book, Chapter 12. 

  Royle, J. A., Chandler, R. B., Gazenski, K. D., & Graves,
  T. A. (2013). Spatial capture-recapture models for jointly estimating
  population density and landscape connectivity. Ecology, 94(2),
  287-294.

  Sutherland, C., A.K. Fuller and J.A. Royle. Modelling non-Euclidean
  movement and landscape connectivity in highly-structured ecological
  networks. (in review)
  
  }
\author{
J. Andy Royle, aroyle@usgs.gov
}
\examples{

library("scrbook")
library("gdistance")
library("sp")
library("rgeos")
library("raster")

out<-make.EDcovariates()
covariate<-out$covariate.patchy
set.seed(2013)

N<-200
alpha0<- -2
sigma<- .5
K<- 5

alpha1<- 1/(2*sigma*sigma)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
alpha2<-1
cost<- exp(alpha2*covariate)

tr1<-transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

# make up some trap locations
xg<-seq(1,4,1)
yg<-4:1
pts<-cbind( sort(rep(xg,4)),rep(yg,4))

traplocs<-pts
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))
D<-costDistance(tr1CorrC,S,traplocs)
probcap<-plogis(alpha0)*exp(-alpha1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}
Y<-Y[apply(Y,1,sum)>0,]

n0<- N-nrow(Y)
frog1<-nlm(intlik3ed,c(alpha0,alpha1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="euclid",covariate=covariate,alpha2=1)

frog2<-nlm(intlik3ed,c(alpha0,alpha1,log(n0),-.3),hessian=TRUE,y=Y,K=K,
               X=traplocs,distmet="ecol",covariate=covariate,alpha2=NA)

}
