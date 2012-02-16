\name{intlik3ed}
\alias{intlik3ed}
\title{
Computes the integrated likelihood for binomial SCR model with
ecological
distance .
}
\description{
  Computes integrated likelihood for binomial SCR model with ecological
  distance metric.
}
\usage{
intlik3ed(start,y,K,X,distmet="ecol",covariate,theta2=NA)
}
\arguments{
  \item{start}{
    values of parameters to evaluate likelihood with.
    order is theta0, theta1, n0, theta2. If distmet="euclid" then
    dont put theta2 here
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
\item{distmet}{
either "ecol" or "euclid"
  }
  \item{covariate}{
    the covariate to be used in the cost function.
    cost(x) = exp(theta2*covariate(x))
  }
  \item{theta2}{
    you can fix the value of theta2 using this argument. In that case
    start should be length 3 and theta2 is not estimated.
    if dist="euclid" then you should have theta2=NA
  }
}
\value{
returns -1*loglikelihood
}
\author{
J. Andy Royle
}
\examples{

library("scrbook")
out<-make.EDcovariates()
covariate<-out$covariate.patchy
set.seed(2013)

N<-200
theta0<- -2
sigma<- .5
K<- 5

theta1<- 1/(2*sigma*sigma)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
theta2<-1
cost<- exp(theta2*covariate)

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
probcap<-plogis(theta0)*exp(-theta1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}
Y<-Y[apply(Y,1,sum)>0,]

n0<- N-nrow(Y)
frog1<-nlm(intlik3ed,c(theta0,theta1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="euclid",covariate=covariate,theta2=1)

frog2<-nlm(intlik3ed,c(theta0,theta1,log(n0),-.3),hessian=TRUE,y=Y,K=K,
               X=traplocs,distmet="ecol",covariate=covariate,theta2=NA)

}
