\name{geesedata}
\alias{geesedata}
\docType{data}
\title{
Canada geese mark-resight data from Chapter 19
}
\description{
Spatial mark-resight dataset of 654 Canada geese (Rutledge 2012), collected in North Carolina between mid July to the end of October 2008. 
}
\usage{data(gd)}
\format{
  The format is:
List of 7
 $ n      : num [1:57, 1:21] 0 23 0 0 0 0 0 0 0 0 ...
 $ X      : num [1:57, 1:2] 3988 3997 3993 3995 3993 ...
 $ trapmat: num [1:87, 1:2] 3997 3988 3997 3993 3995 ...
 $ y      : num [1:654, 1:57, 1:21] 0 0 0 0 0 0 0 0 0 0 ...
 $ M      : num 5700
 $ EffAr  : num [1:5700, 1:57, 1:21] 0 0 0 0 0 0 0 0 0 0 ...
 $ Sex    : num [1:654] 1 2 1 2 1 1 1 2 1 1 ...

}

\source{
Rutledge, M. 2012. GET TITLE
}
\references{
Rutledge, M. 2012. GET TITLE
}
\examples{
#load data and required packages for analysis
data(geesedata)
library(spatstat)
library(maptools)
library(coda)

#set up data and initial values
trapmat<-geesedata$trapmat
 xxx<-convexhull.xy(vertices(dilation.owin(convexhull.xy(trapmat), 3)))
xxy<-as(xxx, "SpatialPolygons")
M=geesedata$M

initS<-function(){list(theta=runif(2, 0.9,1.1), lam0=runif(1,0.28,0.35), w=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
S=runifpoint(M,xxx)   )}

#run 5000 iterations - will take several hours
mod1<-geeseSMR(y=geesedata$n, X=geesedata$X, Zknown=geesedata$y, M=M, EffAr=geesedata$EffAr, niters=5000,inits=initS(),
		Sex=geesedata$Sex, tune=list(sig=c(0.015, 0.02), lam0=0.02, S=2),SSp=xxy)

#lookat summary output
summary(window(mcmc(mod1), start=1001))
}
