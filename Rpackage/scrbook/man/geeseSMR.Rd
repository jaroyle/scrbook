\name{geeseSMR}
\alias{geeseSMR}
\title{
MCMC for Canada geese mark-resight analysis from Chapter 19
}
\description{
MCMC sampler for analysis of spatial mark-resight dataset of 654 Canada geese (Rutledge 2012), 
collected in North Carolina between mid July to the end of October 2008. Total number of marked is unknown. Sigma is modeled sex-specific.
}
\usage{geeseSMR<-function(y, X, Zknown, M, EffAr, niters,inits, Sex, tune,SSp)}

\arguments{
  \item{y}{
	Matrix with spatial counts (marked and unmarked), trap by occasion
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all J traps
}
  \item{Zknown}{
    3-dimensional array with individual encounter histories of marked individuals (individual by trap by occasion)
}
  \item{M}{
    Size of the augmented data set
}
  \item{EffAr}{
   3-D array, M by J by K with 0/1 coding to indicate whether a given resighting location was visited on a particular occasion
}
  \item{nmarked}{
   is the number of marked individuals known ("known") or unknown ("unknown") 
}

  \item{niters}{
    Number of iterations to run the algorithm for
}
  \item{inits}{
   function to generate initial parameter values
}
  \item{Sex}{
   Vector with coding for sex of resighted individuals; 1=female, 2=male, NA=unknown
}
  \item{tune}{
    list with tuning parameter values for sigma (female and male), lam0 and S 
}
  \item{SSp}{
   SpatialPolygon of the state-space (resighting locations plus 3 km buffer)
}

}

\value{
  A matrix with one row per iteration and fix columns, for sigma female, sigma male, lam0, psi, phi and N.
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
