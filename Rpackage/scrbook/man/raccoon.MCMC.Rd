\name{raccoon.MCMC}
\alias{raccoon.MCMC}
\title{
MCMC algorithm to analyze raccoon data set (ch15racdata)
}
\description{
Runs the MCMC algorithm for spatial mark-resight model includign telemetry data for raccoon example described in Chapter 15 
}
\usage{
raccoon.MCMC <- function(n,X,Eff,y,M, niters,SSp,collar, locs,mi,mall, inits) 
}
\arguments{
  \item{n}{
	Matrix with spatial counts (marked and unmarked), trap by occasion
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all J traps
}
  \item{Eff}{
    3-dimensional array with camera trap effort for each trap and occasion, expanded to all individuals (for coding ease)
}
  \item{y}{
    3-dimensional array with individual encounter histories of marked individuals (individual by trap by occasion)
}
  \item{M}{
    Size of the augmented data set
}
  \item{niters}{
    Number of iterations to run the algorithm for
}
  \item{SSp}{
    shapefile polygon of shore lines of the island that forms the entire state space
}
  \item{collar}{
    vector of length M with TRUE/FALSE depending on whether individual wears a radio-collar or not 
}
  \item{locs}{
    list, length = number of individuals with telemetry tags; each entry is a data.frame with telemetry locations
}
  \item{mi}{
   number of pictures of marked individuals that could be individually identified, per occasion
}
  \item{mall}{
   number of pictures of marked individuals, both identified and unidentified, per occasion
}
  \item{a,b}{
   optional; if Gamma prior is used for lam0, specify shape (a) and scale/rate (b) parameters RICHARD WHICH ONE
}
  \item{inits}{
   function to generate initial parameter values
}


}
\value{
XXXXX
}

\example{
#use raccoon.MCMC to run a single MCMC chain for raccoon data from Ch15

data(ch15racdata)
library(spatstat)
library(mvtnorm)
library(coda)

#create spatial window to draw random points within state space
newSSp <- as(racdata$SSp, "SpatialPolygons")
nSSp<-as.owin(newSSp)

#set up inits function
inits=function(){list(theta=runif(1,0.2,0.8), lam0=c( rep(runif(1, .05, .1), 2), 
			    rep(runif(1, .05, .1), 2), rep(runif(1, .05, .1), 2)  ),
                      xz=runifpoint(racdata$M, win=nSSp), w=rbinom(racdata$M,1,.5),
			    psi=runif(1,.2,.8) )}  

#run raccoon model
out<-raccoon.MCMC(n=racdata$n,y=racdata$y, X=racdata$X, mi=racdata$mi, 
                  mall=racdata$mall, M=racdata$M,Eff=racdata$Eff, 
                  collar=racdata$telID, locs=racdata$locs, niters=200,
                  inits=inits(), SSp=racdata$SSp)

#look at output, remove XX iterations as burn-in
mod<-mcmc(out)
summary(window(out, start=10001))
}
\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}

