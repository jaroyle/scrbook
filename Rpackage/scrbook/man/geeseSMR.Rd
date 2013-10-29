\name{geeseSMR}
\alias{geeseSMR}
\title{
MCMC for Canada geese mark-resight analysis from Chapter 19
}
\description{
MCMC sampler for analysis of spatial mark-resight dataset of 654 Canada geese (Rutledge 2013), 
collected in North Carolina between mid July to the end of October 2008. Total number of marked is unknown. Sigma is modeled sex-specific.
}
\usage{geeseSMR(n, X, y, M, mmax, obsmod = c("pois", "bern"),niters, npics, Sex,EffAr,
    inits, delta, SSp )}

\arguments{
  \item{n}{
	Matrix with spatial counts (marked and unmarked), trap by occasion
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all J traps
}
  \item{y}{
    3-dimensional array with individual encounter histories of marked individuals (individual by trap by occasion)
}
  \item{M}{
    Size of the augmented unmarked data set
}
  \item{mmax}{
    Size of the augmented data set, maximum number ever marked before resighting
}
  \item{obsmod}{
   Observation model; either "pois" for Poisson (also called with replacement) or "bern" for Bernoulli (without replacement)
}
  \item{niters}{
    Number of iterations to run the algorithm for
}
  \item{npics}{
    optional; if rate of individual identification of marked individuals <1, provide a vector with the number of identified pictures 
    of marked individuals and the total number of pictures of marked individuals (both identified and unidentified), in this order! If not provided, 
    the identification rate is assumed to be 100%
}
  \item{Sex}{
   Vector with coding for sex of resighted individuals; 1=female, 2=male, NA=unknown
}
  \item{EffAr}{
   3-D array, M+mmax by J by K with 0/1 coding to indicate whether a given resighting location was visited on a particular occasion;
	values are the same across all M+mmax
}
  \item{inits}{
   function to generate initial parameter values for S, sigma (female and male), lam0, psi, psim and phi 
}

  \item{delta}{
    list with tuning parameter values for sigma (female and male), lam0 and S 
}
  \item{SSp}{
   SpatialPolygon of the state-space (resighting locations plus 4.5 km buffer)
}

}

\value{
  A matrix with one row per iteration and eight columns, for sigma female, sigma male, lam0, c,  psi, phi, m and N (sum of marked and unmarked).
}

\references{
Rutledge, M.E., 2013. Impacts of Resident Canada Geese in a Suburban Environment. PhD Dissertation, North Carolina State University,
Raleigh, NC, USA.
}

\examples{
#load data and required packages for analysis
data(geesedata)
library(spatstat)
library(maptools)
library(coda)

#set up data and initial values
X<-newgeesedata$X
 xxx<-convexhull.xy(vertices(dilation.owin(convexhull.xy(X), 4.5)))
xxy<-as(xxx, "SpatialPolygons")

M<-newgeesedata$M
mmax<-newgeesedata$mmax
EffAr<-aperm(array(newgeesedata$Eff, c(dim(X)[1], dim(newgeesedata$y)[3],M+mmax)),c(3,1,2) )
delta=c(0.015, 0.008,0.8)

inits<-function(){list(sigma=runif(2, 0.9,1.1), lam0=runif(1,0.28,0.35), z=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
 psim=runif(1, 0.3, 0.8),  phi=runif(1, 0.3, 0.8), S=runifpoint(M+mmax,xxx)   )}



#run 5000 iterations - will take several hours

out<-geeseSMR(n=newgeesedata$n, X=X,y=newgeesedata$y, M=M, mmax=mmax, obsmod="bern", niters=5000, Sex=newgeesedata$Sex,
			EffAr=EffAr,inits=inits(),delta=delta, SSp=xxy  )

#lookat summary output
summary(window(mcmc(out), start=1001))
}
