\name{SCR0poisSSp}
\alias{SCR0poisSSp}
\title{
MCMC algorithm for basic spatial capture-recapture model with state-space defined by a shape file
}
\description{
Runs the MCMC algorithm developed in Ch17 for a spatial capture-recapture model with Poisson encounter process, but with a state-space that is defined by a shapefile, not a rectangle
}
\usage{
SCR0poisSSp(y, X,M,  SSp, delta,niter)
}
\arguments{
  \item{y}{
	Matrix with observation data, M x J
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all J traps, in km spacing
}
  \item{M}{
    Size of the augmented data set
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
  \item{delta}{
    vector with tuning parameter values for sigma (1), lam0 (2) and S (3)
}
  \item{SSp}{
    state-space polygon, type "SpatialPolygons", in UTM
}

}
\value{
  A matrix with one row per iteration and four columns, for sigma, lam0, psi and N.
}
\details{
  Starting values are chosen at random within the function as:
S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
sigma<-runif(1,0.5,5)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
Function requires the packages sp, spatstat...
The code is currently set up for the format "SpatialPolygons" of SSp
and will run into problems if SSp is a "SpatialPolygonsDataFrame". The
code can easily be modified to fit a
\code{
"SpatialPolygonsDataFrame" by substituting "
if(is.na(SinPoly[i])==FALSE) { " with 	" if(is.na(SinPoly[i,1])==FALSE)
  { "
    }
}
\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}
\examples{

}
