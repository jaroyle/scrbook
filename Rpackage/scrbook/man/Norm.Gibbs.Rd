\name{SCR0pois}
\alias{SCR0pois}
\title{
MCMC algorithm for basic spatial capture-recapture model
}
\description{
Runs the MCMC algorithm described in Ch7. section 7.5 for spatial capture-recapture model with Poisson encounter process
}
\usage{
SCR0pois(y,X,M, xl,xu,yl,yu,delta, niter)
}
\arguments{
  \item{y}{
	Matrix containing observation data, dimensions: M x J
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all Jtraps
}
  \item{M}{
    size of the augmented data set
}
  \item{xl, xu, yl, yu}{
    outer coordinates of the state-space 
}
  \item{delta}{
    vector containing tuning parameters for sigma (1), lam0 (2) and S (3)
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
}
\value{
  A matrix with one row per iteration and four columns, for sigma, lam0, psi and N. 
}
\details{
  Starting values are chosen at random within the function:
S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
sigma<-runif(1,0.5,5)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
}
\author{
Rahel Sollmann
}
\examples{

}
