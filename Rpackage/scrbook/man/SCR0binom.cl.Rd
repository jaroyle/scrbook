\name{SCR0binom.cl}
\alias{SCR0binom.cl}
\title{
MCMC algorithm for basic spatial capture-recapture model
}
\description{
Runs the MCMC algorithm developed in Ch7 section 7.5.1 for a spatial capture-recapture model with Binomial encounter process and clog-log link
}
\usage{
SCR0pois(y, X,M, xl,xu,yl,yu,K,delta, niter)
}
\arguments{
  \item{y}{
	Matrix with observation data, M x J
}
  \item{X}{
    2-dimensional matrix with X and Y coordinates for all J traps
}
  \item{M}{
    Size of the augmented data set
}
  \item{xl,xu,yl,yu}{
    Outer coordinates of the rectangular state-space
}
  \item{K}{
    number of occasions
}
  \item{delta}{
    vector with tuning parameter values for sigma (1), lam0 (2) and S (3)
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
}
\value{
  A matrix with one row per iteration and four columns, for sigmaa, lam0, psi and N.
}
\details{
  Starting values are chosen at random within the function as:
S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
sigma<-runif(1,5,10)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
}
\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}
\examples{

}
