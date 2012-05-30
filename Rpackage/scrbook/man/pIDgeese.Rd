\name{pIDgeese}
\alias{pIDgeese}
\title{
MCMC algorithm for analyzing Canada geese data set (Ruthledge 2012) from Chapter 15 
}
\description{
Runs the MCMC algorithm for the Canada geese example presented in Chapter 15 (spatial mark-resight model)
}
\usage{
SCR0pois(y, X, Zknown, M, Eff, niter,xl,xu,yl,yu,inits, Sex, delta)
}
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
  \item{Eff}{
    3-dimensional array coding if a given site was sampled on a given occasion (1) or not (0) (individual by site by occasion)
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
  \item{xl,xu,yl,yu}{
    Outer coordinates of the rectangular state-space
}
  \item{inits}{
   function to generate initial parameter values
}
  \item{Sex}{
   vector with coding for sex of marked individuals (females=1, males=2, unknown=NA)
}
  \item{delta}{
    list with tuning parameter values for sigma (separately for males and females), lam0 and S 
}

}
\value{
  A matrix with one row per iteration and fsix columns, for sigma female, sigma male, lam0, psi, phi and N.
}

\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}
\examples{
### run 5000 iterations of model for Canada geese
###load data
data('geesedata')  #loads a list gd with all data needed

# extract M and trap location matrix 
M=gd$M
trapmat=gd$X

#write function to generate initial values

initS<-function(){list(theta=runif(2, 0.7,1.2), lam0=runif(1,0.3,0.8), w=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
S=cbind( runif( M, min(trapmat[,1]),max(trapmat[,1] ) ), runif( M,min(trapmat[,2]),max(trapmat[,2]) ) ))}

# run model
mod<-pIDgeese (y=y, X=X, Zknown=Zknown, M=M, EffAr=EffAr, niters=5000, xl=xl,xu=xu,yl=yl,yu=yu, inits=initS(), Sex=Sex, delta=list(sig=c(0.015, 0.02), lam0=0.02, S=2))

#look at output
library(coda)
out<-mcmc(mod)
summary(out)
}
