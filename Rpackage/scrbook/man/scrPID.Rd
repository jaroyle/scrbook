\name{scrPID}
\alias{scrPID}
\title{
MCMC algorithm for spatial mark-resight models/ SCR models with partial individual identifiability
}
\description{
Runs the MCMC algorithm for spatial mark-resight models described in Chapter 15 
}
\usage{
scrPID<-function (n, X, y, M, obsmod = c("pois", "bern"),nmarked=c("known", "unknown"), npics,niters, 
                  xlims, ylims, a, b, inits, delta ) 
}
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
    Size of the augmented data set
}
  \item{obsmod}{
   Observation mode; either "pois" for Poisson (also called with replacement) or "bern" for Bernoulli (without replacement)
}
  \item{nmarked}{
   is the number of marked individuals known ("known") or unknown ("unknown") 
}
  \item{npics}{
    optional; if rate of individual identification of marked individuals <1, provide a vector with the number of identified pictures 
    of marked individuals and the total number of pictures ofmarked individuals (both identified and unidentified), in this order! If not provided, 
    the identification rate is assumed to be 100%
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
  \item{xlims, ylims}{
    Outer coordinates of the rectangular state-space, in the form of two vectors (x coordinates and y coordinates)
}
  \item{a,b}{
   optional; if Gamma prior is used for lam0, specify shape (a) and scale/rate (b) parameters RICHARD WHICH ONE
}
  \item{inits}{
   function to generate initial parameter values
}

  \item{delta}{
    list with tuning parameter values for sigma, lam0 and S 
}

}
\value{
  A matrix with one row per iteration and fsix columns, for sigma female, sigma male, lam0, psi, phi and N.
}

\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}

