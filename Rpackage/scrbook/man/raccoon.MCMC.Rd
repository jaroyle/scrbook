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

\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}

