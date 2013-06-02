\name{int2d}
\alias{int2d}
\title{
Crude two-dimensional integration over the unit square.
}
\description{
This is meant to be used in Ch9
}
\usage{
int2d(alpha, delta = 0.02)
}
\arguments{
  \item{alpha}{
    The coefficient for a simple log-linear model with one spatial covariate
}
  \item{delta}{
    The grid size to use for approximating the integral
}
}
\details{
  This is meant to illustrate the proceedure, but it is not efficient
  when high levels of accuracy are needed.
}
\value{
  The value of the integral
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(alpha, delta=0.02) {
  z <- seq(-1+delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(alpha*elev.fn(S)) * cell.area)
  }
}
