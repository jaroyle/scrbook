\name{Norm.Gibbs}
\alias{Norm.Gibbs}
\title{
Gibbs sampler for Normal model
}
\description{
Gibbs sampler for Normal model y~Normal(mu, sigma)
}
\usage{
Norm.Gibbs(y,mu_0,sigma2_0,a,b,niter)
}
\arguments{
  \item{y}{
	vector of observations
}
  \item{mu_0}{
    mean for Normal prior on mu
}
  \item{sigma2_0}{
    variance for normal prior on mu
}
  \item{a}{
    shape parameter for InvGamma prior for sigma2
}
  \item{b}{
    scale parameter for InvGamma prior for sigma2
}
  \item{niter}{
    Number of iterations to run the algorithm for
}
}
\value{
  A matrix with one row per iteration and two columns, for mu and sigma. Note that the algorithm outputs the standard deviation, not sigma2.
}

\author{
Rahel Sollmann
}
\examples{
#generate normally distributed data
y<-rnorm(100, 1,2) 

#set mean and variance for vague normal prior on mu
mu_0=0
sigma2_0=100

#set shape and scale for vague inverse-gamma prior on variance
a=1
b=1

#set number of iterations
niter=10000

#run Gibbs sampler
out<-Norm.Gibbs (y=y,mu_0=mu_0,sigma2_0=sigma2_0,a=a,b=b,niter=niter)

#use package coda to look at output
library(coda)
summary(window(mcmc(out), start=501))
}
