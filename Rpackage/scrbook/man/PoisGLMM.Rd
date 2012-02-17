\name{PoisGLMM}
\alias{PoisGLMM}
\title{
MCMC algorithm for Poisson regression with random intercept
}
\description{
MCMC algorithm for Poisson regression with random intercept
}
\usage{
PoisGLMM(y,site,mu0,sig0,a0,b0, mu_beta, sig_beta, delta_alpha, delta_beta, niter)
}
\arguments{
  \item{y}{
	Vector with observed counts
}
  \item{site}{
    vector of length(y) relating each count to a site
}
  \item{mu0}{
    size of the augmented data set
}
  \item{sig0}{
    outer coordinates of the state-space
}
  \item{a0}{
    vector containing tuning parameters for sigma (1), lam0 (2) and S (3)
}
  \item{b0}{
    Number of iterations to run the algorithm for
}
  \item{mu_beta}{
    Number of iterations to run the algorithm for
}
  \item{sig_beta}{
    Number of iterations to run the algorithm for
}
  \item{dela_alpha}{
    Number of iterations to run the algorithm for
}
  \item{delta_beta}{
    Number of iterations to run the algorithm for
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
