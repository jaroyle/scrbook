<<<<<<< HEAD
\name{SCR0pois}
\alias{SCR0pois}
\title{
MCMC algorithm for basic spatial capture-recapture model
}
\description{
Runs the MCMC algorithm developed in Ch7 section 7.5 for a spatial capture-recapture model with Poisson encounter process
}
\usage{
SCR0pois(y, X,M, xl,xu,yl,yu,delta, niter)
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
sigma<-runif(1,0.5,5)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
=======
\name{Norm.Gibbs}
\alias{Norm.Gibbs}
\title{
Gibbs sampler for a Normal model 
}
\description{
Run a Gibbs sampler for a Normal model, y ~ Norm(mu, sigma2) with unknown mean and variance, using conjugate priors
}
\usage{
Norm.Gibbs(y, mu_0, sigma2_0, a,b,niter)
}
\arguments{
  \item{y}{
	A vector containing data
}
  \item{mu_0}{
    Mean of the Normal prior for mu
}
  \item{sigma2_0}{
    Variance of the Normal prior for mu
}
  \item{a}{
    Shape parameter of the Inverse Gamma prior for sigma2
}
  \item{mu_0}{
    Scale parameter of the Inverse Gamma prior for sigma2
}
  \item{niter}{
    Number of iterations to run the Gibbs sampler for
}
}
\value{
  A matrix with one row per iteration and two columns, for mu and sigma. Note that the output returns sqrt(sigma2) = SD.
}
\details{
  Starting values for mu and sigma are set at 1
>>>>>>> bee3f60a1d05e29e48b3940cf8627dd57c62779d
}
\author{
Rahel Sollmann, Richard Chandler, J. Andy Royle
}
\examples{

}
