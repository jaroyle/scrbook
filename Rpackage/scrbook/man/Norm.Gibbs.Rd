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

}
