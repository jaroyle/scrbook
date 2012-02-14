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
}
\author{
Rahel Sollmann
}
\examples{

}
