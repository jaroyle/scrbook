\name{scrIPP}
\alias{scrIPP}
\title{
  A Metropolis-within-Gibbs sampler to obtain posterior distributions
  for parameters of an SCR model with an inhomogeneous binomial point
  process. 
}
\description{
  This is not meant to be generic. It assumes the state-space is rectangular.
}
\usage{
scrIPP(Z, X, M, niters, xlims, ylims, tune = rep(0.1, 4))
}
\arguments{
  \item{Z}{
    The 3-d array of encounter histories. The order of the dimensions is
    animal, trap, occassion.
}
  \item{X}{
    The matrix of trap coordinates. Should have 2 columns for easting
    and northing.
}
  \item{M}{
    The number of animals in the a
}
  \item{niters}{
    Number of MCMC iterations.
}
  \item{xlims}{
    The x-coordinates defining the state-space
}
  \item{ylims}{
    The y-coordinates defining the state-space
}
  \item{tune}{
    A vector of tuning parameters for the Metropolis steps in the
    algorithm. This should contain 3 values for the sigma, lambda_0,
    activity center proposals.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{out}{A matrix of posterior samples with one row per iteration
    and a column for the parameters: sigma, lambda_0, alpha, psi, N}
  \item{last}{The values of the latent parameters w and S for the last
    iteration. These can be used to update the model if it needs to be
    run longer}
}
\note{
  n
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
1
}
