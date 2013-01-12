\name{bear.JAGS}
\alias{bear.JAGS}
\title{
analysis of the Ft. Drum black bear data using JAGS
}
\description{
Function to analyze of the Ft. Drum black bear data with a range of possible covariate models using JAGS
}
\usage{
bear.JAGS(model=c('SCR0', 'SCR0exp', 'SCRt','SCRB','SCRb', 'SCRsex', 'SCRh'), n.chains, n.adapt, n.iter)
}
\details{
bear.JAGS allows you to run a spatial the null model (SCR0), the null model with negative exponential detection function (SCR0exp),
time as factor on detection (SCRt), a global and a local trap response model (SCRB and SCRb, respectively), 
a model with sex specific detection and movement (SCRsex) 
and a heterogeneity model with 2-class mixtureswith a logit-normal mixture on detection (SCRh). 
Model results are displayed and discussed in Ch9 of Royle et al. 2013. Spatial capture Recapture.
}
\value{
JAGS model output (for details, see the rjags manual, under coda.samples())
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rahel Sollmann, Beth Gardner
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

##runs SCR model with global behavioral trap response on the Ft Drum bear data with 3 chains, 500 adaptive iterations an 2000 iterations
out<-bear.JAGS(model='SCRB', n.chains=3, n.adapt=500, n.iter=2000)

##look at summary output
summary(out)

## returns:
Iterations = 501:2500
Thinning interval = 1 
Number of chains = 3 
Sample size per chain = 2000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

           Mean       SD  Naive SE Time-series SE
D        0.1907  0.01863 0.0002405       0.001370
N      578.3960 56.49120 0.7292982       4.155007
alpha0  -2.8001  0.23204 0.0029956       0.013962
alpha2   0.8887  0.23009 0.0029704       0.012137
psi      0.8886  0.08766 0.0011317       0.006407
sigma    1.9985  0.12876 0.0016622       0.006311

2. Quantiles for each variable:

           2.5%      25%      50%      75%    97.5%
D        0.1467   0.1781   0.1952   0.2061   0.2137
N      445.0000 540.0000 592.0000 625.0000 648.0000
alpha0  -3.2728  -2.9526  -2.7945  -2.6427  -2.3526
alpha2   0.4491   0.7264   0.8867   1.0441   1.3512
psi      0.6820   0.8300   0.9101   0.9600   0.9960
sigma    1.7654   1.9101   1.9918   2.0767   2.2753


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JAGS }
\keyword{ covariate modeling }% __ONLY ONE__ keyword per line