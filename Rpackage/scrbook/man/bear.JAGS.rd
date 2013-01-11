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

##runs the spatial null model on the Ft Drum bear data with 3 chains, 500 adaptive iterations an 2000 iterations
out<-bear.JAGS(model='SCR0', n.chains=3, n.adapt=500, n.iter=2000)

##look at summary output
summary(out)

## returns:
 XXXXXX

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JAGS }
\keyword{ covariate modeling }% __ONLY ONE__ keyword per line