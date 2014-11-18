\name{make.EDcovariates}
\alias{make.EDcovariates}
\title{
Generate covariates for the ecological distance simulations. The
covariates
are as described in Royle, Chandler, Gazenski and Graves (2013; Ecology).
}
\description{
  ......more here....
  requires libraries gdistance and raster
}
\usage{
make.EDcovariates()
}
\arguments{
  \item{--}{
there are no arguments 
}
}
\value{
  .... a list with arguments "covariate.patchy" and
  "covariate.trend"
}
\author{
Andy Royle, aroyle@usgs.gov 
}
\examples{
library("scrbook")
out<-make.EDcovariates()
par(mfrow=c(1,2))
plot(out$covariate.patchy)
plot(out$covariate.trend)

}
