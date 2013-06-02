\name{make.EDcovariates}
\alias{make.EDcovariates}
\title{
generate covariates for the ecological distance simulations
}
\description{
......more here
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
Andy Royle
}
\examples{
library("scrbook")
out<-make.EDcovariates()
par(mfrow=c(1,2))
plot(out$covariate.patchy)
plot(out$covariate.trend)

}
