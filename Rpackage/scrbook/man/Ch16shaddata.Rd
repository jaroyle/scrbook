\name{Ch16shaddata}
\alias{Ch16shaddata}
\docType{data}
\title{
  American shad data used in Chapter 16
}
\description{
  Data collected by Raabe and colleagues in the Little River, North Carolina, described in Raabe et al. (2014). American shad (Alosa sapidissima) were tagged 
with passive integrated transponders (PIT) at a weir near the mouth the river and passively monitored with an upstream array of PIT antennas in 2010.  In this example,
the weir and PIT antennas are treated the same; however, Raabe et al. (2014) implement a more specific model for the data.
}

\usage{data(Ch16shaddata)}
\format{
  The format is:
 List of 2
 $ shad :'data.frame':  3406 obs. of  4 variables:
  ..$ TagID: int [1:3406] 161210451 161210451 161210451 161210451 161210451 161210451 161210451 161210451 161210451 161210451 ...
  ..$ Week : int [1:3406] 1 1 1 1 1 1 1 4 4 5 ...
  ..$ RKM  : num [1:3406] 7.7 7.7 13.4 13.4 13.4 13.4 13.4 45.3 56.4 13.4 ...
  ..$ state: int [1:3406] 1 1 1 1 1 1 1 2 2 1 ...
 $ first: num [1:315] 1 7 6 6 6 6 6 6 5 5 ...


}
\details{
 TagID describes the fish id number, Week is the week (1-12) that the fish was captured, RKM is the river kilometer that the fish was captured in and corresponds to
the antenna location, state indicates the geographic state (1 or 2 for upper or lower reach) for the multi-state model, 
and first is the week in which the fish was first captured (length is 315 because there are 315 unique fish in this dataset).
}

\references{
Raabe, J. K., Gardner, B., & Hightower, J. E. (2014). A spatial capture-recapture model to estimate fish 
survival and location from linear continuous monitoring arrays. Canadian Journal of Fisheries and Aquatic Sciences, 71(1), 120-130.
}
}
\examples{
data(Ch16shaddata)
}
\keyword{datasets}
