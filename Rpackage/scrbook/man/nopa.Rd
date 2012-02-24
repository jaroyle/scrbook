\name{nopa}
\alias{nopa}
\docType{data}
\title{
  The Northern Parula data.
}
\description{
  Point count data from 105 points surveyed 3 times each during 2006.
}
\usage{data(nopa)}
\format{
  The format is:
  List of 7
  $ y     : num [1:105, 1:3] 0 0 0 0 0 0 0 0 0 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:105] "4005" "3055" "3005" "2055" ...
  .. ..$ : chr [1:3] "v1" "v2" "v3"
  $ X     : num [1:105, 1:2] 5 5 5 5 5 5 5 6 6 6 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "x" "y"
  $ xSide : num 24
  $ ySide : num 16
  $ M     : num 100
  $ nTraps: int 105
  $ nReps : int 3
}
\details{
  \code{y} is the nTraps x nReps (J x K) matrix of counts.
  \code{X} is the Jx2 matrix of point count locations.
  \code{xSide,ySite} describe the state-space (the lengths of the sides
  of the rectangle containing the point count locations).
  \code{M} is the maximum possible value of N used to for data augmentation.
  }
\note{
  Special thanks to Deanna Dawson for collecting the data.
}
\examples{
data(nopa)
## maybe str(nopa) ; plot(nopa) ...
}
\keyword{datasets}
