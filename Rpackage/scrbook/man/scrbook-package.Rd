\name{scrbook-package}
\alias{scrbook-package}
\alias{scrbook}
\docType{package}
\title{
  Companion to the book "Spatial Capture Recapture"
}
\description{
This package contains a large number of R and BUGS/JAGS cripts, data
sets, and utility functions for conducting analyses contained in the
book Spatial Capture-Recapture (Elsevier/Academic Press, 2013). 
}
\details{
\tabular{ll}{
Package: \tab scrbook\cr
Type: \tab Package\cr
Version: \tab 0.xx\cr
Date: \tab 2013-10-26\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

Data format: most of the functions and scripts contained in the scrbook
package operate on standard trap deployment file (TDF) and encounter
data file (EDF) formats or 2-d or 3-d arrays. These are easily produced
from ordinary spreadsheet files with help from a number of utility
functions
including:

SCR23darray  will convert EDF and TDF information to 3-d array useful
for BUGS/JAGS

array2secr will convert a 3-d array into the EDF format useful in secr
and SCRbayes

scr2secr converts the TDF into a secr traps file.


~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
Andy Royle, Richard Chandler, Rahel Sollmann and Beth Gardner. 

Maintainer: Andy Royle, aroyle@usgs.gov
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
