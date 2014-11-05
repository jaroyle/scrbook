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


I. Data formatting and manipulation:

Most of the functions and scripts contained in the scrbook
package operate on standard trap deployment file (?TDF) and encounter
data file (?EDF) formats or 2-d or 3-d arrays. These are easily produced
from ordinary spreadsheet files with help from a number of utility
functions which include:

SCR23darray: will convert EDF and TDF information to 3-d array
  useful for BUGS/JAGS 

array2secr: will convert a 3-d array into the EDF format useful in secr
and SCRbayes 

scr2secr: converts the TDF into a secr traps file. 

II. Fitting some basic models in WinBUGS/JAGS:

\tab Simulating data and fitting a model: SCR0bayes 

\tab Wolverine camera trapping data: wolvSCR0 

\tab Multi-session models: SCRovenbird

III. Fitting some based models in the R package secr (M. Efford's package):

\tab Basic model fitted to wolverine data: secr_wolv

\tab Multi-session models: SCRovenbird

III. Advanced models

\tab Space-usage/RSF/telemetry: RSF_example

\tab Ecological distance

\tab Open models: ovenbirds.js, SCRovenbird

\tab Mark-resight SCR: raccoon.MCMC

IV. Summary analysis of output:

Effective sample area

Density plots


}
\author{
Andy Royle, Richard Chandler, Rahel Sollmann and Beth Gardner. 

Maintainer: Andy Royle, aroyle@usgs.gov
}
\references{
  Royle, J. A., R. B.  Chandler, R. Sollmann, and B. Gardner. 2014.
  Spatial capture-recapture. Academic Press.
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
## 
## Demo to run a basic SCR model fast.
##
library(scrbook)
# The wolverine data is a basic camera trapping data set from an array
# of 37 camera traps that were run over 165 days. You an load the
# data and inspect the structure as follows:
 
data(wolverine)

# "wcaps" is an encounter data file (EDF), see ?EDF for format details
> wolverine$wcaps[1:4,]
     year individual day trap
[1,]    1          2 127    1
[2,]    1          2 128    1
[3,]    1          2 129    1
[4,]    1         18 130    1

# "wtraps" is a trap deployment file (TDF), see ?TDF for format details
> wolverine$wtraps[1:4,1:10]
  1:37 Easting Northing 1 2 3 4 5 6 7
1    1  632538  6316012 0 0 0 0 0 0 0
2    2  634822  6316568 1 1 1 1 1 1 1
3    3  638455  6309781 0 0 0 0 0 0 0
4    4  634649  6320016 0 0 0 0 0 0 0

# From this data structure it is easy to run models either in
#  WinBUGS/JAGS or in the R package 'secr' by M. Efford.
# we have an R script 'wolvSCR0' to do this (see ?wolvSCR0).
# To fit a basic SCR model do this:

tdf <- as.matrix(wolverine$wtraps)
# convert the data to the 3-d array
y3d <-SCR23darray(wolverine$wcaps,tdf)
toad<-wolvSCR0(y3d,tdf,nb=1000,ni=2000,buffer=1,M=100,coord.scale=10000)

###
### Now let's fit a similar model in 'secr'
###
library(secr)
# Provide suitable names required by 'secr'
dimnames(tdf) <- list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

# Convert to a 'secr' trapfile
trapfile<-scr2secr(scrtraps=tdf,type="proximity")

# grab the encounter data file 
wolv.edf <- wolverine$wcaps

# provide suitable names required by 'secr' and convert to a dataframe
dimnames(wolv.edf)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.edf <- as.data.frame(wolv.edf)

# convert to a capture history file usable by 'secr' 
wolvcapt<-make.capthist(wolv.edf, trapfile,fmt="trapID",noccasions=165)


~~ simple examples of the most important functions ~~
}
