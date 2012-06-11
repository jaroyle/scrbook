\name{bbsdata.rda}
\alias{bbsdata.rda}
\docType{data}
\title{
Muntjac line transect survey data from Samba
}
\description{
Mourning dove counts
made along North American Breeding Bird Survey (BBS) routes in
Pennsylvania, USA. Generally, a route consists of 50 stops separated by 0.5
miles. This data set contains route total counts and the sample location is the center point of the BBS
route.
}
\usage{data(bbsdata)}
\format{
  List containing two data frames. bbsdata$counts: Data frame with 103 rows and 37 columns; column 1-4 contain route ID, longitude and latitude, and forest cover (standardized); 
  column 5-37 contain yearly counts (1966-1998); bbsdata$habitat: data frame with 1200 rows and 19 columns) with raster habitat data for the state of Pennsylvania
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://www.pwrc.usgs.gov/bbs
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(bbsdata)
#look at mourning dove counts
y<-bbsdata$counts[,29]  # pick out 1990
notna<-!is.na(y)
y<-y[notna]
#histogram of mourning dove counts in yr. 1990
hist(y)

#produce spatial plot of forest cover for state of PA
library(maps)
habdata<-bbsdata$habitat
map('state',regions="penn",lwd=2)
spatial.plot(habdata[,2:3],habdata[,"dfor"],cx=2)
map('state',regions="penn",lwd=2,add=TRUE
}
\keyword{datasets}
