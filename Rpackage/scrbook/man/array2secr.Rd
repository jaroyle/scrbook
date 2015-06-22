\name{array2secr}
\alias{array2secr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts detection histories from BUGS format to secr format
}
\description{
Converts 3-dimensional arrays with individual 0/1 detection histories (BUGS format) to 2-dimensional secr format.
Presently does not work for 4-d arrays (multi-sessions)
}
\usage{
array2secr(y, session, detector="proximity")
}

\arguments{
  \item{y}{
 3-dimensional array of individuals x traps x occasions with 0/1 coding for captures
}
  \item{session}{
Vector with session information for each individual; optional; if left empty all individuals are assumed to be captured in session 1. See details for use of session.
}
  \item{detector}{
Character string specifying whether detection data in y is binary ("proximity") or count data ("count"). 
}
}
\details{
Session in secr can be used for robust design type analyses. If you have several sessions in your data set, 
the argument 'session' allows you to specify those by providing a vector of length dim(y)[1] with session information for each individual in y. 
Alternatively, sessions can be used in secr to include individual categorical covariates such as sex in the model (see http://www.otago.ac.nz/density/pdfs/secr-overview.pdf).
See the example below on how to use 'session' in array2secr to specify sex of the detected individuals.
array2secr currently only supports the detector types "count" and "proximity" (secr naming convention).For more information on detector types used 
in secr, see secr documentation (help(secr)).
}
\value{
2-dimensional secr detection matrix with one row for each capture and the columns (using secr naming conventions) "Session","ID","Occasion" and "trapID" 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rahel Sollmann  %%  ~~who you are~~
}
\note{
Currently array2secr only supports binary and count detection data for the detector types "proximity" and "count", respectively.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
scr2secr, SCR23darray, SCR24darray
}
\examples{

library("secr")
library("scrbook")
data("beardata")

y<-beardata$bearArray
sex<-beardata$sex

caps<-array2secr(y, detector="proximity") #no sex specification, all individuals get session=1
caps.sex<-array2secr(y, session=sex, detector="proximity") #each individual has sex assigned 


## Now reverse the steps, use SCR23darray to produce the 3darray from "caps"

## trapmat is not a trap deployment file format so make one
traps<-beardata$trapmat
bear.tdf<-cbind(1:nrow(traps), traps, matrix(1,ncol=8,nrow=nrow(traps)))
bear.edf<- caps

# should look exactly like beardata$bearArray
y3d<- SCR23darray(bear.edf,bear.tdf)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
