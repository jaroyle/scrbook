\name{EDF}
\alias{EDF}
\docType{data}
\title{
Encounter Data File 
}
\description{
The trap data file describes the locations and sample occasions where uniquely identified
individuals were captured. 
}
\usage{}
\format{
}
\details{
The encounter data file is a matrix containing the essential information to fit 
a large number of SCR models. 
The data structure is also compatible with the data structure of 'secr' 
(see examples below).

The EDF should have nrows = the total number of encounter events, having 1 row
for each unique encounter.  The number of columns should be 4. The column
variables must be:  stratum (group, year or session), 
individual ID, day or occasion, trap ID.

The names of the columns don't matter (but they do for usage by 'secr' functions)


}
\source{
}
\references{
SCR book, section 5.9
}
\seealso{spiderplot, TDF, wolvSCR0}
\examples{

library("scrbook")

data(wolverine)

tdf<-wolverine$wtraps
edf<-wolverine$wcaps

str(edf)
 num [1:115, 1:4] 1 1 1 1 1 1 1 1 1 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:4] "year" "individual" "day" "trap"

dimnames(edf)<-list(NULL,NULL)  # Names don't matter

## Convert to a 3-d array suitable for analysis by BUGS or JAGS
y3d <-SCR23darray(edf,tdf)

## Make a spiderplot
## Spiderplot takes the 3-d matrix and the trap coordinates only
spiderplot(y3d, tdf[,2:3])

# Fit a basic model 
toad<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,coord.scale=10000)

}
\keyword{datasets}
