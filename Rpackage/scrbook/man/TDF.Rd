\name{TDF}
\alias{TDF}
\title{
Trap deployment file 
}
\description{
The trap deployment file describes the locations and operation schedule of traps or data collection locations which will be used in the fitting of SCR models
}
\usage{}
\format{
}
\details{
The trap deployment file (TDF) is a matrix or dataframe which describes the locations and operation schedule of a trap array.  
For running models in WinBUGS or JAGS it is best to coerce a TDF into a matrix.

A typical TDF should have dimension ntraps x (noccasions + 3)  where ntraps = the number of traps ever operated and 
noccasions = K = the number of sampling occasions in the study. e.g., for a camera trapping study K might be the number of nights
the cameras were operated. For a hair snag study of bears K might be the number of weekly checks of the hair snags.

Column 1 of the TDF should be the trap identity (an integer)

Column 2 of the TDF should be the x-coordinates of the traps. It is best to have these as UTMs divided by 1000

Column 3 of the TDF should be the y-coordinates of the traps.

Each column from 4 to (noccasions + 3) should be a vector of 0's and 1's indicating whether the trap was operational for that
sampling occasion (value = 1) or not operational (value = 0). 

The TDF is compatible with trap information file used by Efford's 'secr' package (see example below)
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
section 5.9 of the SCR book
}
\seealso{secr_wolverine}
\author{Andy Royle, aroyle@usgs.gov}
\examples{

library(scrbook)
data(wolverine)
TDF<-wolverine$wtraps

## To be safe always convert this to a matrix
TDF<-as.matrix(TDF)

## Let's take a look. The wolverine study involved 37 camera traps and 
##   they were operated over a 168 day period
str(TDF)

 num [1:37, 1:168] 1 2 3 4 5 6 7 8 9 10 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:168] "1:37" "Easting" "Northing" "1" ...


## Do a Bayesian analysis in WinBUGS 

y3d <-SCR23darray(wolverine$wcaps,TDF)

# produce 'nind' x 'ntraps' 2-d matrix

y2d <- apply(y3d,c(1,2),sum)

# Execute the R/BUGS script 
toad<-wolvSCR0(y3d, TDF ,nb=1000,ni=2000,buffer=1,M=100,coord.scale=10000)


# Part of the example from secr_wolverine:

library("secr")
traps<-as.matrix(wolverine$wtraps)   
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

trapfile2<-scr2secr(scrtraps=traps,type="proximity")

wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)
wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)

# now use wolvcapt2 which has trap operation information
wolv.secr0<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000)
 




}
\keyword{datasets}
