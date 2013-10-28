\name{scr2secr}
\alias{scr2secr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts trap deployment file with operation history to secr format.
}
\description{
Converts trap deployment file with operation history to secr
format. Requires the function read.traps() from the 'secr' package. 

}
\usage{
scr2secr(scrtraps, type = "proximity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrtraps}{
Standard trap deployment file (TDF) with columns 'trap id' 'x-coord'
'y-coord' and columns [,4:(noccasions + 3)] being the trap operation
information. %%     ~~Describe \code{scrtraps} here~~
}
  \item{type}{
The type of encounter device recobnized by 'secr'. e.g., proximity
detectors correspond to type="proximity"  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book, sec. 6.5.2 %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# required package 'secr'
library("secr")

# load package 'scrbook' and the wolverine data
library("scrbook")
data("wolverine")

# grab the trap deployment file
traps<-as.matrix(wolverine$wtraps)   

# provide suitable names required by 'secr'
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

# convert to a 'secr' trapfile
trapfile<-scr2secr(scrtraps=traps,type="proximity")

# grab the encounter data file 
wolv.dat<-wolverine$wcaps

# provide suitable names required by 'secr' and convert to a dataframe
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)

# convert to a capture history file usable by 'secr' 
wolvcapt<-make.capthist(wolv.dat,trapfile,fmt="trapID",noccasions=165)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
