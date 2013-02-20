\name{SCRdensity}
\alias{SCRdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a spatial capture-recapture density plot %%  ~~function to do ... ~~
}
\description{
requires MCMC output in the form of histories of activity center
coordinates and data augmentation variables
}
\usage{
SCRdensity(obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL, scalein = 100, scaleout = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
list with elements Sx Sy and w each is "niters x M" matrices
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{Xl}{
%%     ~~Describe \code{Xl} here~~
}
  \item{Xu}{
%%     ~~Describe \code{Xu} here~~
}
  \item{Yl}{
%%     ~~Describe \code{Yl} here~~
}
  \item{Yu}{
%%     ~~Describe \code{Yu} here~~
}
  \item{scalein}{
%%     ~~Describe \code{scalein} here~~
}
  \item{scaleout}{
%%     ~~Describe \code{scaleout} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
# this takes 341 seconds on a standard CPU circa 2011
unix.time(out<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,keepz=TRUE))

Sx<-out$sims.list$s[,,1]
Sy<-out$sims.list$s[,,2]
w<- out$sims.list$z
obj<-list(Sx=Sx,Sy=Sy,z=z)
SCRdensity(obj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
