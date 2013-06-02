\name{SCRgof}
\alias{SCRgof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gof analysis of the point process model%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRgof(out, nx = 6, ny = 6, traplocs = NULL, buffer = 2, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
%%     ~~Describe \code{out} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{traplocs}{
%%     ~~Describe \code{traplocs} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

library("R2jags")
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
wsex<-wolverine$wsex

toad1<-wolvSCR0ms(y3d,traps,wsex=wsex,nb=1000,ni=2000,delta=2,M=200,model=1)
bleen<- toad1$BUGSoutput$sims.list

traplocs<-wolverine$wtraps[,2:3]
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 


a<-SCRgof(bleen,5,5,traplocs=traplocs,buffer=.4)











}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
