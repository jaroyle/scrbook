\name{SCRgof}
\alias{SCRgof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gof analysis of the point process model%%  ~~function to do ... ~~
}
\description{
Computes a Bayesian p-value for the binomial point process of the basic
SCR model. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRgof(out, nx = 6, ny = 6, traplocs = NULL, buffer = 2, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{ This is the output list from WinBUGS or the "BUGSoutput"
    argument from the jags function in the R2jags package.
%%     ~~Describe \code{out} here~~
}
  \item{nx}{
number of breakpoints in the x-direction for creating (nx-1)*(ny-1)
total bins %%     ~~Describe \code{nx} here~~
}
  \item{ny}{
number of breakpoints in the y-direction for creating (nx-1)*(ny-1)
total bins %%     ~~Describe \code{ny} here~~
}
  \item{traplocs}{
the trap locations%%     ~~Describe \code{traplocs} here~~
}
  \item{buffer}{
size of the buffer around the trap array for computing the expected bin
counts %%     ~~Describe \code{buffer} here~~
}
  \item{Xl}{
Alternative to "buffer", specifiy the limits of the state-space%%     ~~Describe \code{Xl} here~~
}
  \item{Xu}{
Alternative to "buffer", specifiy the limits of the state-space%%     ~~Describe \code{Xl} here~~
}
  \item{Yl}{
Alternative to "buffer", specifiy the limits of the state-space%%     
}
  \item{Yu}{
Alternative to "buffer", specifiy the limits of the state-space%%     
}
}
\details{
See section 8.4.1 of the SCR book%%  ~~ If necessary, more details than the description above ~~
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

library("R2jags")
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
wsex<-wolverine$wsex

toad1<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=1)
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
