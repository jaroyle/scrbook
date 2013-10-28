\name{wolvSCR0ms}
\alias{wolvSCR0ms.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An R script that fits multiple models to the wolverine data.  }
\description{ 
Fits sex-specific models to wolverine data using the
jags() function of the R2jags package.  This is a one-trick pony -- it
works only for the wolverine data and is provided as an example. It is
not meant to be a general-purpose function.
} 
\usage{
wolvSCR0ms(nb = 1000, ni = 2000, buffer = 2, M = 200, model = 1,parm="a")
 }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{
number of MCMC iterations to be discarded as burn-in
}
  \item{ni}{
total number of MCMC iterations 
}
  \item{buffer}{
size of the buffer in scaled units (relevant to the wolverine example)%%     ~~Describe \code{delta} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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

# library("R2jags")
# library("scrbook")
# data(wolverine)
#
# These functions load libraries and data -- please read through them
# Default to run with jags()

## Run here with too FEW iterations -- increase to at least 10000 total
##
## This takes several hours to run
## 
 toad0<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=0)
 toad1<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=1)
 toad2<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=2)
 toad3<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=3)
 toad4<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=4)
 toad5<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=5)

# summary function to extract DIC info
dev.smy<-function(a){
  b<-print(a)
  dev<-mean(a$BUGSoutput$sims.list$deviance)
  c(dev,b$pD,b$DIC,a$BUGSoutput$DICbyR)
}

# make a table
rbind(dev.smy(toad0),dev.smy(toad1),dev.smy(toad2),dev.smy(toad3),dev.smy(toad4))

# summary function for parameters

parms.smy<-function(a){
  b<-print(a)
  b<- b$summary[,1:2]
  nam<-c("D","N","alpha0","alpha.sex","beta","beta.sex","sigma","psi","psi.sex","deviance")
  round(b[nam,],2)
}

# make a table
cbind(parms.smy(toad0),parms.smy(toad1),parms.smy(toad2),parms.smy(toad3),parms.smy(toad4))







}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
