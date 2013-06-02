\name{wolvSCR0gof}
\alias{wolvSCR0gof}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gof for wolverine %%  ~~function to do ... ~~
}
\description{
fits a model but saves s and w output%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0gof(nb = 1000, ni = 2000, buffer = 2, M = 200, model = 1, parm = "a")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{
%%     ~~Describe \code{nb} here~~
}
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{parm}{
%%     ~~Describe \code{parm} here~~
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

toad1<-wolvSCR0gof(nb=1000,ni=6000,buffer=2,M=200,model=0)

bugsout <- toad1$BUGSoutput$sims.list

traplocs<-wolverine$wtraps[,2:3]
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 

set.seed(2013)
# crop tight to trap grid. Buffer = .4 (Remember the buffer for fitting
# was 2 units.
SCRgof(bugsout.h,5,5,traplocs=traplocs,buffer=.4)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
