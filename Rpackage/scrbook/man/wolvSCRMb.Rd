\name{wolvSCRMb}
\alias{wolvSCRMb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fit Mb (trap level) to wolverine data %%  ~~function to do ... ~~
}
\description{
fit Mb (trap level) to wolverine data  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCRMb(nb = 1000, ni = 2000, buffer = 2, M = 200)
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

wolv.Mb<-wolvSCRMb(nb=1000,ni=6000,buffer=2,M=200)

mean(wolv.Mb$sims.list$X1new>wolv.Mb$sims.list$X1obs)
mean(wolv.Mb$sims.list$X2new>wolv.Mb$sims.list$X2obs)
mean(wolv.Mb$sims.list$X3new>wolv.Mb$sims.list$X3obs)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
