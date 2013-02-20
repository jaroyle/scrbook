\name{wolvSCR0pois.fn}
\alias{wolvSCR0pois.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
wolverine SCR model with Poisson observation model
}
\description{
... more here ...%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0pois.fn(y3d, traps, nb = 1000, ni = 2000, delta = 2, M = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{
%%     ~~Describe \code{y3d} here~~
}
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
}
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
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
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
toad<-wolvSCR0.fn(y3d,traps,nb=1000,ni=6000,delta=2,M=200)
toad2<-wolvSCR0pois.fn(y3d,traps,nb=1000,ni=6000,delta=2,M=200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
