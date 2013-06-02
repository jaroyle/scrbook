\name{modelMhBUGS}
\alias{modelMhBUGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits Model Mh using WinBUGS or JAGS
}
\description{
.....%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
modelMhBUGS(Yarr, engine = "winbugs", M = 500, ni = 10000, nb = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yarr}{
%%     ~~Describe \code{Yarr} here~~
}
  \item{engine}{
%%     ~~Describe \code{engine} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
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
data(beardata)

set.seed(2013)
jout<-modelMhBUGS(beardata$bearArray,engine="jags",ni=2000,nb=1000)

set.seed(2013)
wbout<-modelMhBUGS(beardata$bearArray,engine="winbugs",ni=2000,nb=1000)
summary(as.mcmc.list(wbout))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{heterogeneity}
\keyword{model Mh}% __ONLY ONE__ keyword per line
