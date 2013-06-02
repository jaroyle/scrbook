\name{make.statespace}
\alias{make.statespace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a grid of points over some prescribed region
}
\description{
Makes a grid of points of a certain size, density, and so forth. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.statespace(ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA, nx = 40, ny = NULL, buffer = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ll}{
%%     ~~Describe \code{ll} here~~
}
  \item{minx}{
%%     ~~Describe \code{minx} here~~
}
  \item{maxx}{
%%     ~~Describe \code{maxx} here~~
}
  \item{miny}{
%%     ~~Describe \code{miny} here~~
}
  \item{maxy}{
%%     ~~Describe \code{maxy} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
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
# This example is from Chapter XXXX on modeling space usage
set.seed(1234)
gr<-make.statespace(minx=1,maxx=40,miny=1,maxy=40,nx=40,ny=40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))%*%rnorm(1600)
spatial.plot(gr,z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
