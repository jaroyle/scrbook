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
(optional) coordinates around which to locate the state-space %%     ~~Describe \code{ll} here~~
}
  \item{minx}{
lower x bound  %%     ~~Describe \code{minx} here~~
}
  \item{maxx}{
upper x bound %%     ~~Describe \code{maxx} here~~
}
  \item{miny}{
lower y bound %%     ~~Describe \code{miny} here~~
}
  \item{maxy}{
upper y bound %%     ~~Describe \code{maxy} here~~
}
  \item{nx}{
number of points in x-direction %%     ~~Describe \code{nx} here~~
}
  \item{ny}{
(should be missing) ny is chosen to produce y-spacing consistent with
x-spacing implied by choice of nx %%     ~~Describe \code{ny} here~~
}
  \item{buffer}{
optional buffer in relative terms (e.g., .1 chooses a buffer 10 percent
of (maxx-minx) ) %%     ~~Describe \code{buffer} here~~
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
  Andy Royle, aroyle@usgs.gov 
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
z<-t(chol(V))\%*\%rnorm(1600)
spatial.plot(gr,z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
