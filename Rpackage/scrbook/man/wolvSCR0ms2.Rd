\name{wolvSCR0ms2}
\alias{wolvSCR0ms2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model selection of encounter models%%  ~~function to do ... ~~
}
\description{
uses the Kuo Mallick idea applied to 3 models%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0ms2(nb = 1000, ni = 2000, buffer = 2, M = 200,engine="winbugs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{
%%     ~~Describe \code{y3d} here~~
}
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{wsex}{
%%     ~~Describe \code{wsex} here~~
}
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


toad<-wolvSCR0ms2(nb=500,ni=1000,buffer=2,M=121,engine="winbugs")

mprobs<- table(toad$sims.list$mod)
mprobs<-mprobs/sum(mprobs)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
