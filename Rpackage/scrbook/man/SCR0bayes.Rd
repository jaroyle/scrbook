\name{SCR0bayes}
\alias{SCR0bayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit SCR model using JAGS or WinBUGS %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR0bayes(dataobj, M = 200, engine = "jags", ni = 2000, nb = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataobj}{
%%     ~~Describe \code{dataobj} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{engine}{
%%     ~~Describe \code{engine} here~~
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
data<-simSCR0.fn(discard0=TRUE,sd=2013)
out1<-SCR0bayes(data,M=200,engine="jags",ni=2000,nb=1000)

out2<-SCR0bayes(data,M=200,engine="winbugs",ni=2000,nb=1000)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
