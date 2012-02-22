\name{SCR0bayesDss}
\alias{SCR0bayesDss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit SCR0 with a discrete state-space%%  ~~function to do ... ~~
}
\description{
  same as above.
  requires simSCR0.fn
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR0bayesDss(data, ng = 10, engine = "jags", ni = 2000, nb = 1000, M = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{ng}{
%%     ~~Describe \code{ng} here~~
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
# test with a small run -- this takes a few minutes
data<-simSCR0.fn(discard0=TRUE,sd=2013)
out1<-SCR0bayesDss(data,ng=8,M=200,engine="jags",ni=2000,nb=1000)
summary(out1$out)
out2<-SCR0bayesDss(data,ng=8,M=200,engine="winbugs",ni=2000,nb=1000)
summary(as.mcmc.list(out2$out))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
