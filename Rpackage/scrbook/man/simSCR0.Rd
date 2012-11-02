\name{simSCR0}
\alias{simSCR0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate some SCR data %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simSCR0(discard0 = TRUE, array3d = FALSE, rnd = 2013)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{discard0}{
%%     ~~Describe \code{discard0} here~~
}
  \item{array3d}{
%%     ~~Describe \code{array3d} here~~
}
  \item{rnd}{
%%     ~~Describe \code{sd} here~~
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
data<-simSCR0(discard0=TRUE,rnd=2013)


##out1<-SCR0bayes(data,M=200,engine="jags",ni=2000,nb=1000)
# Fit the model using WinBUGS:
out2<-SCR0bayes(data,M=200,engine="winbugs",ni=2000,nb=1000)






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
