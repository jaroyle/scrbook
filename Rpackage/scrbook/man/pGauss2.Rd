\name{pGauss2}
\alias{pGauss2}
\title{
pFauss2
}
\description{
blah blah
}
\usage{
pGauss2(parms, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parms}{
%%     ~~Describe \code{parms} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
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

## The function is currently defined as
function (parms, D) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    lp <- parms[1] - (1/(2 * parms[2] * parms[2])) * D * D
    p <- 1 - exp(-exp(lp))
    p
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
