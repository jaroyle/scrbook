\name{area}
\alias{area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
computes area of a polygon
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
area(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (x, y) 
{
    if (missing(y)) {
        if (is.matrix(x) && ncol(x) == 2) {
            y <- x[, 2]
            x <- x[, 1]
        }
        else if (!is.null(x$x) && !is.null(x$y)) {
            y <- x$y
            x <- x$x
        }
    }
    x <- c(x, x[1])
    y <- c(y, y[1])
    i <- 2:length(x)
    return(0.5 * sum(x[i] * y[i - 1] - x[i - 1] * y[i]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
