\name{wolvESA}
\alias{wolvESA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
computes effective sample area plot for wolverine data %%  ~~function to do ... ~~
}
\description{
computes effective sample area plot for wolverine data
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvESA(noargs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noargs}{
%%     ~~Describe \code{noargs} here~~
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
function (noargs = TRUE) 
{
    library("scrbook")
    data(wolverine)
    traps <- wolverine$wtraps
    y3d <- SCR23darray(wolverine$wcaps, wolverine$wtraps)
    g <- as.matrix(wolverine$grid2)
    SCRsmy(y3d)
    Sgrid <- g
    traplocs <- as.matrix(traps[, 2:3])
    mingridx <- min(traplocs[, 1])
    mingridy <- min(traplocs[, 2])
    traplocs[, 1] <- traplocs[, 1] - min(traplocs[, 1])
    traplocs[, 2] <- traplocs[, 2] - min(traplocs[, 2])
    traplocs <- traplocs/10000
    ntraps <- nrow(traplocs)
    MASK <- traps[, 4:ncol(traps)]
    ndays <- apply(MASK, 1, sum)
    Sgrid[, 1] <- Sgrid[, 1] - mingridx
    Sgrid[, 2] <- Sgrid[, 2] - mingridy
    Sgrid <- Sgrid/10000
    Dmat <- e2dist(traplocs, Sgrid)
    probs <- rep(1/nrow(Sgrid), nrow(Sgrid))
    p0 <- 0.05
    sigma <- 0.62
    D <- e2dist(traplocs, Sgrid)
    netp <- rep(NA, nrow(Sgrid))
    for (i in 1:nrow(Sgrid)) {
        d2 <- (traplocs[, 1] - Sgrid[i, 1])^2 + (traplocs[, 2] - 
            Sgrid[i, 2])^2
        pvec <- p0 * exp(-(1/(2 * sigma * sigma)) * d2)
        netp[i] <- 1 - prod((1 - pvec)^ndays)
    }
    y <- netp
    x <- Sgrid
    par(mar = c(3, 3, 3, 6))
    plot(x, pch = " ")
    nc <- as.numeric(cut(y, 10))
    cc <- topo.colors(10)
    points(x, pch = 20, col = cc[nc], cex = 1)
    image.scale(y, col = cc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
