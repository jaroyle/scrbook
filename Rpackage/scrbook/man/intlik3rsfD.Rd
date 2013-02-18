\name{intlik3rsfD}
\alias{intlik3rsfD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
scr+rsf likelihood%%  ~~function to do ... ~~
}
\description{
does as described in Royle et al. mee paper%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlik3rsfD(start = NULL, y = y, K = NULL, X = traplocs, ztrap, G, ntel = NULL, zall = NULL, stel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{ztrap}{
%%     ~~Describe \code{ztrap} here~~
}
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{ntel}{
%%     ~~Describe \code{ntel} here~~
}
  \item{zall}{
%%     ~~Describe \code{zall} here~~
}
  \item{stel}{
%%     ~~Describe \code{stel} here~~
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
function (start = NULL, y = y, K = NULL, X = traplocs, ztrap, 
    G, ntel = NULL, zall = NULL, stel = NULL) 
{
    nG <- nrow(G)
    D <- e2dist(X, G)
    alpha0 <- start[1]
    sigma <- exp(start[2])
    alpha2 <- start[3]
    n0 <- exp(start[4])
    beta <- start[5]
    a0 <- 1
    if (!is.null(zall)) {
        psi <- exp(beta * zall)
        psi <- psi/sum(psi)
    }
    else {
        psi <- rep(1/nG, nG)
    }
    if (!is.null(y)) {
        loglam <- alpha0 - (1/(2 * sigma * sigma)) * D * D + 
            alpha2 * ztrap
        probcap <- 1 - exp(-exp(loglam))
        Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
        ymat <- y
        ymat <- rbind(y, rep(0, ncol(y)))
        lik.marg <- rep(NA, nrow(ymat))
        for (i in 1:nrow(ymat)) {
            Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, 
                nG), probcap[1:length(Pm)], log = TRUE))
            lik.cond <- exp(colSums(Pm))
            lik.marg[i] <- sum(lik.cond * psi)
        }
        nv <- c(rep(1, length(lik.marg) - 1), n0)
        part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
        part2 <- sum(nv * log(lik.marg))
        out <- -1 * (part1 + part2)
    }
    else {
        out <- 0
    }
    if (!is.null(ntel) & !is.null(stel)) {
        D2 <- e2dist(stel, G)^2
        lam <- t(exp(a0 - (1/(2 * sigma * sigma)) * t(D2) + alpha2 * 
            zall))
        denom <- rowSums(lam)
        probs <- lam/denom
        tel.loglik <- -1 * sum(ntel * log(probs))
        out <- out + tel.loglik
    }
    if (!is.null(ntel) & is.null(stel)) {
        D2 <- e2dist(G, G)^2
        lam <- t(exp(a0 - (1/(2 * sigma * sigma)) * t(D2) + alpha2 * 
            zall))
        denom <- rowSums(lam)
        probs <- t(lam/denom)
        temp <- exp(ntel \%*\% log(probs))
        marg <- as.vector(rowSums(temp * psi))
        tel.loglik <- -1 * sum(log(marg))
        out <- out + tel.loglik
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
