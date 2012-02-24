\name{scrUN}
\alias{scrUN}
\title{
  Fit SCR model to data collected on unmarked animals.
}
\description{
  This funciton uses an MCMC algorithm written in pure R code.
}
\usage{
scrUN(y, X, Zknown, M, niters, xlims, ylims, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Zknown}{
%%     ~~Describe \code{Zknown} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{niters}{
%%     ~~Describe \code{niters} here~~
}
  \item{xlims}{
%%     ~~Describe \code{xlims} here~~
}
  \item{ylims}{
%%     ~~Describe \code{ylims} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
library(coda)

# Analysis of the Northern Parula data
data(nopa)
fm1 <- scrUN(nopaDat$y, nopaDat$X, M=100, obsmod="pois", niters=1000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.5,0.2,1))
plot(mcmc(fm1))

fm2 <- scrUN(nopaDat$y, nopaDat$X, M=80, obsmod="bern", niters=12000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.4, 0.05, 15))
mc2 <- mcmc(fm2)
plot(mc2)
rejectionRate(mc2)


# Simulate data using dimensions of parula data

N <- 50
X <- nopaDat$X
J <- nrow(X)
K <- ncol(nopaDat$y)
z <- array(NA, c(N, J, K))
s <- cbind(runif(N, 0, 24), runif(N, 0, 16))
p0 <- 0.3
sigma <- 2
for(j in 1:J) {
    dist <- sqrt((s[,1]-X[j,1])^2 + (s[,2]-X[j,2])^2)
    p <- p0 * exp(-dist^2/(2*sigma^2))
    probs <- p / sum(p)
    for(k in 1:K) {
        z[,j,k] <- rbinom(N, 1, p) # latent encounter histories
    }
}
n <- colSums(z) # observed counts

fm3 <- scrUN(n, X, M=100, obsmod="bern", niters=5000,
    xlims=c(0,24), ylim=c(0,16), tune=c(0.3, 0.05, 10))

mc3 <- mcmc(fm3)
plot(mc3)
rejectionRate(mc3)






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
