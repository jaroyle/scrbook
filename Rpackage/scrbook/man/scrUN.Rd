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
ntraps x noccasions matrix of counts 
}
  \item{X}{
ntraps x 2 matrix of trap coordinates
}
  \item{Zknown}{
%%     ~~Describe \code{Zknown} here~~
}
  \item{M}{
Integer larger than maximum possible population size %%     ~~Describe \code{M} here~~
}
  \item{niters}{
Number of MCMC iterations %%     ~~Describe \code{niters} here~~
}
  \item{xlims}{
xlimits of state-space %%     ~~Describe \code{xlims} here~~
}
  \item{ylims}{
ylimits of state-space %%     ~~Describe \code{ylims} here~~
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
Chandler and Royle (2013) AOAS

SCR book Ch. 18

}
\author{
Richard Chandler 
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
fm1 <- scrUN(nopaDat$y, nopaDat$X, M=100, niters=1000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.5,0.2,1))
plot(mcmc(fm1))


% not working quite right
fm2 <- scrUN(nopaDat$y, nopaDat$X, M=80, obsmod="bern", niters=12000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.4, 0.05, 15))
mc2 <- mcmc(fm2)
plot(mc2)
rejectionRate(mc2)


# Simulate data using dimensions of parula data

set.seed(34123)
N <- 20
X <- nopaDat$X
J <- nrow(X)
K <- 10 # ncol(nopaDat$y)
z <- array(NA, c(N, J, K))
#xlims <- c(0, 24)
#ylims <- c(0, 16)
xlims <- c(2, 22)
ylims <- c(2, 14)
s <- cbind(runif(N, xlims[1], xlims[2]),
           runif(N, ylims[1], ylims[2]))
p0 <- 0.8
sigma <- 1
for(j in 1:J) {
    dist <- sqrt((s[,1]-X[j,1])^2 + (s[,2]-X[j,2])^2)
    p <- p0 * exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        z[,j,k] <- rbinom(N, 1, p) # latent encounter histories
#        z[,j,k] <- rpois(N, p) # latent encounter histories
    }
}
n <- colSums(z) # observed counts = sum z over individuals
table(n)

plot(s, pch=16, col=4, asp=1, xlim=xlims)
rect(xlims[1], ylims[1], xlims[2], ylims[2])
points(X, pch="+")

fm3 <- scrUN(n, X, M=100, niters=5000,
    xlims=xlims, ylim=ylims, tune=c(0.1, 0.1, 5))

mc3 <- mcmc(fm3)
plot(mc3)
rejectionRate(mc3)
summary(window(mc3, start=10001))


# Use M-H to update z
fm4 <- scrUN(n, X, M=40, obsmod="bern", niters=50000,
    xlims=xlims, ylim=ylims, tune=c(0.03, 0.04, 2), zGibbs=FALSE)

mc4 <- mcmc(fm4)
plot(mc4)
rejectionRate(mc4)
summary(window(mc4, start=101001))


# Use Gibbs to update z (Why doesn't this work?????)
fm5 <- scrUN(n, X, M=40, obsmod="bern", niters=150000,
    xlims=xlims, ylim=ylims, tune=c(0.03, 0.05, 0.6), zGibbs=TRUE)

mc5 <- mcmc(fm5)
plot(mc5)
rejectionRate(mc5)
summary(window(mc5, start=3001))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
