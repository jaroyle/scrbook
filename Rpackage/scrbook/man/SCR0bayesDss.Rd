\name{SCR0bayesDss}
\alias{SCR0bayesDss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a basic SCR model using JAGS or WinBUGS, with a discrete state-space.
}
\description{
Fits "model SCR0" -- Bernoulli encounter model with a Gaussian
(half-normal) encounter probability model, using a discrete state-space
represented by a grid of points. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCR0bayesDss(data, ng = 10, engine = "jags", ni = 2000, nb = 1000, M = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
    A list having elements

     Y = encounter history data, 2-dimensional array, nind x ntraps

     traplocs = ntraps x 2 matrix of trap coordinates

     xlim = vector of length 2, lower and upper x-boundary of state-space

     ylim = vector of length 2, lower and upper y-boundary of
     state-space

     K = number of sampling occasions.

     See help file for simSCR0.
    %%     ~~Describe \code{dataobj} here~~
}
  \item{ng}{
Dimension of one side of the state-space grid. e.g., ng=10 will produce
    a 10 x 10 state-space grid having 100 points. %%     ~~Describe \code{ng} here~~
}

  \item{engine}{
    Uses either "jags" or "winbugs" via the rjags or R2WinBUGS
    packages. Obviously at least one of those must be installed.
    %%     ~~Describe \code{engine} here~~
}
  \item{ni}{
Total number of MCMC iterations. %%     ~~Describe \code{ni} here~~
}
  \item{nb}{
Number of burn-in iterations. %%     ~~Describe \code{nb} here~~
}
\item{M}{
    Pertains to data augmentation. Size of the augmented data
    set. %%     ~~Describe \code{M} here~~
}
}
\details{
   This R functions is a script that organizes a data set, creates the
  BUGS model file, and fits the model to the provided data set. It is
  meant as a learning tool and as a template for developing custom SCR
  models.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The object returned is that provided by the bugs function (R2WinBUGS
package) or the coda.samples function (rjags package). See the help
files for those functions.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  SCR book, Chapter 5.
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
 !!!!!!!!!!!!!!WARNING!!!!!!!!!!!! We believe JAGS does not usually work effectively when
  a discrete state-space is used. WinBUGS seems to always work, but
  runtime increases about linearly with the size of the state-space. 
}
  \seealso{
SCR0bayes, simSCR0, SCR0pois
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("scrbook")
# Test with a small run -- this takes a few minutes
data<-simSCR0.fn(N=100, K=20, discard0=TRUE, sd=2013)

# Fit the model using WinBUGS
out2<-SCR0bayesDss(data,ng=8,M=200,engine="winbugs",ni=2000,nb=1000)
summary(as.mcmc.list(out2$out))

### WARNING: JAGS doesn't work sometimes.
out1<-SCR0bayesDss(data,ng=8,M=200,engine="jags",ni=2000,nb=1000)
summary(out1$out)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
