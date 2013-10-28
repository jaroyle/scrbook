\name{wolvSCRMb}
\alias{wolvSCRMb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fit local behavioral response (model Mb, trap level) to wolverine data %%  ~~function to do ... ~~
}
\description{
  This is a complete R script to load and organize the wolverine data,
  produce a BUGS model script and fit the behavioral response model
  using WinBUGS. The model is the "local behavioral" response which
  operates at the trap level. 
}
\usage{
wolvSCRMb(nb = 1000, ni = 2000, buffer = 2, M = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nb}{
MCMC iterations to be discarded as burn-in. 
}
  \item{ni}{
Total MCMC iterations. 
}
  \item{buffer}{
    Buffer width of a rectangular state-space enclosing the trap array.
    Standardized units of 10 km (e.g., buffer = 2 is 20 km). 
    %%     ~~Describe \code{buffer} here~~
}
  \item{M}{
 Size of the augmented data set.%%     ~~Describe \code{M} here~~
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
SCR book section 8.4.3 %% ~put references to the literature/web site here ~
}
\author{
  Andy Royle, aroyle@usgs.gov 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# This takes 20-30 hours to run on a very fast desktop.
# To calibrate on your computer set
# nb=100, ni=200 and M=100.

wolv.Mb<-wolvSCRMb(nb=1000,ni=6000,buffer=2,M=200)

mean(wolv.Mb$sims.list$X1new>wolv.Mb$sims.list$X1obs)
mean(wolv.Mb$sims.list$X2new>wolv.Mb$sims.list$X2obs)
mean(wolv.Mb$sims.list$X3new>wolv.Mb$sims.list$X3obs)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
