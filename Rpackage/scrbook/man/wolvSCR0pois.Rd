\name{wolvSCR0pois}
\alias{wolvSCR0pois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Poisson encounter model to wolverine data
}
\description{
This is an R script that processes the wolverine data,
produces a BUGS model file, and then executes WinBUGS to obtain the
posterior distribution of model parameters including density (wolverines
per 1000 km^2). SCR model is the Poisson encounter model. 

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0pois(y3d, traps, nb = 1000, ni = 2000, buffer = 2, M = 200)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{y3d}{
3-dimensional encounter history array having dimensions nind x
  ntraps x noccasions.  Can be produced using the encounter data file
  (EDF) format and the function SCR23darray (see example). 
  }
  \item{traps}{
Standard trap deployment file (TDF). %%     ~~Describe \code{traps} here~~
}
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
Returns the output file (a list object) produced by the bugs() function from the
  R2WinBUGS package. 
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book section 9.1 %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov %%  ~~who you are~~
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
data(wolverine)

# convert the data to a 3-d array
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)

## Compare the Bernoulli and Poisson encounter models. This takes an
## hour or two to run
traps<-wolverine$wtraps
toad<-wolvSCR0(y3d,traps,nb=1000,ni=6000, buffer=2, M=200)
toad2<-wolvSCR0pois(y3d,traps,nb=1000,ni=6000, buffer=2, M=200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
