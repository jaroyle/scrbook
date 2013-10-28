\name{wolvSCR0}
\alias{wolvSCR0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit model SCR0 to Wolverine data %%  ~~function to do ... ~~
}
\description{
  This is an R script to fit model SCR0 to the wolverine (?wolverine)
  data. Some flexibility to modify the state-space buffer and size of
  the augmented data set is provided. This is a learning tool that is
  meant to illustrate basic data handling and fitting a model in WinBUGS
  and could serve as a template for other analyses. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0(y3d, traps, nb = 1000, ni = 2000,buffer = 2, M = 200,keepz=FALSE,coord.scale=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{
3-dimensional encounter history array: nind x ntraps x noccasions 
}
  \item{traps}{
    Standard trap deployment file (TDF) having at least 3 + noccasions columns:
    trapid, xcoord, ycoord, trap_operation_information
%%     ~~Describe \code{traps} here~~
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

  \item{keepz}{
If TRUE the data-augmentation variables are retained. Necessary for
computing certain summaries of density.
    %%     ~~Describe \code{M} here~~
}
  \item{coord.scale}{
Scaling of the coordinate system. Note suitable values of argument
"buffer" depend on this value!
coord.scale=10000 is used for the wolverine data so that a standard unit
is 10 km. 
    %%     ~~Describe \code{M} here~~
}
}
\details{
The model is fitted using the bugs() function in the R2WinBUGS package. 
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the output list from executing "bugs" (package R2WinBUGS).
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book, Section 5.9%% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
data(wolverine), wolvESA, secr_wolverine, wolvSCR0Dss %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
toad<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,coord.scale=10000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
