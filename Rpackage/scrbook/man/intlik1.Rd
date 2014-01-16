\name{intlik1}
\alias{intlik1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the integrated likelihood for model SCR0 with known N %%   ~~ data name/kind ... ~~
%%  ~~function to do ... ~~
}
\description{
This is an R function to illustrate the calculations necessary for
evaluating the SCR0 likelihood for the special case where N is
known. This uses the binomial encounter model and assumes the complete
data set is provided (including the 0
observations). %% %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlik1(parm, y = y, delta = 0.2, X = traplocs, ssbuffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
    Starting values for the optimization procedure. This is a vector of
    length 2: (alpha0, alpha1) for the encounter probability model:
    probcap<- plogis(alpha0)*exp(-alpha1*D*D)

    Note: numerical problems may arise if alpha1 is given a negative
    value. We should restrict it to be positive but haven't done so here.
%%     ~~Describe \code{parm} here~~
}
  \item{y}{
    2-dimensional N (population size) x J (traps) matrix of encounter
    frequencies. This is a "complete" data set which should include
    individuals never captured. 
    %%     ~~Describe \code{y} here~~
}
  \item{delta}{
    Spacing of the grid points for the numerical integration grid. 
    %%     ~~Describe \code{delta} here~~
}
  \item{X}{
Trap locations, a J x 2 matrix where J = number of traps. %%     ~~Describe \code{X} here~~
}
  \item{ssbuffer}{
The width of the buffer around the minimum area rectangle that encloses
    the trap locations. %%     ~~Describe \code{ssbuffer} here~~
}
}
\details{
  For each row in the matrix y the binomial likelihood is computed
    conditional on each possible value of s (corresponding to a point in
    the integration grid) and the average over all possible values of s
    is done to produce the marginal likelihood.

    The point of this function is to show that essential calculation.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns -1* the log-likelihood.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book, Chapter 6 %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  intlik2, intlik3, intlik4 and related functions
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# simulate a data set. In this case, discard0=FALSE returns a "complete"
#  data set, i.e., including the all-0 encounter histories. 
data<-simSCR0(discard0=FALSE,rnd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K

# pick some starting values and use the R function "nlm" to obtain the
# MLE for the case where N is known. 
starts<-c(-2,2)
frog<-nlm(intlik1,starts,y=y,delta=.1,X=traplocs,ssbuffer=2,hessian=TRUE)

frog
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
