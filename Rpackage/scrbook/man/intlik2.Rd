\name{intlik2}
\alias{intlik2}
\title{
  Computes the marginal (or integrated) likelihood for model SCR0 under
  the Gaussian (half-normal) encounter probability model. 
}
\description{
Evaluates the marginal likelihood for model SCR0 for a given set of parameter
values, trap locations and data. This extends the intlik1 function by
regarding N as a variable also.
}
\usage{
intlik2(parm, y = y, delta = 0.3, X = traplocs, ssbuffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
     Starting values for the optimization procedure. This is a vector of
    length 3: (alpha0, log(alpha1),log(n0)) for the encounter probability model:
    probcap<- plogis(alpha0)*exp(-alpha1*D*D).

    
    Note the 2nd argument (alpha1) is assumed to be the log-transformed
    alpha1 which keeps alpha1 positive during optimization.
    
    alpha1 = (1/(2*sigma*sigma)) where sigma is the usual scale
    parameter of a Gaussian pdf. 

    n0 is the number of uncaptured individuals, i.e., N = n + n0
%%     ~~Describe \code{parm} here~~
}
  \item{y}{
An observed data set which is a matrix of dimension n x J where n =
number of encounter individuals and J = number of traps. Thus, the
encounters are summarized over sampling occasions. 
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
}
\value{
  Retunrs -1* the log-likelihood
}
\references{
  SCR book, Chapter 6
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
intlik1, intlik3, intlik4 and related functions %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Simulate a data set. Note by setting discard0=TRUE a data set with
##  no "all-0" encounter histories is generated. 
data<-simSCR0(discard0=TRUE,rnd=2013) 
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
J<-nrow(traplocs)
K<-data$K

## Specify some starting values and use the R function "nlm" to optimize
##   the likelihood in order to obtain the MLEs of model parameters.
starts<-c(-2,0,4)
frog<-nlm(intlik2,c(-2.5,2,log(4)),hessian=TRUE,y=y,X=traplocs,delta=.2,ssbuffer=2)

frog

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
