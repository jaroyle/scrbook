\name{hra}
\alias{hra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hra: compute homerange area
}
\description{
This functions computes the implied homerange area of an encounter
probability model under the assumption that individuals are making
independent movements about their home range center.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hra(func, parms, plot = TRUE, xlim, ylim, ng = 100, target.area = NULL, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
A user-defined encounter probability function that takes arguments "parms" and
distance (vector-valued) D. Internally, D is the distance vector from a
central point s0 to a raster of points defining the state-space about
the point s0. %%     ~~Describe \code{func} here~~
}
  \item{parms}{
A scalar or vector of parameters used by "func"%%     ~~Describe \code{parms} here~~
}
  \item{plot}{
If TRUE a contour plot of the encounter probability model is produced. %%     ~~Describe \code{plot} here~~
}
  \item{xlim}{
Vector extent (lower and upper bound) of the state-space in the x-direction.%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
Vector extent (lower and upper bound) of the state-space in the y-direction.%%     ~~Describe \code{ylim} here~~
}
\item{ng}{
  Size of one-dimension of the state-space grid used to evaluate the
  encounter probability model. ng=100 produces 10000 points. 
%%     ~~Describe \code{ng} here~~
}
\item{target.area}{
  If specified, the "hra" function finds the value of parms[2] which
  produces the target homerange area. 
%%     ~~Describe \code{target.area} here~~
}
\item{tol}{
  How close do we want to get to the target area? (in absolute
  terms). This function needs revised so that "tol" is expressed
  relative to the parameter..... TO DO......
%%     ~~Describe \code{tol} here~~
}
}
\details{
The algorithm computes the probability of encounter on each pixel in the
ng x ng raster and takes the weighted area over the whole raster as a
statement of implied homerange area.
This is the homerange area under a space usage or resource selection
model that is proportional to the encounter probability model, and in
which animal use decisions are independent multinomial outcomes. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the implied homerange area if target.area=FALSE, otherwise
  returns the value of parms[2] that produces the desired target.area. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book, Section 5.4%% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov%%  ~~who you are~~
}
\note{
    Function assumes parms[2] is the
parameter controling the rate-of-decrease in encounter probability. You
will get funny results if that is not the case. %%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

### Define an encounter probability model. Here we use the Gaussian
###  kernel model

pGauss1<-function(parms,Dmat){
   a0<-parms[1]
   sigma<-parms[2]
   p<-  (exp(parms[1])/(1+exp(parms[1])))*exp( -(1/(2*parms[2]*parms[2]))*Dmat*Dmat )
   p
}

## Set boundary of state-space
xlim<-c(0,6)
ylim<-c(0,6)


# Compute home range area for sigma = 0.3993
hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,tol=.0005)

# This should produce the result 3.00735
#
# Now what is the value of sigma that produces targeta area of 3.00735?

hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,target.area=3.00735,tol=.0005)

##
# Under the Gaussian kernel model, the 95\% homerange radius is
#   r95 = sqrt(5.99)*parms[2] and so the area95 = 2*pi*r95*r95
#   solving area95 = pi*r95*r95 for parms[2] yields
#   sigma95= sqrt(area95/(pi*5.99))
# Therefore, the value of sigma to produce homerange area of 3 is this:

> sqrt(3.00/pi)/sqrt(5.99)

[1] 0.3992751


## Similar results are produced using the Gaussian hazard model:

pGauss2<-function (parms, Dmat) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    lp <- parms[1] - (1/(2 * parms[2] * parms[2])) * Dmat * Dmat
    p <- 1 - exp(-exp(lp))
    p
}


}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
