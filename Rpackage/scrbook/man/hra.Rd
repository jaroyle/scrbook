\name{hra}
\alias{hra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HRA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hra(func, parms, plot = TRUE, xlim, ylim, ng = 100, target.area = NULL, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
%%     ~~Describe \code{func} here~~
}
  \item{parms}{
%%     ~~Describe \code{parms} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{ng}{
%%     ~~Describe \code{ng} here~~
}
  \item{target.area}{
%%     ~~Describe \code{target.area} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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

pGauss1<-function(parms,Dmat){
a0<-parms[1]
sigma<-parms[2]
p<-  expit(parms[1])*exp( -(1/(2*parms[2]*parms[2]))*Dmat*Dmat )
p
}

xlim<-c(0,6)
ylim<-c(0,6)

# compute home range area for sigma = 0.3993
hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,tol=.0005)

# This should produce around 3.00735
#
# now what is the value of sigma that produces targeta area of 3.00735?
hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,target.area=3.00735,tol=.0005)

##
# true sigma to produce area of 3
sqrt(3.007/pi)/sqrt(5.99)
}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
