\name{hra}
\alias{hra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
home range area
}
\description{
computes implied home range area model 
}
\usage{
hra(func, parms, plot = TRUE, xlim, ylim)
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

# to produce results from Chapter 4:

parms<-matrix(c(
 -2.565099, -2.598714, -1.514959, -0.4701133,
  2.462154,  2.558201,  3.590701,  3.8634103),nrow=2,byrow=TRUE)


pcloglog<-function(parms,D){
lp<-  parms[1] -parms[2]*D*D
p<- 1-exp(-exp(lp))
p
}

pGauss<-function(parms,D){
lp<-  parms[1] -parms[2]*D*D
p<- plogis(parms[1])*exp(-parms[2]*D*D)
p
}

pexp<-function(parms,D){
p<- plogis(parms[1])*exp(-parms[2]*D)
p
}

plogit<-function(parms,D){
p<- plogis(parms[1] - parms[2]*D*D)
p
}

hra(pGauss,parms[,1],xlim=c(-1,7),ylim=c(-1,7))

hra(pcloglog,parms[,2],xlim=c(-1,7),ylim=c(-1,7))

hra(pexp,parms[,3],xlim=c(-1,7),ylim=c(-1,7))

hra(plogit,parms[,4],xlim=c(-1,7),ylim=c(-1,7))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
