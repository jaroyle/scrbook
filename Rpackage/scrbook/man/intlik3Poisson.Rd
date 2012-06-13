\name{intlik3Poisson}
\alias{intlik3Poisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate the SCR Poisson-integrated likelihood from Borchers and Efford
2008. 
}
\description{
This computes the likelihood for a given set of parameters and data and
trap locations. The likelihood is the Poisson-integrated likelihood from
Borchers and Efford 2008. 
}
\usage{
intlik3Poisson(start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs, ssbuffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{ssbuffer}{
%%     ~~Describe \code{ssbuffer} here~~
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


### do a simulation study to compare the binomial form of the
### likelihood with the Poisson-integrated likelihood

tmp<-matrix(NA,nrow=100,ncol=7)
for(i in 1:100){
data<-simSCR0.fn(discard0=TRUE,sd=2013+i) 
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
starts<-c(-2,2,log(4))
frog2<-nlm(intlik3,c(-2.5,2,log(2)),hessian=TRUE,y=y,K=K,X=traplocs,delta=.2,ssbuffer=2)
frog2b<-nlm(intlik3Poisson,c(-2.5,2,log(2)),hessian=TRUE,y=y,K=K,X=traplocs,delta=.2,ssbuffer=2)
tmp[i,]<-c(nind,frog2$estimate,frog2b$estimate)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
