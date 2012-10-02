\name{intlik2}
\alias{intlik2}
\title{
Maximum likelihood estimation of SCR0 
}
\description{
Evaluates the likelihood for model SCR0 for a given set of parameter
values, trap locations and data. 
}
\usage{
intlik2(parm, y = y, delta = 0.3, X = traplocs, ssbuffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parm}{
%%     ~~Describe \code{parm} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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


data<-simSCR0.fn(discard0=TRUE,sd=2013) 
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
J<-nrow(traplocs)
K<-data$K


starts<-c(-2,0,4)
frog<-nlm(intlik2,c(-2.5,2,log(4)),hessian=TRUE,y=y,X=traplocs,delta=.2,ssbuffer=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
