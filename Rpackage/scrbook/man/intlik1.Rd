\name{intlik1}
\alias{intlik1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated likelihood for SCR0 with known N %%   ~~ data name/kind ... ~~
%%  ~~function to do ... ~~
}
\description{
Fits a binomial SCR model to a complete data set (including the 0
observations). %% %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlik1(parm, y = y, delta = 0.2, X = traplocs, ssbuffer = 2)
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

data<-simSCR0.fn(discard0=FALSE,sd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
starts<-c(-2,2)
frog<-nlm(intlik1,starts,y=y,delta=.1,X=traplocs,ssbuffer=2,hessian=TRUE)
frog
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
