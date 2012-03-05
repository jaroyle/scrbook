\name{intlik3}
\alias{intlik3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computes marginal likelihood for model SCR0 
  %%  ~~function to do ... ~~
}
\description{
Use this to obtain MLEs for the binomial SCR model -- based on the 2-d
array of individual x trap encounter frequencies
}
\usage{
intlik3(start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs, ssbuffer = 2)
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

library("scrbook")
data("wolverine")

traps<-wolverine$wtraps
traplocs<-traps[,1:2]/10000
K.wolv<-apply(traps[,3:ncol(traps)],1,sum)
traps<-cbind(1:nrow(traps),traps)  # pad an ID variable
y3d<-SCR23darray.fn(wolverine$wcaps,traps)
y2d<-apply(y3d,c(1,3),sum)


starts<-c(-1.5,1.2,log(4))
frog<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)

#
#test out different integration grid densities. This takes a long time
#
#unix.time(a1<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.3,ssbuffer=2)$estimate)
#unix.time(a2<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)$estimate)
#unix.time(a3<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.1,ssbuffer=2)$estimate)
#unix.time(a4<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.05,ssbuffer=2)$estimate)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
