\name{spcov}
\alias{spcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spcov(B = 1, v = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(B=1, v=20) {
    R <- data.frame(x=rep(seq(0, B, length=v), each=v),
                    y=rep(seq(0, B, length=v), times=v))
    elev <- function(x, y) (x-.5)+(y-.5) # Elevation is a function of x,y
    D<-e2dist1(R,R)
    V<-exp(-D/2)
    Vi<-solve(V)

    cov1<-t(chol(V))\%*\%rnorm(nrow(R))
    image(matrix(cov1,20,20))
    R$elev <- apply(R[,1:2], 1, function(x) elev(x[1], x[2]))
    R$cov1<-cov1-mean(cov1)
    cov1.fn<-function(newpt,cov1,cov1.coords=cbind(R$x,R$y),Vi){
        newpt<-matrix(newpt,ncol=2)
        k<- exp(-e2dist1(newpt,cov1.coords)/2)
        pred<-k\%*\%Vi\%*\%cov1
        as.numeric(pred)
    }
    return(cov1.fn)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
