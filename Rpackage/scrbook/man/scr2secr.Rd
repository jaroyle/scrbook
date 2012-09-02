\name{scr2secr}
\alias{scr2secr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts trap deployment file with operation history to secr format
}
\description{
yea does that%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scr2secr(scrtraps, type = "proximity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrtraps}{
%%     ~~Describe \code{scrtraps} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
library("secr")
library("scrbook")
data("wolverine")
traps<-as.matrix(wolverine$wtraps)   #[,1:3]
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

trapfile2<-scr2secr(scrtraps=traps,type="proximity")

wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)

wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
