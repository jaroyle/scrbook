\name{wolvSCR0Dssv2}
\alias{wolvSCR0Dssv2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SCR0 with discrete state-space for wolverine data%%  ~~function to do ... ~~
}
\description{
same thing%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0Dssv2(y3d, traps, nb = 1000, ni = 2000, M = 200, Sgrid = wolverine$grid8, engine = "bugs", area = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{
%%     ~~Describe \code{y3d} here~~
}
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
}
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{Sgrid}{
%%     ~~Describe \code{Sgrid} here~~
}
  \item{engine}{
%%     ~~Describe \code{engine} here~~
}
  \item{area}{
%%     ~~Describe \code{area} here~~
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
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
g<-as.matrix(wolverine$grid8)
bln<-wolvSCR0Dss(y3d,traps,nb=1000,ni=2000,M=160,Sgrid=g,engine="bugs",area=8)
# time trial with/without distance matrix calculation in BUGS//JAGS compare result
# run regular grid with mask?
bln2<-wolvSCR0Dss(y3d,traps,nb=1000,ni=2000,M=160,Sgrid=g,engine="jags",area=8)




}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
