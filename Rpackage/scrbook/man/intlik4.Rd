\name{intlik4}
\alias{intlik4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
evaluates likelihood for SCR0 given a restricted state-space %%  ~~function to do ... ~~
}
\description{
evaluates the likelihood for SCR0 given a restricted state space or, if
none is provided, it puts a buffer around the trap-array and does the
same thing as function intlik3
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlik4(start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs, G = NULL, ssbuffer = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
parameter starting values: logit(p0), log(1/(2*sigma*sigma)), log(n0)
}
  \item{y}{
2-d encounter frequency matrix: individuals x traps 
}
  \item{K}{
vector of trap-specific sample occasions. 
}
  \item{delta}{
resolution of integration grid (spacing)%%     ~~Describe \code{delta} here~~
}
  \item{X}{
trap locations%%     ~~Describe \code{X} here~~
}
  \item{G}{
state-space grid %%     ~~Describe \code{G} here~~
}
  \item{ssbuffer}{
width of buffer around trap array to define state-space. If G is not provided%%     ~~Describe \code{ssbuffer} here~~
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
Andy Royle, aroyle@usgs.gov 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# fake data with a shapefile

data<-simSCR0(discard0=TRUE,rnd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
Xl<-data$xlim[1]
Yl<-data$ylim[1]
Xu<-data$xlim[2]
Yu<-data$ylim[2]

library(maptools)
library(sp)

delta<-.3
ssbuffer<-2
Xl<-min(X[,1]) -ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])- ssbuffer
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))

data("fakeshapefile")
#### replaces this:
#####SSp<-readShapeSpatial('Sim_Polygon.shp')
Pcoord<-SpatialPoints(G)
PinPoly<-over(Pcoord,SSp)
Pin<-as.numeric(!is.na(PinPoly[,1]))
### over() returns NA when the point is not within any polygon, or the ID number of the polygon the point is in; 
### so it works for more complex multi-polygons, too. 
G<-G[Pin==1,]
frog<-nlm(intlik4,c(-2.5,0,3),hessian=TRUE,y=y,K=K,delta=.3,X=traplocs,G=G)


####
#### obtain MLEs for the wolverine camera trapping data
####
####


library("scrbook")
data("wolverine")
traps<-wolverine$wtraps
traplocs<-traps[,2:3]/10000   # change to columns 2:3
K.wolv<-apply(traps[,4:ncol(traps)],1,sum)

y3d<-SCR23darray(wolverine$wcaps,traps)
y2d<-apply(y3d,c(1,2),sum)

G<-wolverine$grid2/10000    ## Use the 2km grid
plot(G,pch=20)
points(traplocs,pch=20,col="red")
starts<-c(-1.5,0,3)
frog<-nlm(intlik4,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,G=G,ssbuffer=2)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
