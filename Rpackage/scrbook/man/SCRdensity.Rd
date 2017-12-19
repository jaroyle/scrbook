\name{SCRdensity}
\alias{SCRdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a spatial capture-recapture density plot %%  ~~function to do ... ~~
}
\description{
requires MCMC output in the form of histories of activity center
coordinates and data augmentation variables. Basically this makes a 2-d
histogram plot. It bins the points and makes an image plot out of it. 
}
\usage{
SCRdensity(obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, Yu = NULL, scalein = 100, scaleout = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
list with elements Sx Sy and w each of which is "niters x M" matrices
where M is the size of the augmented data set and w are the data
augmentation variables (binary)
}
  \item{nx}{
number of bins in the x direction%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
number of bins in the y direction%%     ~~Describe \code{ny} here~~
}
  \item{Xl}{
Lower x limit of the statespace %%     ~~Describe \code{Xl} here~~
}
  \item{Xu}{
Upper x limit of the statespace%%     ~~Describe \code{Xu} here~~
}
  \item{Yl}{
Lower y limit of the statespace %%     ~~Describe \code{Yl} here~~
}
  \item{Yu}{
Upper y limit of the statespace %%     ~~Describe \code{Yu} here~~
}
  \item{scalein}{
Units of the input coordinate system (of Sx and Sy)%%     ~~Describe \code{scalein} here~~
}
  \item{scaleout}{
Units to scale density by in the output%%     ~~Describe \code{scaleout} here~~
}
}
\details{
scalein and scaleout need some explanation. If inputs coordinate system
is in "km" (i.e., UTM divided by 1000) then input scalein=1. If you want
output in "individuals per 100 km^2" then scaleout = 100. 
}
\value{
A matrix of density values. Also produces a figure. 
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Royle et al., (2014), ch. 5. 
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
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
# this takes 341 seconds on a standard CPU circa 2011
unix.time(out<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,keepz=TRUE))

Sx<-out$sims.list$s[,,1]
Sy<-out$sims.list$s[,,2]
z<- out$sims.list$z
obj<-list(Sx=Sx,Sy=Sy,z=z)
## Input scale 1 unit = 100 km^2
##  express output in wolverines per 1000 km^2 (note: this is different
##  from the map in Fig 5.6 of the SCR book)
SCRdensity(obj,scalein=100,scaleout=1000)



#png("wolv_post_s.png",width=5.25,height=5.25, units="in", res=400)
# Fig. 5.7 (more or less)

par(mfrow=c(1,1),mar=c(2,2,2,6))
X<-wolverine$wtraps[,2:3]
X[,1]<-(X[,1]-min(X[,1]))/10000
X[,2]<-(X[,2]-min(X[,2]))/10000

SCRdensity(obj,whichguy=1)
points(X,pch=20,cex=1,col="white")
points(X[30,],pch=1,cex=2,col="black")
#dev.off()











}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
