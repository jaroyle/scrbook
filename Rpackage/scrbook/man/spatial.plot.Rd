\name{spatial.plot}
\alias{spatial.plot}
\title{
Makes a primative spatial plot using a grid of points shaded by the
attribute of interest. 
}
\description{
An easy way to get a quick and dirty look at whats going on
}
\usage{
spatial.plot(x,y,add=FALSE,cx=1)
}
\arguments{
  \item{x}{
    Coordinates of points to plot 
}
  \item{y}{
    Attribute to provide color for
}
\item{add}{
FALSE to generate a new plot, TRUE to add points to an existing plot
}
}
\value{
 NULL
}
\author{
J. Andy Royle
}
\examples{

spatial.plot <-
function(x,y,add=TRUE,cx=1){
 nc<-as.numeric(cut(y,20))
if(!add) plot(x,pch=" ")
 points(x,pch=20,col=terrain.colors(20)[nc],cex=cx)
image.scale(y,col=terrain.colors(20))

}

}
