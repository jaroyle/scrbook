\name{nybears}
\alias{nybears}
\docType{data}
\title{
ny bear data from Cat Sun and Angela Fuller%%   ~~ data name/kind ... ~~
}
\description{
the data from the MEE paper  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(nybears)}
\format{
  The format is:
List of 6
 $ y2d        : num [1:50, 1:49] 0 0 0 0 0 0 0 0 0 0 ...
 $ K          : num 10
 $ traplocs   : num [1:49, 1:2] 5 5 5 5 5 5 5 10 10 10 ...
 $ elevation  : num [1:1600, 1] -1.207 -0.829 -0.056 -1.393 -0.894 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1600] "1" "2" "3" "4" ...
  .. ..$ : NULL
 $ trap2raster: int [1:49] 165 365 565 765 965 1165 1365 170 370 570 ...
 $ ssgrid     : int [1:1600, 1:2] 1 2 3 4 5 6 7 8 9 10 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Var1" "Var2"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(nybears)
## maybe str(nybears) ; plot(nybears) ...
}
\keyword{datasets}
