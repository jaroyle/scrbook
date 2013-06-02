\name{ch15racdata}
\alias{ch15racdata}
\docType{data}
\title{
Raccoon spatial mark-resight and telemetry data from Chapter 15
}
\description{
Spatial mark-resight dataset of 104 raccoons and telemetry location data for 38 of these individuals (Sollmann et al in prep), collected between October 2007 to January 2008 on South Core Banks, 
a barrier island within Cape Lookout National Seashore, North Carolina. 
}
\usage{data(ch15racdata)}
\format{
list of 11:
racdata$n    : num [1:20, 1:6] 
$X    : num [1:20, 1:2]
$y    : num [1:104, 1:20, 1:6] 0 0 0 0 0 0 0 0 0 0 ...
$M    : num 500
$Eff  : num [1:500, 1:20, 1:6] 1.86 1.86 1.86 1.86 1.86 ...
$telID: chr [1:500] "FALSE" "FALSE" "TRUE" "TRUE" ...
$locs : List of 38
$mi   : num [1:6] 17 6 4 2 2 2
$mall : num [1:6] 
$SSp  : Formal class 'SpatialPolygonsDataFrame' [package "sp"] with 5 slots 
}
\details{
need package sp
}

\source{
XXXX
}
\references{
Sollmann et al in prep
}

