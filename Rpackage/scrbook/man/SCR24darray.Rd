\name{SCR24darray}
\alias{SCR24darray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SCR encounter data file (EDF) into a 4-d array for fitting
models in BUGS/JAGS
}
\description{
Convert SCR encounter data file (EDF) into a 4-d array useful for
fitting models in BUGS/JAGS. The output array is
'nsessions' x 'nind' x 'ntraps' x 'noccassions'. 


}
\usage{
SCR24darray(edf, tdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
Encounter data file (EDF) with assumed order of columns: session,
individual, occasion, trap. Ideall all of these are integers but the
function tries to convert them if they are not. It is not possible to
make this completely general given all types of naming conventions that
people have so why not just make them integers?
}
  \item{tdf}{
Trap deployment file (TDF). The columns are standard: tdf[,1] = trap ID
(ideally an integer 1:nrow(tdf)); tdf[,2] = x coordinate of trap;
tdf[,3] = y coordinate of trap; tdf[, 4] to tdf[, ncol(tdf)] trap
deployment information (0 if trap not operational in that period, 1 if operational).
}
}
\details{
Important notes:

1. Ideally you should make sure that individual ID (edf[,2]) are sequential
integers although the function will try to make that conversion if they
are not. Beware.

2. Ideally you should make sure that trap ID (edf[,4]) corresponds to
integer traps numbered sequentially 1:nrow(tdf). The function tries to
make the conversion if they are not. Beware.

3. The number of distinct sampling occasions is necessarily equal to
(ncol(tdf) - 3). The sampling occasion ID (edf[,3]) should be an integer
that is not greater than the total number of sampling occasions. Note
however if the column names of the TDF are integers then they are
treated as days and there may be gaps. In this case max(edf[,3]) should
not be greater than the maximum column name of the TDF[,4:ncol(TDF)]. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a 4-d array having dimensions: 
'nsessions' x 'nind' x 'ntraps' x 'noccassions'. The result is suitable
for modeling in BUGS/JAGS.  
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
SCR23darray, array2secr
}
\examples{
data(wolverine)

#> names(wolverine)
#[1] "wcaps"  "wtraps" "grid8"  "grid4"  "grid2"  "wsex"  

 test4d<- SCR24darray(wolverine$wcaps, wolverine$wtraps)

# wolverine data has 1 session only so you get this:
 dim(test4d)

#[1]   1  21  37 165

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
