\name{geesedata}
\alias{geesedata}
\docType{data}
\title{
Canada geese mark-resight data from Chapter 19
}
\description{
Spatial mark-resight dataset of 654 Canada geese (Rutledge 2013), collected in North Carolina between mid July to the end of October 2008. 
}
\usage{data(geesedata)}
\format{
  The format is:

List of 7

 $ n   : num [1:57, 1:21] 0 8 0 0 0 0 0 0 0 0 ...        #unmarked sightings, site by occasion

 $ X   : num [1:57, 1:2] 3988 3997 3993 3995 3993 ...    #matrix of resighting location X and Y coordinates

 $ y   : num [1:654, 1:57, 1:21] 0 0 0 0 0 0 0 0 0 0 ... #array of indidividual resighting data of marked geese

 $ M   : num 7000					 #augmentation limit, unmarked geese

 $ mmax: num 751					 #augmentation limit, marked geese (total number ever marked)

 $ Eff : num [1:57, 1:21] 0 1 0 0 0 0 0 0 0 0 ...	 #site by occasion matrix with survey effort (surveyed=1, not surveyed=0)
 
$ Sex : num [1:654] 1 2 1 2 1 1 1 2 1 1 ...		 #se of resighted marked geese

}

\source{
Rutledge, M.E., 2013. Impacts of Resident Canada Geese in a Suburban Environment. PhD Dissertation, North Carolina State University,
Raleigh, NC, USA.
}

\references{
Rutledge, M.E., 2013. Impacts of Resident Canada Geese in a Suburban Environment. PhD Dissertation, North Carolina State University,
Raleigh, NC, USA.
}

\examples{
See \code{\link{geeseSMR}} for an example analysis of this data set.
}

\seealso{
\code{\link{geeseSMR}}
}