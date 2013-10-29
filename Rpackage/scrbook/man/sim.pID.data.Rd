\name{sim.pID.data }
\alias{sim.pID.data }
\title{
Simulation function for spatial mark-resight data
}
\description{
Simulates spatial mark-resight data, with options to include number of marked individuals known or unknown, 
identification rate of marked indidivuals < 100 %, and telemetry location data.
}
\usage{
sim.pID.data(N=N, K=K, sigma=sigma, lam0=lam0, knownID=knownID,X=X,
                     xlims=xlims, ylims=ylims,  obsmod= c("pois", "bern"), nmarked=c("known", "unknown"),rat=1, tel =0, nlocs=0)
}
\arguments{
  \item{N}{
	Number of individuals in the state space
}
  \item{K}{
    Number of sampling occasions
}
  \item{sigma}{
    scale parameter of bivariate normal model of space usage
}
  \item{lam0}{
    Baseline encounter probability or rate, depending on observation model
}
  \item{knownID}{
    Number of marked animals in state space
}
  \item{X}{
    matrix with trap/resighting locations; X and Y coordinates
}
  \item{xlims, ylims}{
    Outer coordinates of the rectangular state-space, in the form of two vectors (x coordinates and y coordinates)
}
  \item{obsmod}{
   Observation mode; either "pois" for Poisson (also called with replacement) or "bern" for Bernoulli (without replacement)
}
  \item{nmarked}{
   is the number of marked individuals known ("known") or unknown ("unknown") 
}
  \item{rat}{
    rate of individual identification of records of marked animals; defaults to 1 (marked individuals can always be identified individually)
}
  \item{tel}{
    Number of animals with telemetry transmitters; can be equal to or below knownID; defaults to 0;
}
  \item{nlocs}{
   number of telemetry locations to be simulated for each telemetered individual; defaults to 0
}


}
\value{
  A list with simulated data
	\item{n}{all records per sampling location and occasion; includes marked and unmarked}
	\item{Y}{array with N individual encounter histories of all individuals}
	\item{Yknown}{array with individual encounter histories of marked individuals }

	\item{Yobs}{array with individual encounter histories; 
		equal to Yknown if number of marked is known, otherwise only contains encounter histories of animals observed at least once}
	\item{YknownR}{array with dimensions equal to Yobs, but with number of records adjusted for imperfect marked identification rate (rat)}
	\item{counter}{number of records of marked individuals that are not identifiable to individual level}
	\item{locs}{list of telemetry locations, order corresponds to telID}
	\item{telID}{ID number (out of 1:N) of animals with telemetry transmitter}
}

\examples{
####example from Royle et al. 2014; Chapt. 19 XXX
set.seed(2501)

#set input values
N=80
lam0=0.5
knownID=40
rat=0.8
sigma=0.5
K=5

#create grid and state space
coords<-seq(0,7, 1)
grid<-expand.grid(coords, coords)
trapmat<-as.matrix(grid)
buff<- 3*sigma
xl<-min(trapmat[,1])-buff
xu<-max(trapmat[,1])+buff
yl<-min(trapmat[,2])-buff
yu<-max(trapmat[,2])+buff
xlims=c(xl, xu)
ylims=c(yl,yu)
area<-(xu-xl)*(yu-yl)

#simulate data
dat<-sim.pID.data(N=N, K=K, sigma=sigma, lam0=lam0, knownID=knownID,
		X=trapmat, xlims=xlims, ylims=ylims,  obsmod= "pois", 
		nmarked="unknown",rat=1, tel =0, nlocs=0)

#create initial values function for scrPID, set M and tuning parameters
inits<-function(){list(S=cbind(runif(M, xlims[1], xlims[2]), 
		runif(M, ylims[1], ylims[2])), lam0=runif(1, 0.4, 0.6), 
		sigma=runif(1, 0.4, 0.6), psi=runif(1, 0.4, 0.6))}
M<-160
delta=c(0.1, 0.01, 2)

#run model, first m=unknown, then m=known
mod<-scrPID(n=dat$n, X=trapmat, y=dat$Yobs, M=M, obsmod = "pois", 
		nmarked="unknown", niters=20000, xlims=xlims, ylims=ylims, 
		inits=inits(), delta=delta ) )
mod2<-scrPID(n=dat$n, X=trapmat, y=dat$Yobs, M=M, obsmod = "pois", 
		nmarked="known", niters=20000, xlims=xlims, ylims=ylims, 
		inits=inits(), delta=delta ) )	
}

\author{
Rahel Sollmann, Richard Chandler
}

\seealso{
\code{\link{scrPID, scrPID.tel, scrPID.um}}
}
