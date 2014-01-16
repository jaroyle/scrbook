\name{simSCR0}
\alias{simSCR0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate some SCR data under a simple Bernoulli encounter model.
  %%  ~~function to do ... ~~
}
\description{
  Right now this is a simple R script that provides very little
  flexibility. At some point, we'll generalize the function but, until
  then, you should be able to easily modify this script to suit
  your own needs. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simSCR0(N=100,K=20,alpha0=-2.5, sigma = .5, discard0 = TRUE, array3d = FALSE, rnd = 2013)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ Population size }
  \item{K}{ Number of sampling occasions (nights of camera trapping effort)}
  \item{alpha0}
  { Intercept of the Gaussian encounter probability model:
    probcap<-plogis(alpha0)*exp(-alpha1*D*D)
   }
   \item{sigma}{
   The scale parameter of the Gaussian encounter probability model,
   i.e., alpha1 = (1/(2*sigma*sigma)). 
     }
   \item{discard0}{
     If TRUE only the encounter individuals are returned, if FALSE it
     returns the encounter history array or matrix for all N individuals
     whether encountered or not. 
   }
   \item{array3d}{
  If TRUE, returns the full 3-d encounter history array having
  dimensions nind x ntraps x K (n occasions) 
  If FALSE, summarizes over noccasions and returns the total number of
  encounters out of K.
  }
  \item{rnd}{
   Random number seed.  rnd=2013 should reproduce the results in the SCR
   book. 
}
}
\details{
  This function is meant as a learning tool and a template to modify to
  suit your own needs. It is not a general purpose function.
  
  It simulates Bernoulli encounters of N individuals by a 5x5 array of
  camera traps or similar. Trap coordinates are the integers
  (1,2,3,4,5)x(1,2,3,4,5). The state-space is defined by a buffer
  extending 2 units beyond the trap array, for a total size of 8x8
  units.  i.e., the square [-1,7]x[-1,7].

  The encounter probability model is the Gaussian ("half-normal") model in which,  for p[i,j] = encounter probability of individual i at trap j:
  p[i,j] = p0*exp(-alpha1*d[i,j]^2).

  Default values are logit(p0) = -2.5 and alpha1 = 2 which corresponds
  to sigma = .5 where alpha1=(1/(2*sigma*sigma))
  
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a data list suitable for use in SCR0bayes and other functions. It
  has elements:

  Y = nind x ntraps (array3d=FALSE) of encounter frequencies or N x
  ntraps if discard0=FALSE. If array3d=TRUE, returns nind x K x ntraps
  array which must be summarized to nind x ntraps to use SCR0bayes and
  some other functions.

  traplocs = ntraps x 2 matrix of trap coordinates

  xlim=ylim = vectors of lower and upper x- and y-coordinates of the state-space

  N = true value of N

  alpha0 = true value of logit(p0)

  alpha1 = true value of alpha1 = 1/(2*sigma*sigma)

  sigma = true value of sigma

  K = number of sampling occasions
}
\references{
SCR book, Chapter 5. %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
SCR0bayes, SCR0bayesDss%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("scrbook")
## Simulate a data set for a population of N = 100 individuals
data<-simSCR0(N=100, K=20, alpha0= -2.5, sigma = 0.5, discard0=TRUE,rnd=2013)

##out1<-SCR0bayes(data,M=200,engine="jags",ni=2000,nb=1000)
# Fit the model using WinBUGS:

out2<-SCR0bayes(data,M=200,engine="winbugs",ni=2000,nb=1000)

## Inspect the results....

print(out2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
