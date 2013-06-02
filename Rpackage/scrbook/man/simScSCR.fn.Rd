\name{simScSCR.fn}
\alias{simScSCR.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bleen
}
\description{
Simulates single catch data. Simulates data on a 10 x 10 grid, unit spacing.
}
\usage{
simScSCR.fn(parms, K = 7, ssbuff = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parms}{
%%     ~~Describe \code{parms} here~~
}
  \item{K}{
number of replicate samples
}
  \item{ssbuff}{
%%     ~~Describe \code{ssbuff} here~~
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

### illustration of fitting independent multinomial model to
### single-catch trap data
##

set.seed(2013)
# simulate some data
parms<-list(N=100,alpha0= -1.50, sigma=0.5,alpha2=0)
data<-simScSCR.fn(parms,K=7,ssbuff=2)

nind<-nrow(data$Ycat)
M<-150

## data augmentation
Ycat<-rbind(data$Ycat,matrix(nrow(data$X)+1,nrow=(M-nind),ncol=data$K))
Sst<-rbind(data$S[1:nind,],cbind(runif(M-nind,data$xlim[1],data$xlim[2]),
                        runif(M-nind,data$ylim[1],data$ylim[2])))
zst<-c(rep(1,M*.6),rep(0,M*.4))
# starting values and bundle data
inits <- function(){list (z=zst,sigma=runif(1,.5,1) ,S=Sst,psi=.5) }
data <- list (X=data$X,K=data$K,
Ycat=Ycat,M=M,ntraps=nrow(data$X),ylim=data$ylim,xlim=data$xlim)

# write out the model

cat("
model {
psi ~ dunif(0,1)
alpha0 ~ dnorm(0,10)
sigma ~dunif(0,10)
alpha1<- 1/(2*sigma*sigma)

for(i in 1:M){
  z[i] ~ dbern(psi)
  S[i,1] ~ dunif(xlim[1],xlim[2])
  S[i,2] ~ dunif(ylim[1],ylim[2])
  for(j in 1:ntraps){
    #distance from capture to the center of the home range
    d[i,j] <- pow(pow(S[i,1]-X[j,1],2) + pow(S[i,2]-X[j,2],2),1)
  }
  for(k in 1:K){
    for(j in 1:ntraps){
      lp[i,k,j] <- exp(alpha0 - alpha1*d[i,j])*z[i]
      cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,]))
    }
    cp[i,k,ntraps+1] <- 1-sum(cp[i,k,1:ntraps])  # last cell = not captured
    Ycat[i,k] ~ dcat(cp[i,k,])
  }
}
N <- sum(z[1:M])
A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
D <- N/A
}
",file="model.txt")


# parameters to monitor
parameters <- c("psi","alpha0","alpha1","sigma","N","D")

#library("R2WinBUGS")
#out <- bugs (data, inits, parameters, "model.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE)
#Uncomment and use this bugs statement rather than the one above if using OpenBUGS
#out <- bugs (data, inits, parameters, "model.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, program=c("OpenBUGS"))

library("rjags")
out1 <- jags.model("model.txt", data, inits, n.chains=3, n.adapt=500)
out2 <- coda.samples(out1,parameters,n.iter=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
