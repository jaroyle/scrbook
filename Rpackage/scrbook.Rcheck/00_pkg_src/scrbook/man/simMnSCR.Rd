\name{simMnSCR}
\alias{simMnSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bleen
}
\description{
Simulate some independent multinomial SCR data..... %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simMnSCR(parms, K = 7, ssbuff = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parms}{
a vector of parameter values.........%%     ~~Describe \code{parms} here~~
}
  \item{K}{
number of sample periods
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

set.seed(2013)

parms<-list(N=100,alpha0= -.40, sigma=0.5,alpha2=0)

# simulate some data
data<-simMnSCR(parms,K=7,ssbuff=2)

nind<-nrow(data$Ycat)

# data augmentation
M<-200
Ycat<-rbind(data$Ycat,matrix(nrow(data$X)+1,nrow=(M-nind),ncol=data$K))
Sst<-rbind(data$S,cbind(runif(M-nind,data$xlim[1],data$xlim[2]),
                        runif(M-nind,data$ylim[1],data$ylim[2])))

# starting values
zst<-c(rep(1,160),rep(0,40))
inits <- function(){list (z=zst,sigma=runif(1,.5,1) ,S=Sst) }

# parameters to monitor
parameters <- c("psi","alpha0","alpha1","sigma","N","D")

# bundle the data. Note this reuses "data"
data <- list (X=data$X,K=data$K, trap.space=1,Ycat=Ycat,M=M,ntraps=nrow(data$X),ylim=data$ylim,xlim=data$xlim)

cat("
model {
psi ~ dunif(0,1)
alpha0 ~ dnorm(0,10)
sigma ~dunif(0,10)
alpha1<-  1/(2*sigma*sigma)

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
A <- ((xlim[2]-xlim[1])*trap.space)*((ylim[2]-ylim[1])*trap.space)
D <- N/A
}
",file="model.txt")


# commands to use WinBUGS or OpenBUGS
#library("R2WinBUGS")
#out <- bugs (data, inits, parameters, "model.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE)
#Uncomment and use this bugs statement rather than the one above if using OpenBUGS
#out <- bugs (data, inits, parameters, "model.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, program=c("OpenBUGS"))

library("R2jags")
out <- jags (data, inits, parameters, "model.txt", n.thin=1, n.chains=3, n.burnin=1000, n.iter=2000)


# fit the model by jags. This takes less time.
library("rjags")
out1 <- jags.model("model.txt", data, inits, n.chains=3, n.adapt=500)
out2 <- coda.samples(out1,parameters,n.iter=2000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

