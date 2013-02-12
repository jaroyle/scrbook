\name{simPoissonSCR}
\alias{simPoissonSCR.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulatesion SCR data with Poisson observation model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simPoissonSCR(discard0 = TRUE, array3d = FALSE, rnd = 2013)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{discard0}{
%%     ~~Describe \code{discard0} here~~
}
  \item{array3d}{
%%     ~~Describe \code{array3d} here~~
}
  \item{rnd}{
random number seed
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





data<-simPoissonSCR(discard0=TRUE,rnd=2013)
y<-data$Y
nind<-nrow(y)
X<-data$traplocs
K<-data$K
J<-nrow(X)
xlim<-data$xlim
ylim<-data$ylim

## Data augmentation stuff
M<-200
y<-rbind(y,matrix(0,nrow=M-nind,ncol=ncol(y)))
z<-c(rep(1,nind),rep(0,M-nind))

cat("
model {
alpha0~dnorm(0,.1)
alpha1~dnorm(0,.1)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(xlim[1],xlim[2])
 s[i,2]~dunif(ylim[1],ylim[2])
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dpois(lam[i,j])
lam[i,j]<- z[i]*K*exp(alpha0)*exp(- alpha1*d[i,j]*d[i,j])
}
}
N<-sum(z[])
D<- N/64
}
",file = "SCR-Poisson.txt")

sst<-X[sample(1:J,M,replace=TRUE),] # starting values for s
for(i in 1:nind){
if(sum(y[i,])==0) next
sst[i,1]<- mean( X[y[i,]>0,1] )
sst[i,2]<- mean( X[y[i,]>0,2] )
}
sst<-sst + runif(nrow(sst)*2,0,1)/8
data <- list (y=y,X=X,K=K,M=M,J=J,xlim=xlim,ylim=ylim)
inits <- function(){
  list (alpha0=rnorm(1,-2,.4),alpha1=runif(1,1,2),s=sst,z=z,psi=.5)
}
library("R2WinBUGS")
parameters <- c("alpha0","alpha1","N","D")
nthin<-1
nc<-3
nb<-1000
ni<-2000
out1 <- bugs (data, inits, parameters, "SCR-Poisson.txt", n.thin=nthin,n.chains=nc, n.burnin=nb,n.iter=ni,working.dir=getwd(),debug=TRUE)

library(rjags)
jm <- jags.model("SCR-Poisson.txt", data=data, inits=inits, n.chains=nc,
n.adapt=nb)
out2 <- coda.samples(jm, parameters, n.iter=ni, thin=nthin)
summary(out2)
plot(out2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
