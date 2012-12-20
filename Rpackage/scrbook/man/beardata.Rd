\name{beardata}
\alias{beardata}
\docType{data}
\title{
The bear data
}
\description{
Data from the Fort Drum black bear hair snare study
}
\usage{data(beardata)}
\format{
  The format is:
 chr "beardata"
}
\details{
loads the Fort Drum beardata.
A list with "trapmat", "bearArray" (3 d array), "flat" (2 d matrix, secr format, with sessions coding for sex) and "sex" (sex of the individuals, 1=female, 2=male)
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

###
### Fit Model M0 to the bear data using WinBUGS
###

library(scrbook)
data(beardata)
trapmat<-beardata$trapmat
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=175
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum) # summarize by ind x rep
y[y>1]<- 1                 # toss out multiple encounters b/c
                           #    traditional CR models ignore space



set.seed(2013)               # to obtain the same results each time
library(R2WinBUGS)
data0<-list(y=y,M=M,K=K)
params0<-list('psi','p','N')
zst=c(rep(1,nind),rbinom(M-nind, 1, .5))
inits =  function() {  list(z=zst, psi=runif(1), p=runif(1)) }

cat("
model {

psi~dunif(0, 1)
p~dunif(0,1)

for (i in 1:M){
   z[i]~dbern(psi)
   for(k in 1:K){
     tmp[i,k]<-p*z[i]
     y[i,k]~dbin(tmp[i,k],1)
      }
     }
N<-sum(z[1:M])
}
",file="modelM0.txt")

fit0 = bugs(data0, inits, params0, model.file="modelM0.txt",n.chains=3, 
       n.iter=2000, n.burnin=1000, n.thin=1,debug=TRUE,working.directory=getwd())

###
### Fit model Mh to the bear data using MCMC algorithm written in R
###

data(beardata)
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=500
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum)
y[y>1]<-1
y<-apply(y,1,sum)   # total encounters out of K

set.seed(2013)

out<-modelMh(y,K,nsim=11000)


### more here
###
###
## maybe str(beardata) ; plot(beardata) ...
}
\keyword{datasets}
