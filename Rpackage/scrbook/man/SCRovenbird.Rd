\name{SCRovenbird}
\alias{SCRovenbird}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fit SCR model to ovenbird data using JAGS%%  ~~function to do ... ~~
}
\description{
This R script will do a Bayesian analysis of the ovenbird data set from
M. Efford's secr package (see below for references).
}
\usage{
SCRovenbird()
}
%- maybe also 'usage' for other objects documented here.
\details{
The approach uses a brute-force "T-fold" data augmentation in which each
of the T years of data has it's own data augmentation parameter,
psi[t]. The formulation here requires a year indicator variable to keep
track of which year or stratum each record belongs to.

The data from each year are stacked into a single large data set and
things are indexed using a fixed "year" index variable which is input as
data.

It is also easy to do this type of T-fold data augmentation using a 3-D
array.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
See the helpfile ?ovenbird from the R package 'secr' 
  
Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture-recapture studies. Biometrics 64, 377-385. 

Dawson, D. K. and Efford, M. G. (2009) Bird population density estimated from acoustic signals. Journal of Applied Ecology 46, 1201-1209. 

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software for analysing capture-recapture data from passive detector arrays. Animal Biodiversity and Conservation 27, 217-228. 

Stamm, D. D., Davis, D. E. and Robbins, C. S. (1960) A method of studying wild bird populations by mist-netting and banding. Bird-Banding 31, 115-130. 
}
\author{
Andy Royle, aroyle@usgs.gov 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("scrbook")
library("secr")
data(ovenbird)
set.seed(2013)


png("Ovenbird_traps.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,3))
 plot(ovenCH[["2005"]])
 plot(ovenCH[["2007"]])
 plot(ovenCH[["2009"]])
dev.off()

## extract the trap locations and create a state-space by adding 150 m
X<-traps<-traps(ovenCH)
xlim<-c(min(X[[1]][,1])-150,max(X[[1]][,1])+150)
ylim<-c(min(X[[1]][,2])-150,max(X[[1]][,2])+150)
ntraps<- nrow(traps[[1]])

## Y are the encounter history data
Y<-ovenCH
K<-10  # number of sample occasions in each year
M<-100 # do constant data augmentation to all years

## starting values for each individual's activity centers
Sst0<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
Sst<-NULL
Ycat<-NULL
died<-NULL

####
####  This look stacks all of the data matrices into one large data matrix.
####
for(i in 1:5){
   tmp<-Y[[i]]
   nind<-nrow(tmp)
   nrep<-ncol(tmp)
 ## make an indicator variable of "died"
   D<-matrix(0,nrow=M,ncol=10)
   D[1:nind,1:nrep]<-tmp
   D[D>0]<- 0
   D[D<0]<- 1
   died<-rbind(died,D)

  tmp[tmp<0]<- 0 ## dead guy set to 0 b/c "died" created above
  tmp[tmp==0]<-ntraps+1

  # In year 1 there were 9 sample occasions so we pad the 10th with NA for that year
  # We indicate non-capture by "J+1" as the convention for specification
  # of the categorical distribution in BUGS
  tmp2<-matrix(NA,nrow=M,ncol=10)
  tmp2[,1:nrep]<-ntraps+1
  tmp2[1:nind,1:nrep]<-tmp

  Ycat<- rbind(Ycat, tmp2)

  sout<-spiderplot(tmp2[1:nind,1:nrep],as.matrix(X[[i]]))$avg.s
  Stmp<-Sst0
  Stmp[1:nind,1:2]<-sout
  Sst<-rbind(Sst,Stmp)
}
 ### This loop will determine whether an individual is known to be dead so that it's encounter probability
 ### is zeroed out
for(i in 1:nrow(died)){
xx<-died[i,]
if(sum(xx)>0){
first<-(1:length(xx))[xx==1]
died[i,first:ncol(died)]<-1
died[i,first]<-0
}
}

##
## This bit of code dumps out the BUGS model file
##
cat("
model {
 # year-specific N parameterized in DA parameter

alpha0 ~ dnorm(0,.1)
sigma ~dunif(0,200)
alpha1<- 1/(2*sigma*sigma)

A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
for(t in 1:5){
N[t] <- inprod(z[1:bigM],yrdummy[,t])
D[t] <- (N[t]/A)*10000  # put in units of per ha
psi[t] ~ dunif(0,1)
}

for(i in 1:bigM){

  z[i] ~ dbern(psi[year[i]])
  S[i,1] ~ dunif(xlim[1],xlim[2])
  S[i,2] ~ dunif(ylim[1],ylim[2])

 for(j in 1:ntraps){
    d2[i,j] <- pow(pow(S[i,1]-X[j,1],2) + pow(S[i,2]-X[j,2],2),1)
  }
for(k in 1:K){
 Ycat[i,k] ~ dcat(cp[i,k,])

  for(j in 1:ntraps){
    ## This is where the died[,] indicator is used to zero-out encounter probability for
    ## dead individuals
    lp[i,k,j] <- exp(alpha0 - alpha1*d2[i,j])*z[i]*(1-died[i,k])
    cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,1:ntraps]))
  }
  cp[i,k,ntraps+1] <- 1-sum(cp[i,k,1:ntraps])  # last cell = not captured

}
}

}
",file="model.txt")
###
###
###

nind<-c(nrow(Y[[1]]),nrow(Y[[2]]),nrow(Y[[3]]),nrow(Y[[4]]),nrow(Y[[5]]))
yrid<-sort(rep(1:5,M))
yrdummy<-as.numeric(c(yrid==1,yrid==2,yrid==3,yrid==4,yrid==5))
yrdummy<-matrix(yrdummy,ncol=5,byrow=FALSE)

bigM<- 5*M
## starting values for data augmentation parameters
zst<-NULL
for(i in 1:5){
zst<-c(zst,rep(1,nind[i]),rep(0,M-nind[i]))
}

inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst,alpha0=runif(1,-2,-1))}
## parameters to monitor
parameters <- c("psi","alpha0","alpha1","sigma","N","D")
## data used in BUGS model
data <- list (died=died,yrdummy=yrdummy,year=yrid,X=as.matrix(X[[1]]),K=10,Ycat=Ycat,bigM=bigM,ntraps=ntraps,ylim=ylim,xlim=xlim)

##
### This takes ~1 hour or so to run
##
library("R2jags")
out <- jags(data, inits, parameters, "model.txt", n.thin=1,n.chains=3,
n.burnin=1000,n.iter=2000,DIC=FALSE)



###
###
### Now fit some models in secr
###
###
###

## fit constant-density model
ovenbird.model.1 <- secr.fit(ovenCH)
## fit net avoidance model
ovenbird.model.1b <- secr.fit(ovenCH, model =   list(g0~b))
## fit model with time trend in detection
ovenbird.model.1T <- secr.fit(ovenCH, model =    list(g0 ~ T))
## fit model with 2-class mixture for g0
ovenbird.model.h2 <- secr.fit(ovenCH, model =    list(g0~h2))
## fit a model with session (year)-specific Density
ovenbird.model.D <- secr.fit(ovenCH, model =    list(D~session))

## compare & average pre-fitted models
AIC (ovenbird.model.1, ovenbird.model.1b, ovenbird.model.1T,
    ovenbird.model.h2,ovenbird.model.D)




}




% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
