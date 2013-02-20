pkgname <- "scrbook"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('scrbook')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Muntjac.csv")
### * Muntjac.csv

flush(stderr()); flush(stdout())

### Name: Muntjac.csv
### Title: Samba's distance sampling data
### Aliases: Muntjac.csv
### Keywords: datasets

### ** Examples

data(Muntjac.csv)
## maybe str(Muntjac.csv) ; plot(Muntjac.csv) ...



cleanEx()
nameEx("SCR0bayes")
### * SCR0bayes

flush(stderr()); flush(stdout())

### Name: SCR0bayes
### Title: Fit SCR model using JAGS or WinBUGS
### Aliases: SCR0bayes
### Keywords: ~kwd1 ~kwd2

### ** Examples



library("scrbook")
data<-simSCR0(discard0=TRUE,rnd=2013)
out1<-SCR0bayes(data,M=200,engine="jags",ni=2000,nb=1000)

out2<-SCR0bayes(data,M=200,engine="winbugs",ni=2000,nb=1000)






cleanEx()
nameEx("SCR0bayesDss")
### * SCR0bayesDss

flush(stderr()); flush(stdout())

### Name: SCR0bayesDss
### Title: Fit SCR0 with a discrete state-space
### Aliases: SCR0bayesDss
### Keywords: ~kwd1 ~kwd2

### ** Examples



library("scrbook")
# test with a small run -- this takes a few minutes
data<-simSCR0.fn(discard0=TRUE,sd=2013)
out1<-SCR0bayesDss(data,ng=8,M=200,engine="jags",ni=2000,nb=1000)
summary(out1$out)
out2<-SCR0bayesDss(data,ng=8,M=200,engine="winbugs",ni=2000,nb=1000)
summary(as.mcmc.list(out2$out))





cleanEx()
nameEx("SCR0pois")
### * SCR0pois

flush(stderr()); flush(stdout())

### Name: SCR0pois
### Title: MCMC algorithm for basic spatial capture-recapture model
### Aliases: SCR0pois

### ** Examples

1



cleanEx()
nameEx("SCR0poisSSp")
### * SCR0poisSSp

flush(stderr()); flush(stdout())

### Name: SCR0poisSSp
### Title: MCMC algorithm for basic spatial capture-recapture model with
###   state-space defined by a shape file
### Aliases: SCR0poisSSp

### ** Examples

1



cleanEx()
nameEx("SCR23darray")
### * SCR23darray

flush(stderr()); flush(stdout())

### Name: SCR23darray
### Title: converts SCR flat format to a 3 d array of dim nind x ndays x
###   ntraps
### Aliases: SCR23darray.fn
### Keywords: datasets

### ** Examples

data(SCR23darray)
## maybe str(SCR23darray.fn) ; plot(SCR23darray.fn) ...



cleanEx()
nameEx("SCRdensity")
### * SCRdensity

flush(stderr()); flush(stdout())

### Name: SCRdensity
### Title: Makes a spatial capture-recapture density plot
### Aliases: SCRdensity
### Keywords: ~kwd1 ~kwd2

### ** Examples

library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
# this takes 341 seconds on a standard CPU circa 2011
unix.time(out<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,keepz=TRUE))

Sx<-out$sims.list$s[,,1]
Sy<-out$sims.list$s[,,2]
w<- out$sims.list$z
obj<-list(Sx=Sx,Sy=Sy,z=z)
SCRdensity(obj)



cleanEx()
nameEx("SCRgof")
### * SCRgof

flush(stderr()); flush(stdout())

### Name: SCRgof
### Title: Gof analysis of the point process model
### Aliases: SCRgof
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

library("R2jags")
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
wsex<-wolverine$wsex

toad1<-wolvSCR0ms(y3d,traps,wsex=wsex,nb=1000,ni=2000,delta=2,M=200,model=1)
bleen<- toad1$BUGSoutput$sims.list

traplocs<-wolverine$wtraps[,2:3]
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 


a<-SCRgof(bleen,5,5,traplocs=traplocs,buffer=.4)














cleanEx()
nameEx("SCRovenbird")
### * SCRovenbird

flush(stderr()); flush(stdout())

### Name: SCRovenbird
### Title: fit SCR model to ovenbird data using JAGS
### Aliases: SCRovenbird
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
K<-10  # number of samples in each year
M<-100 # do constant data augmentation to all years

## starting values for each individual's activity centers
Sst0<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
Sst<-NULL
Ymat<-NULL
died<-NULL
# make the data into a 3-d array

for(i in 1:5){
tmp<-Y[[i]]
nind<-nrow(tmp)
nrep<-ncol(tmp)
D<-matrix(0,nrow=M,ncol=10)
D[1:nind,1:nrep]<-tmp
D[D>0]<- 0
D[D<0]<- 1
died<-rbind(died,D)

tmp[tmp<0]<- 0 ## dead guy set to 0 b/c "died" created above
tmp[tmp==0]<-ntraps+1
tmp2<-matrix(NA,nrow=M,ncol=10)  # pad last col with NA for year 1
tmp2[,1:nrep]<-ntraps+1
tmp2[1:nind,1:nrep]<-tmp
Ymat<- rbind(Ymat, tmp2)

sout<-spiderplot(tmp2[1:nind,1:nrep],as.matrix(X[[i]]))$avg.s
Stmp<-Sst0
Stmp[1:nind,1:2]<-sout
Sst<-rbind(Sst,Stmp)
}
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

inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst,alpha0=runif(1,-2,-1) 
,alpha2=-2,alpha3=-2) }              
## parameters to monitor
parameters <- c("psi","alpha0","alpha1","sigma","N","D")
## data used in BUGS model                                                                  
data <- list (died=died,yrdummy=yrdummy,year=yrid,X=as.matrix(X[[1]]),K=10,Ycat=Ymat,bigM=bigM,ntraps=ntraps,ylim=ylim,xlim=xlim)         

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







graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("SCRpossum")
### * SCRpossum

flush(stderr()); flush(stdout())

### Name: SCRpossum
### Title: JAGS analysis of possum data
### Aliases: SCRpossum
### Keywords: ~kwd1 ~kwd2

### ** Examples


# load the data
library("secr")
data(possum)

x<-as.matrix( traps(possumCH))
# we don't use the mask in this analysis
mask<-possummask
mask[,1]<-mask[,1]-min(x[,1])
mask[,2]<-mask[,2]-min(x[,2])
x[,1]<-x[,1]-min(x[,1])
x[,2]<-x[,2]-min(x[,2])
ntraps<-nrow(x)
y<-as.matrix(possumCH)
# categorical trap variable, set to ntraps+1 for "not captured"
# secr uses 0 to indicate "not captured"
y[y==0]<-  ntraps+1




cat("
model {
psi ~ dunif(0,1)
alpha0 ~ dnorm(0,.1)
sigma ~dunif(0,1000)
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
Dha<- D*10000
    
}
",file="model.txt")

nind<-nrow(y)
K<-ncol(y)
M<-300
# describe the state-space 
buff<-100
xlim<-c(min(x[,1])-buff,max(x[,1])+buff)
ylim<-c(min(x[,2])-buff,max(x[,2])+buff)
# data augmentation
yaug<-matrix(ntraps+1,nrow=(M-nind),ncol=K)
Ycat<-rbind(y[1:nind,],yaug[1:(M-nind),])

# this shows a small amount of movement among trapping grids
spider<-spiderplot(y,x)
Sst<-rbind(spider$avg.s,cbind(runif(M-nind,xlim[1],xlim[2]),runif(M-nind,ylim[1],ylim[2])))

zst<-c(rep(1,nind),rep(0,M-nind))
inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst) }              

parameters <- c("psi","alpha0","alpha1","sigma","N","D","Dha")
                                                                   
data <- list (X=x,K=K,Ycat=Ycat,M=M,ntraps=ntraps,ylim=ylim,xlim=xlim)         

# this takes 1-2 hours on a high-end desktop
library("R2jags")
out <- jags(data, inits, parameters, "model.txt", n.thin=1, n.chains=3, 
n.burnin=1000,n.iter=2000,DIC=FALSE)

# print a summary
out






cleanEx()
nameEx("SCRsmy")
### * SCRsmy

flush(stderr()); flush(stdout())

### Name: SCRsmy
### Title: Compute summary statistics of SCR data
### Aliases: SCRsmy
### Keywords: ~kwd1 ~kwd2

### ** Examples


library("scrbook")
data(wolverine)
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
SCRsmy(y3d)




cleanEx()
nameEx("Sim_Polygon")
### * Sim_Polygon

flush(stderr()); flush(stdout())

### Name: Sim_Polygon
### Title: Shapefile for ch5 example
### Aliases: 'fake shapefile'
### Keywords: datasets

### ** Examples

none yet



cleanEx()
nameEx("area")
### * area

flush(stderr()); flush(stdout())

### Name: area
### Title: computes area of a polygon
### Aliases: area
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y) 
{
    if (missing(y)) {
        if (is.matrix(x) && ncol(x) == 2) {
            y <- x[, 2]
            x <- x[, 1]
        }
        else if (!is.null(x$x) && !is.null(x$y)) {
            y <- x$y
            x <- x$x
        }
    }
    x <- c(x, x[1])
    y <- c(y, y[1])
    i <- 2:length(x)
    return(0.5 * sum(x[i] * y[i - 1] - x[i - 1] * y[i]))
  }



cleanEx()
nameEx("array2secr")
### * array2secr

flush(stderr()); flush(stdout())

### Name: array2secr
### Title: Converts detection histories from BUGS format to secr format
### Aliases: array2secr
### Keywords: ~kwd1 ~kwd2

### ** Examples

library("secr")
library("scrbook")
data("beardata")
.... need to fill in with final bear data format
caps<-array2secr(y, detector="proximity") #no sex specification, all individuals get session=1
caps.sex<-array2secr(y, session=sex, detector="proximity") #each individual has sex assigned 




cleanEx()
nameEx("bbsdata")
### * bbsdata

flush(stderr()); flush(stdout())

### Name: bbsdata.rda
### Title: Muntjac line transect survey data from Samba
### Aliases: bbsdata.rda
### Keywords: datasets

### ** Examples

data(bbsdata)
#look at mourning dove counts
y<-bbsdata$counts[,29]  # pick out 1990
notna<-!is.na(y)
y<-y[notna]
#histogram of mourning dove counts in yr. 1990
hist(y)

#produce spatial plot of forest cover for state of PA
library(maps)
habdata<-bbsdata$habitat
map('state',regions="penn",lwd=2)
spatial.plot(habdata[,2:3],habdata[,"dfor"],cx=2)
map('state',regions="penn",lwd=2,add=TRUE



cleanEx()
nameEx("bcharea")
### * bcharea

flush(stderr()); flush(stdout())

### Name: bcharea
### Title: buffer a convex hull and compute area
### Aliases: bcharea
### Keywords: datasets

### ** Examples

library("rgeos")
library("scrbook")
data(beardata)
trapmat<-beardata$trapmat

B.MMDM<- 2.373

aa<-chull(trapmat)
buffered.area<- area(as.matrix(trapmat[aa,]))

bcharea(B.MMDM,traplocs=trapmat)
cat("area of buffered convex hull: ",buffered.area,fill=TRUE)




cleanEx()
nameEx("bear.JAGS")
### * bear.JAGS

flush(stderr()); flush(stdout())

### Name: bear.JAGS
### Title: analysis of the Ft. Drum black bear data using JAGS
### Aliases: bear.JAGS
### Keywords: JAGS covariate modeling

### ** Examples


##runs SCR model with global behavioral trap response on the Ft Drum bear data with 3 chains, 500 adaptive iterations an 2000 iterations
out<-bear.JAGS(model='SCRB', n.chains=3, n.adapt=500, n.iter=2000)

##look at summary output
summary(out)

## returns:
Iterations = 501:2500
Thinning interval = 1
Number of chains = 3
Sample size per chain = 2000

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

           Mean       SD  Naive SE Time-series SE
D        0.1907  0.01863 0.0002405       0.001370
N      578.3960 56.49120 0.7292982       4.155007
alpha0  -2.8001  0.23204 0.0029956       0.013962
alpha2   0.8887  0.23009 0.0029704       0.012137
psi      0.8886  0.08766 0.0011317       0.006407
sigma    1.9985  0.12876 0.0016622       0.006311

2. Quantiles for each variable:

           2.5
D        0.1467   0.1781   0.1952   0.2061   0.2137
N      445.0000 540.0000 592.0000 625.0000 648.0000
alpha0  -3.2728  -2.9526  -2.7945  -2.6427  -2.3526
alpha2   0.4491   0.7264   0.8867   1.0441   1.3512
psi      0.6820   0.8300   0.9101   0.9600   0.9960
sigma    1.7654   1.9101   1.9918   2.0767   2.2753





cleanEx()
nameEx("beardata")
### * beardata

flush(stderr()); flush(stdout())

### Name: beardata
### Title: The bear data
### Aliases: beardata
### Keywords: datasets

### ** Examples


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



cleanEx()
nameEx("ch11secr-jags")
### * ch11secr-jags

flush(stderr()); flush(stdout())

### Name: ch11secr-jags
### Title: Fit IPP using secr and JAGS
### Aliases: ch11secr-jags

### ** Examples


## Not run: 
##D 
##D library(secr)
##D library(rjags)
##D 
##D data(ch11simData)
##D 
##D ch <- ch11simData$ch.secr
##D msk <- ch11simData$spcov.secr
##D 
##D 
##D # SECR analysis
##D 
##D secr1 <- secr.fit(ch, model=D~canht, mask=msk)
##D 
##D region.N(secr1, se.N=TRUE)
##D 
##D 
##D 
##D 
##D 
##D # JAGS analysis
##D 
##D # JAGS model
##D sink("ippDiscrete.txt")
##D cat("
##D model{
##D sigma ~ dunif(0, 20)
##D lam0 ~ dunif(0, 5)
##D beta0 ~ dunif(-10, 10)
##D beta1 ~ dunif(-10, 10)
##D for(j in 1:nPix) {
##D   mu[j] <- exp(beta0 + beta1*CANHT[j])*pixArea
##D   probs[j] <- mu[j]/EN
##D }
##D EN <- sum(mu[])
##D psi <- EN/M
##D for(i in 1:M) {
##D   z[i] ~ dbern(psi)
##D   s[i] ~ dcat(probs[])
##D   x0g[i] <- Sgrid[s[i],1]
##D   y0g[i] <- Sgrid[s[i],2]
##D   for(j in 1:ntraps) {
##D     dist[i,j] <- sqrt((x0g[i]-traps[j,1])^2 +
##D                      (y0g[i]-traps[j,2])^2)
##D     lambda[i,j] <- lam0*exp(-dist[i,j]^2/(2*sigma^2)) * z[i]
##D     y[i,j] ~ dpois(lambda[i,j])
##D     }
##D   }
##D N <- sum(z[])
##D D <- N/1 # 1ha state-space
##D }
##D ", fill=TRUE)
##D sink()
##D 
##D 
##D 
##D modfile <- "ippDiscrete.txt"
##D 
##D jags.data <- with(ch11simData, {
##D     list(y=ch.jags, CANHT=drop(spcov.jags$CANHT),
##D             nPix=nrow(spcov.jags),
##D             M=nrow(ch.jags), ntraps=nrow(traps),
##D             Sgrid=as.matrix(spcov.jags[,1:2]),
##D             pixArea=25,
##D             traps=traps)
##D     })
##D str(jags.data)
##D 
##D init <- function() {
##D     list(sigma=runif(1, 10, 11), lam0=runif(1),
##D          beta0=-8, beta1=1, #rnorm(1, 1),
##D          s=sample.int(jags.data$nPix, jags.data$M, replace=TRUE),
##D          z=rep(1, jags.data$M))
##D }
##D str(init())
##D 
##D pars <- c("sigma", "lam0", "beta0", "beta1", "N", "EN")
##D 
##D # Obtain posterior samples. This takes a few minutes
##D # Compile and adapt
##D set.seed(03453)
##D jm <- jags.model(modfile, jags.data, init, n.chains=2, n.adapt=1000)
##D # MCMC
##D jags1 <- coda.samples(jm, pars, n.iter=6000)
##D 
##D plot(jags1, ask=TRUE)
##D summary(window(jags1, start=1001))
##D 
##D 
##D unlink(modfile)
##D 
##D 
##D 
##D 
##D jags.est <- summary(jags1)
##D jags.r <- cbind(jags.est$stat[,1:2], jags.est$quant[,c(1,5)])
##D jags.r
##D 
##D secr.est <- predict(secr1)
##D secr.r <- cbind(secr.est[2:3,2:5])
##D secr.r <- rbind(beta=as.numeric(coef(secr1)[2,]), secr.r)
##D secr.r <- data.matrix(rbind(region.N(secr1)[,1:4], secr.r))
##D secr.r
##D 
##D jagsVsecr <-
##D data.frame(Par=rep(c("$\lambda_0$", "$\sigma$", "$\beta_1$", "$N$",
##D                    "$\mathbb{E}[N]$"), each=2),
##D            Truth=rep(c(1, 10, 1, 30, 32.3), each=2),
##D            Software = rep(c("\textbf{JAGS}", "\texttt{secr}"), 5),
##D            rbind(jags.r[5,], secr.r[4,],
##D                  jags.r[6,], secr.r[5,],
##D                  jags.r[4,], secr.r[3,],
##D                  jags.r[2,], secr.r[2,],
##D                  jags.r[1,], secr.r[1,]))
##D 
##D 
##D format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE)
##D 
##D write.table(format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE),
##D             file="../../../Ch11-Statespace/R/jagsVsecr2.txt",
##D             row.names=FALSE,
##D             quote=FALSE, sep=" \t& ", eol=" \\\n ")
##D 
##D 
##D 
##D 
## End(Not run)






cleanEx()
nameEx("ch11simData")
### * ch11simData

flush(stderr()); flush(stdout())

### Name: ch11simData
### Title: Simulated data from an inhomogeneous point process model
### Aliases: ch11simData
### Keywords: datasets

### ** Examples

data(ch9simData)

# Here is how the data were generated

# state-space covariate
set.seed(3453)
v <- 21
dat <- spcov(v=v)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","elev")
image(t(matrix(dat$elev, v, v)))

# Multinomial cell probs
set.seed(300225)
N <- 50
alpha <- 2
dat$cp <- exp(alpha*dat$elev) / sum(exp(alpha*dat$elev))
s.tmp <- rmultinom(1, N, dat$cp) # a single realization to be ignored later
image(t(matrix(dat$elev, v, v)))
points(dat[s.tmp>0,c("x","y")])

# Trap locations
xsp <- seq(0.3, 0.7, 0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))


# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps))

nz <- 50 # augmentation
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps))

sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.8   # basal encounter rate
lam <- matrix(NA, N, ntraps)

s <- matrix(NA, N, 3)
colnames(s) <- c("pixID", "x", "y")

set.seed(5588)
for(i in 1:N) {
    s.i <- sample(1:npix, 1, prob=dat$cp)
    sx <- dat[s.i, "x"]
    sy <- dat[s.i, "y"]
    s[i,] <- c(s.i, sx, sy)
    for(j in 1:ntraps) {
        distSq <- (sx-X[j,1])^2 + (sy - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j] <- rpois(1, lam[i,j])
    }
}
yz[1:nrow(y),] <- y # Fill

sum(y)





# Format using secr

# Create a "traps" object
Xs <- data.frame(X)
colnames(Xs) <- c("x","y")
secr.traps <- read.traps(data=Xs, detector="count")

# Create a "capthist" object
secr.caps <- matrix(NA, sum(y), 5)
colnames(secr.caps) <- c("Session", "ID", "Occasion", "X", "Y")
counter <- 0
for(i in 1:nrow(y)) {
    for(j in 1:ncol(y)) {
        y.ij <- y[i,j]
        if(y.ij==0)
            next
        for(v in 1:y.ij) {
            counter <- counter+1
            secr.caps[counter,] <- c(1, i, 1, X[j,1], X[j,2])
        }
    }
}
ch <- make.capthist(secr.caps, secr.traps, fmt="XY")
plot(ch, tol=0.0005) # ouch

# Make mask
msk <- make.mask(secr.traps, buffer=0.325, spacing=.05, nx=v)
summary(msk)
plot(msk)

ssArea <- attr(msk, "area")*nrow(msk)

covariates(msk) <- data.frame(elev=dat$elev[order(dat$y, dat$x)])

ch9simData <- list(ch.secr=ch, ch.jags=yz, spcov.jags=dat, spcov.secr=msk,
                   traps=X)



cleanEx()
nameEx("e2dist")
### * e2dist

flush(stderr()); flush(stdout())

### Name: e2dist
### Title: Compute the distance between trap locations and animal activity
###   centers.
### Aliases: e2dist

### ** Examples

nAnimals <- 10
nTraps <- 5
activity.centers <- cbind(runif(nAnimals), runif(nAnimals))
trap.locs <- cbind(seq(0, 1, length=nTraps), seq(0, 1, length=nTraps))
e2dist(activity.centers, trap.locs)



cleanEx()
nameEx("fakecorridor")
### * fakecorridor

flush(stderr()); flush(stdout())

### Name: fakecorridor
### Title: Fake buffered river corridor for ecological distance stuff
### Aliases: fakecorridor
### Keywords: datasets

### ** Examples

library("sp")
library("rgeos")
library("scrbook")
data("fakecorridor")
buffer<- 0.5
par(mfrow=c(1,1))
aa<-gUnion(l1,l2)
plot(gBuffer(aa,width=buffer),xlim=c(0,10),ylim=c(0,10))
pg<-gBuffer(aa,width=buffer)
pg.coords<- pg@polygons[[1]]@Polygons[[1]]@coords

xg<-seq(0,10,,30)
yg<-seq(10,0,,30)

delta<-mean(diff(xg))
pts<- cbind(sort(rep(xg,30)),rep(yg,30))
points(pts,pch=20)
points(traps$locs,pch=20,cex=2,col="red")

in.pts<-point.in.polygon(pts[,1],pts[,2],pg.coords[,1],pg.coords[,2])
points(pts[in.pts==1,],pch=20,col="red")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fakeshapefile")
### * fakeshapefile

flush(stderr()); flush(stdout())

### Name: fakeshapefile
### Title: this is the objec
### Aliases: fakeshapefile
### Keywords: datasets

### ** Examples

library("scrbook")
data(fakeshapefile)
library(maptools)
library(sp)

data<-simSCR0.fn(discard0=TRUE,sd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
Xl<-data$xlim[1]
Yl<-data$ylim[1]
Xu<-data$xlim[2]
Yu<-data$ylim[2]

delta<-.3
ssbuffer<-2
Xl<-min(X[,1]) -ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])- ssbuffer
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))

data("fakeshapefile")
#### replaces this:
#####SSp<-readShapeSpatial('Sim_Polygon.shp')
Pcoord<-SpatialPoints(G)
PinPoly<-over(Pcoord,SSp)
Pin<-as.numeric(!is.na(PinPoly[,1]))
### over() returns NA when the point is not within any polygon, or the ID number of the polygon the point is in; 
### so it works for more complex multi-polygons, too. 
G<-G[Pin==1,]
frog<-nlm(intlik4,c(-2.5,2,log(4)),hessian=TRUE,y=y,K=K,delta=.3,X=traplocs,G=G)

## maybe str(fakeshapefile) ; plot(fakeshapefile) ...



cleanEx()
nameEx("geeseSMR")
### * geeseSMR

flush(stderr()); flush(stdout())

### Name: geeseSMR
### Title: MCMC for Canada geese mark-resight analysis from Chapter 19
### Aliases: geeseSMR

### ** Examples

#load data and required packages for analysis
data(geesedata)
library(spatstat)
library(maptools)
library(coda)

#set up data and initial values
trapmat<-geesedata$trapmat
 xxx<-convexhull.xy(vertices(dilation.owin(convexhull.xy(trapmat), 3)))
xxy<-as(xxx, "SpatialPolygons")
M=geesedata$M

initS<-function(){list(theta=runif(2, 0.9,1.1), lam0=runif(1,0.28,0.35), w=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
S=runifpoint(M,xxx)   )}

#run 5000 iterations - will take several hours
mod1<-geeseSMR(y=geesedata$n, X=geesedata$X, Zknown=geesedata$y, M=M, EffAr=geesedata$EffAr, niters=5000,inits=initS(),
		Sex=geesedata$Sex, tune=list(sig=c(0.015, 0.02), lam0=0.02, S=2),SSp=xxy)

#lookat summary output
summary(window(mcmc(mod1), start=1001))



cleanEx()
nameEx("geesedata")
### * geesedata

flush(stderr()); flush(stdout())

### Name: geesedata
### Title: Canada geese mark-resight data from Chapter 19
### Aliases: geesedata

### ** Examples

#load data and required packages for analysis
data(geesedata)
library(spatstat)
library(maptools)
library(coda)

#set up data and initial values
trapmat<-geesedata$trapmat
 xxx<-convexhull.xy(vertices(dilation.owin(convexhull.xy(trapmat), 3)))
xxy<-as(xxx, "SpatialPolygons")
M=geesedata$M

initS<-function(){list(theta=runif(2, 0.9,1.1), lam0=runif(1,0.28,0.35), w=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
S=runifpoint(M,xxx)   )}

#run 5000 iterations - will take several hours
mod1<-geeseSMR(y=geesedata$n, X=geesedata$X, Zknown=geesedata$y, M=M, EffAr=geesedata$EffAr, niters=5000,inits=initS(),
		Sex=geesedata$Sex, tune=list(sig=c(0.015, 0.02), lam0=0.02, S=2),SSp=xxy)

#lookat summary output
summary(window(mcmc(mod1), start=1001))



cleanEx()
nameEx("get.traplocs")
### * get.traplocs

flush(stderr()); flush(stdout())

### Name: get.traplocs
### Title: Make a trapping grid by clicking on some points.
### Aliases: get.traplocs

### ** Examples

# make two line segments however you wish them to look
library("scrbook")

xg<-seq(0,10,,40)
yg<-seq(10,0,,40)
pts<- cbind(sort(rep(xg,40)),rep(yg,40))

out<- get.traplocs(ntraps=10,ssgrid=pts)

points(out$loc,pch=20,col="red")





cleanEx()
nameEx("hra")
### * hra

flush(stderr()); flush(stdout())

### Name: hra
### Title: HRA
### Aliases: hra
### Keywords: ~kwd1 ~kwd2

### ** Examples


pGauss1<-function(parms,Dmat){
a0<-parms[1]
sigma<-parms[2]
p<-  expit(parms[1])*exp( -(1/(2*parms[2]*parms[2]))*Dmat*Dmat )
p
}

xlim<-c(0,6)
ylim<-c(0,6)

# compute home range area for sigma = 0.3993
hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,tol=.0005)

# This should produce around 3.00735
#
# now what is the value of sigma that produces targeta area of 3.00735?
hra(pGauss1,c(-2,.3993),plot=FALSE,xlim,ylim,ng=500,target.area=3.00735,tol=.0005)

##
# true sigma to produce area of 3
sqrt(3.007/pi)/sqrt(5.99)



cleanEx()
nameEx("int2d")
### * int2d

flush(stderr()); flush(stdout())

### Name: int2d
### Title: Crude two-dimensional integration over the unit square.
### Aliases: int2d

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(alpha, delta=0.02) {
  z <- seq(-1+delta/2, 1-delta/2, delta)
  len <- length(z)
  cell.area <- delta*delta
  S <- cbind(rep(z, each=len), rep(z, times=len))
  sum(exp(alpha*elev.fn(S)) * cell.area)
  }



cleanEx()
nameEx("intlik1")
### * intlik1

flush(stderr()); flush(stdout())

### Name: intlik1
### Title: Integrated likelihood for SCR0 with known N
### Aliases: intlik1
### Keywords: ~kwd1 ~kwd2

### ** Examples


data<-simSCR0(discard0=FALSE,rnd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
starts<-c(-2,2)
frog<-nlm(intlik1,starts,y=y,delta=.1,X=traplocs,ssbuffer=2,hessian=TRUE)
frog



cleanEx()
nameEx("intlik2")
### * intlik2

flush(stderr()); flush(stdout())

### Name: intlik2
### Title: Maximum likelihood estimation of SCR0
### Aliases: intlik2
### Keywords: ~kwd1 ~kwd2

### ** Examples



data<-simSCR0(discard0=TRUE,sd=2013) 
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
J<-nrow(traplocs)
K<-data$K


starts<-c(-2,0,4)
frog<-nlm(intlik2,c(-2.5,2,log(4)),hessian=TRUE,y=y,X=traplocs,delta=.2,ssbuffer=2)





cleanEx()
nameEx("intlik3")
### * intlik3

flush(stderr()); flush(stdout())

### Name: intlik3
### Title: Computes marginal likelihood for model SCR0
### Aliases: intlik3
### Keywords: ~kwd1 ~kwd2

### ** Examples


library("scrbook")
data("wolverine")

traps<-wolverine$wtraps
traplocs<-traps[,1:2]/10000
K.wolv<-apply(traps[,3:ncol(traps)],1,sum)
traps<-cbind(1:nrow(traps),traps)  # pad an ID variable
y3d<-SCR23darray(wolverine$wcaps,traps)
y2d<-apply(y3d,c(1,3),sum)


starts<-c(-1.5,0,3)
frog<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)

#
#test out different integration grid densities. This takes a long time
#
#unix.time(a1<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.3,ssbuffer=2)$estimate)
#unix.time(a2<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)$estimate)
#unix.time(a3<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.1,ssbuffer=2)$estimate)
#unix.time(a4<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.05,ssbuffer=2)$estimate)





cleanEx()
nameEx("intlik3Poisson")
### * intlik3Poisson

flush(stderr()); flush(stdout())

### Name: intlik3Poisson
### Title: Evaluate the SCR Poisson-integrated likelihood from Borchers and
###   Efford 2008.
### Aliases: intlik3Poisson
### Keywords: ~kwd1 ~kwd2

### ** Examples



### do a simulation study to compare the binomial form of the
### likelihood with the Poisson-integrated likelihood

tmp<-matrix(NA,nrow=100,ncol=7)
for(i in 1:100){
data<-simSCR0.fn(discard0=TRUE,sd=2013+i) 
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
starts<-c(-2,2,log(4))
frog2<-nlm(intlik3,c(-2.5,2,log(2)),hessian=TRUE,y=y,K=K,X=traplocs,delta=.2,ssbuffer=2)
frog2b<-nlm(intlik3Poisson,c(-2.5,2,log(2)),hessian=TRUE,y=y,K=K,X=traplocs,delta=.2,ssbuffer=2)
tmp[i,]<-c(nind,frog2$estimate,frog2b$estimate)
}






cleanEx()
nameEx("intlik3ed")
### * intlik3ed

flush(stderr()); flush(stdout())

### Name: intlik3ed
### Title: Computes the integrated likelihood for binomial SCR model with
###   ecological distance .
### Aliases: intlik3ed

### ** Examples


library("scrbook")
out<-make.EDcovariates()
covariate<-out$covariate.patchy
set.seed(2013)

N<-200
alpha0<- -2
sigma<- .5
K<- 5

alpha1<- 1/(2*sigma*sigma)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
alpha2<-1
cost<- exp(alpha2*covariate)

tr1<-transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

# make up some trap locations
xg<-seq(1,4,1)
yg<-4:1
pts<-cbind( sort(rep(xg,4)),rep(yg,4))

traplocs<-pts
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))
D<-costDistance(tr1CorrC,S,traplocs)
probcap<-plogis(alpha0)*exp(-alpha1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}
Y<-Y[apply(Y,1,sum)>0,]

n0<- N-nrow(Y)
frog1<-nlm(intlik3ed,c(alpha0,alpha1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="euclid",covariate=covariate,alpha2=1)

frog2<-nlm(intlik3ed,c(alpha0,alpha1,log(n0),-.3),hessian=TRUE,y=Y,K=K,
               X=traplocs,distmet="ecol",covariate=covariate,alpha2=NA)




cleanEx()
nameEx("intlik3edv2")
### * intlik3edv2

flush(stderr()); flush(stdout())

### Name: intlik3edv2
### Title: Integrated likelihood of the binomial SCR model
### Aliases: intlik3edv2

### ** Examples


library("sp")
library("rgeos")
library("scrbook")
data("fakecorridor")

## Step 1
## produce the corridor system:
buffer<- 0.5
par(mfrow=c(1,1))
aa<-gUnion(l1,l2)
plot(gBuffer(aa,width=buffer),xlim=c(0,10),ylim=c(0,10))
pg<-gBuffer(aa,width=buffer)
pg.coords<- pg@polygons[[1]]@Polygons[[1]]@coords
 # note: can you believe this shit?

xg<-seq(0,10,,40)
yg<-seq(10,0,,40)
delta<-mean(diff(xg))
pts<- cbind(sort(rep(xg,40)),rep(yg,40))
points(pts,pch=20,cex=.5)
in.pts<-point.in.polygon(pts[,1],pts[,2],pg.coords[,1],pg.coords[,2])
points(pts[in.pts==1,],pch=20,col="red")

### Step 2
### Assign cost matrix
cost<-rep(NA,nrow(pts))
cost[in.pts==1]<-1   # low cost to move among pixels but not 0
cost[in.pts!=1]<-10000   # high cost 

## Stuff this into a raster
library("raster")
r<-raster(nrows=40,ncols=40)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(0-delta/2,10+delta/2,0-delta/2,10+delta/2)
values(r)<-matrix(cost,40,40,byrow=FALSE)
par(mfrow=c(1,1))
plot(r)
points(pts,pch=20,cex=.4)


library("gdistance")
## use max = doesn't count moving through boundary pixel
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
costs1<-costDistance(tr1CorrC,pts)
outD<-as.matrix(costs1)

plot(pts,pch=".")
points(pts[in.pts==1,],pch=20,col="red")

#
# Note "traps" file in scrbook library based on 40 x 40 raster
#
traplocs<-traps$loc
trap.id<-traps$locid
ntraps<-nrow(traplocs)


set.seed(2013)
###
### Step 3: Simulate some SCR data
###
N<-200
# restrict to points in the corridor system
S.possible<- (1:nrow(pts))[in.pts==1]
S.id<-sample(S.possible,N,replace=TRUE)
S<- pts[S.id,]

D<- outD[S.id,trap.id]
eD<- e2dist(S,traplocs)
# the follow two distance matrices will be used in the likelihood
Dtraps<-outD[trap.id,]
Deuclid<-e2dist(pts[trap.id,],pts)

alpha0<- -1.5
sigma<- 1.5
alpha1<- 1/(2*sigma*sigma)
K<-10

probcap<-plogis(alpha0)*exp(-alpha1*D*D)
probcapE<-plogis(alpha0)*exp(-alpha1*eD*eD)

Y<-matrix(NA,nrow=N,ncol=ntraps)
Ye<-Y
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
 Ye[i,]<-rbinom(ntraps,K,probcapE[i,])
}

Y<-Y[apply(Y,1,sum)>0,]

###
###  Step 4: Run SCR model
###
# first with distance-in-corridor
frog1<-nlm(intlik3edv2,c(-2.5,2,log(4)),hessian=TRUE,y=Y,K=K,X=traplocs,S=pts,D=Dtraps,inpoly=in.pts)
# now run Euclidean distance
frog2<-nlm(intlik3edv2,c(-2.5,2,log(4)),hessian=TRUE,y=Y,K=K,X=traplocs,S=pts,D=Deuclid,inpoly=in.pts)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("intlik4")
### * intlik4

flush(stderr()); flush(stdout())

### Name: intlik4
### Title: evaluates likelihood for SCR0 given a restricted state-space
### Aliases: intlik4
### Keywords: ~kwd1 ~kwd2

### ** Examples


# fake data with a shapefile

data<-simSCR0(discard0=TRUE,rnd=2013)
y<-data$Y
traplocs<-data$traplocs
nind<-nrow(y)
X<-data$traplocs
J<-nrow(X)
K<-data$K
Xl<-data$xlim[1]
Yl<-data$ylim[1]
Xu<-data$xlim[2]
Yu<-data$ylim[2]

library(maptools)
library(sp)

delta<-.3
ssbuffer<-2
Xl<-min(X[,1]) -ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])- ssbuffer
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))

data("fakeshapefile")
#### replaces this:
#####SSp<-readShapeSpatial('Sim_Polygon.shp')
Pcoord<-SpatialPoints(G)
PinPoly<-over(Pcoord,SSp)
Pin<-as.numeric(!is.na(PinPoly[,1]))
### over() returns NA when the point is not within any polygon, or the ID number of the polygon the point is in; 
### so it works for more complex multi-polygons, too. 
G<-G[Pin==1,]
frog<-nlm(intlik4,c(-2.5,0,3),hessian=TRUE,y=y,K=K,delta=.3,X=traplocs,G=G)


####
#### obtain MLEs for the wolverine camera trapping data
####
####

library("scrbook")
data("wolverine")
traps<-wolverine$wtraps
traplocs<-traps[,1:2]/10000
K.wolv<-apply(traps[,3:ncol(traps)],1,sum)
traps<-cbind(1:nrow(traps),traps)  # pad an ID variable
y3d<-SCR23darray.fn(wolverine$wcaps,traps)
y2d<-apply(y3d,c(1,3),sum)
G<-wolverine$grid2/10000    ## Use the 2km grid
starts<-c(-1.5,0,3)
frog<-nlm(intlik4,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,G=G,ssbuffer=2)





cleanEx()
nameEx("jaguarDataCh9")
### * jaguarDataCh9

flush(stderr()); flush(stdout())

### Name: jaguarDataCh9
### Title: Jaguar data used in Chapter 9
### Aliases: jaguarDataCh9
### Keywords: datasets

### ** Examples

data(jaguarDataCh9)



cleanEx()
nameEx("make.EDcovariates")
### * make.EDcovariates

flush(stderr()); flush(stdout())

### Name: make.EDcovariates
### Title: generate covariates for the ecological distance simulations
### Aliases: make.EDcovariates

### ** Examples

library("scrbook")
out<-make.EDcovariates()
par(mfrow=c(1,2))
plot(out$covariate.patchy)
plot(out$covariate.trend)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("make.seg")
### * make.seg

flush(stderr()); flush(stdout())

### Name: make.seg
### Title: Make a line segment by laying down a few points.
### Aliases: make.seg

### ** Examples

# make two line segments however you wish them to look
library("scrbook")
plot(NULL,xlim=c(0,10),ylim=c(0,10))
l1<-make.seg(9)
plot(l1)
l2<-make.seg(5)
plot(l1)
lines(l2)



cleanEx()
nameEx("make.statespace")
### * make.statespace

flush(stderr()); flush(stdout())

### Name: make.statespace
### Title: Create a grid of points over some prescribed region
### Aliases: make.statespace
### Keywords: ~kwd1 ~kwd2

### ** Examples

# This example is from Chapter XXXX on modeling space usage
set.seed(1234)
gr<-make.statespace(minx=1,maxx=40,miny=1,maxy=40,nx=40,ny=40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))
spatial.plot(gr,z)



cleanEx()
nameEx("mallard")
### * mallard

flush(stderr()); flush(stdout())

### Name: mallard
### Title: mallard data
### Aliases: mallard
### Keywords: datasets

### ** Examples

library(scrbook)
data(mallard)
library("R2WinBUGS")


sink("model.txt")
cat("
model {
 for(t in 1:5){
    for (i in 1:nobs){
       y[i,t] ~ dbin(p[i,t], B[i,t])
       logit(p[i,t]) <- beta0[t] + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,1]*X[i,2]
     }
}
	beta1~dnorm(0,.001)
	beta2~dnorm(0,.001)
	beta3~dnorm(0,.001)
	for(t in 1:5){
 	beta0[t] ~ dnorm(0,.001)  
 }
}
",fill=TRUE)
sink()

data <- list(B=mallard$bandings, y=mallard$recoveries,
X=mallard$locs,nobs=nrow(mallard$locs))
inits <- 	function(){
list(beta0=rnorm(5),beta1=0,beta2=0,beta3=0)
}
parms <- list('beta0','beta1','beta2','beta3')
out <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=FALSE)
print(out,digits=3)

### now run to get estimates of p at each grid cell
### and make a nice image plot

parms <- list('beta0','beta1','beta2','beta3','p')
out2 <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=FALSE)

pbar<-apply(out2$sims.list$p,2,mean)
ux<-unique(mallard$locs[,1])
uy<-unique(mallard$locs[,2])

X<- mallard$locs
odx<-order(X[,1],X[,2])
I<-matrix(pbar[odx],nrow=length(uy),ncol=length(ux))
par(mfrow=c(2,1))

image(ux,uy,t(I),col=terrain.colors(10))
spatial.plot(X,pbar,col="notgray",cx=4)




## maybe str(mallard) ; plot(mallard) ...



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("model.average2")
### * model.average2

flush(stderr()); flush(stdout())

### Name: model.average2
### Title: model averaging of secr fit objects
### Aliases: model.average2
### Keywords: ~kwd1 ~kwd2

### ** Examples


to appear






cleanEx()
nameEx("modelMh")
### * modelMh

flush(stderr()); flush(stdout())

### Name: modelMh
### Title: fits logit-normal model Mh
### Aliases: modelMh

### ** Examples

library("scrbook")
data("beardata")
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=500
nz<-M-nind
Yaug <- array(0, dim=c(M,ntraps,K))

Yaug[1:nind,,]<-beardata$bearArray
y<- apply(Yaug,c(1,3),sum)
y[y>1]<-1
ytot<-apply(y,1,sum)   # total encounters out of K

set.seed(2013)
out<-modelMh(ytot,K,nsim=11000)



cleanEx()
nameEx("modelMhBUGS")
### * modelMhBUGS

flush(stderr()); flush(stdout())

### Name: modelMhBUGS
### Title: Fits Model Mh using WinBUGS or JAGS
### Aliases: modelMhBUGS
### Keywords: heterogeneity model Mh

### ** Examples

library("scrbook")
data(beardata)

set.seed(2013)
jout<-modelMhBUGS(beardata$bearArray,engine="jags",ni=2000,nb=1000)

set.seed(2013)
wbout<-modelMhBUGS(beardata$bearArray,engine="winbugs",ni=2000,nb=1000)
summary(as.mcmc.list(wbout))



cleanEx()
nameEx("muntjac.rda")
### * muntjac.rda

flush(stderr()); flush(stdout())

### Name: muntjac.rda
### Title: Muntjac line transect survey data from Samba
### Aliases: muntjac.rda
### Keywords: datasets

### ** Examples

data(muntjac.rda)
## maybe str(muntjac.rda) ; plot(muntjac.rda) ...



cleanEx()
nameEx("nopa")
### * nopa

flush(stderr()); flush(stdout())

### Name: nopa
### Title: The Northern Parula data.
### Aliases: nopa
### Keywords: datasets

### ** Examples

data(nopa)
## maybe str(nopa) ; plot(nopa) ...



cleanEx()
nameEx("pGauss1")
### * pGauss1

flush(stderr()); flush(stdout())

### Name: pGauss1
### Title: pGauss1
### Aliases: pGauss1
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (parms, Dmat) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    p <- plogis(parms[1]) * exp(-(1/(2 * parms[2] * parms[2])) * 
        Dmat * Dmat)
    p
  }



cleanEx()
nameEx("pGauss2")
### * pGauss2

flush(stderr()); flush(stdout())

### Name: pGauss2
### Title: pFauss2
### Aliases: pGauss2
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (parms, Dmat) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    lp <- parms[1] - (1/(2 * parms[2] * parms[2])) * Dmat * Dmat
    p <- 1 - exp(-exp(lp))
    p
  }



cleanEx()
nameEx("pIDgeese")
### * pIDgeese

flush(stderr()); flush(stdout())

### Name: pIDgeese
### Title: MCMC algorithm for analyzing Canada geese data set (Ruthledge
###   2012) from Chapter 15
### Aliases: pIDgeese

### ** Examples

### run 5000 iterations of model for Canada geese
###load data
data('geesedata')  #loads a list gd with all data needed

# extract M and trap location matrix 
M=gd$M
trapmat=gd$X

#write function to generate initial values

initS<-function(){list(theta=runif(2, 0.7,1.2), lam0=runif(1,0.3,0.8), w=rbinom(M,1,0.8), phi=runif(1,0.3,0.7), psi=runif(1, 0.3, 0.8), 
S=cbind( runif( M, min(trapmat[,1]),max(trapmat[,1] ) ), runif( M,min(trapmat[,2]),max(trapmat[,2]) ) ))}

# run model
mod<-pIDgeese (y=gd$y, X=gd$X, Zknown=gd$Zknown, M=M, Eff=gd$EffAr, niters=5000, xl=gd$xl,xu=gd$xu,yl=gd$yl,
yu=gd$yu, inits=initS(), Sex=gd$Sex, delta=list(sig=c(0.015, 0.02), lam0=0.02, S=2))

#look at output
library(coda)
out<-mcmc(mod)
summary(out)



cleanEx()
nameEx("raccoon.MCMC")
### * raccoon.MCMC

flush(stderr()); flush(stdout())

### Name: raccoon.MCMC
### Title: MCMC algorithm to analyze raccoon data set (ch15racdata)
### Aliases: raccoon.MCMC

### ** Examples

#use raccoon.MCMC to run a single MCMC chain for raccoon data from Ch15

data(ch15racdata)
library(spatstat)
library(mvtnorm)
library(coda)

#create spatial window to draw random points within state space
newSSp <- as(racdata$SSp, "SpatialPolygons")
nSSp<-as.owin(newSSp)

#set up inits function
inits=function(){list(theta=runif(1,0.2,0.8), lam0=c( rep(runif(1, .05, .1), 2), 
			    rep(runif(1, .05, .1), 2), rep(runif(1, .05, .1), 2)  ),
                      xz=runifpoint(racdata$M, win=nSSp), w=rbinom(racdata$M,1,.5),
			    psi=runif(1,.2,.8) )}  

#run raccoon model
out<-raccoon.MCMC(n=racdata$n,y=racdata$y, X=racdata$X, mi=racdata$mi, 
                  mall=racdata$mall, M=racdata$M,Eff=racdata$Eff, 
                  collar=racdata$telID, locs=racdata$locs, niters=200,
                  inits=inits(), SSp=racdata$SSp)

#look at output, remove XX iterations as burn-in
mod<-mcmc(out)
summary(window(out, start=10001))



cleanEx()
nameEx("rot")
### * rot

flush(stderr()); flush(stdout())

### Name: rot
### Title: rotates a matrix so image() will plot it as you look at it on
###   the page
### Aliases: rot
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m) 
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
  }



cleanEx()
nameEx("scr2secr")
### * scr2secr

flush(stderr()); flush(stdout())

### Name: scr2secr
### Title: Converts trap deployment file with operation history to secr
###   format
### Aliases: scr2secr
### Keywords: ~kwd1 ~kwd2

### ** Examples

library("secr")
library("scrbook")
data("wolverine")
traps<-as.matrix(wolverine$wtraps)   #[,1:3]
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

trapfile2<-scr2secr(scrtraps=traps,type="proximity")

wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)

wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)




cleanEx()
nameEx("scrDED")
### * scrDED

flush(stderr()); flush(stdout())

### Name: scrDED
### Title: Fit SCR model with state-space covariates and covariates of
###   ecological distance.
### Aliases: scrDED

### ** Examples


library(raster)
library(gdistance)

set.seed(353)
pix <- 0.05
dat <- spcov(pix=pix)$R
npix <- nrow(dat)
colnames(dat) <- c("x","y","elev")
cell <- seq(pix/2, 1-pix/2, pix)
image(cell, cell, t(matrix(dat$elev, 1/pix, 1/pix)), ann=FALSE)

head(dat)

# Simulate IPP
set.seed(30275)
N <- 50
alpha <- 2
dat$cp <- exp(alpha*dat$elev) / sum(exp(alpha*dat$elev))
s.tmp <- rmultinom(1, N, dat$cp) # a single realization to be ignored later

# Trap locations
xsp <- seq(0.275, 0.725, by=0.05)
X <- cbind(rep(xsp, each=length(xsp)), rep(xsp, times=length(xsp)))
str(X)

# Elevation covariate as a matrix, then as a raster
elevMat <- t(matrix(dat$elev, 1/pix, 1/pix))
elev <- flip(raster(t(elevMat)), direction="y")
layerNames(elev) <- "elev"

# Simulate capture histories, and augment the data
npix <- nrow(dat)
ntraps <- nrow(X)
T <- 5
y <- array(NA, c(N, ntraps))

# Parameters
sigma <- 0.1  # half-normal scale parameter
lam0 <- 0.8   # basal encounter rate
lam <- matrix(NA, N, ntraps)
theta <- 1

# Activity centers
s <- matrix(NA, N, 3)
colnames(s) <- c("pixID", "x", "y")

set.seed(557828)
elev.tran <- elev-cellStats(elev, min)
elev.tran <- elev.tran/cellStats(elev.tran, max)
cost <- exp(theta*elev.tran)
tr1 <- transition(cost, transitionFunction = function(x) 1/mean(x),
                  directions=8)
tr1CorrC <- geoCorrection(tr1, type="c", multpl=FALSE, scl=FALSE)
for(i in 1:N) {
    s.i <- sample(1:npix, 1, prob=dat$cp)
    sx <- dat[s.i, "x"]
    sy <- dat[s.i, "y"]
    s[i,] <- c(s.i, sx, sy)
    distSq <- costDistance(tr1CorrC, X, s[i,2:3])^2
    lam[i,] <- exp(-distSq/(2*sigma*sigma)) * lam0
    y[i,] <- rpois(nrow(X), lam[i,])
}

sum(y)



y.ded <- y[rowSums(y)>0,]
str(y.ded)


## Not run: 
##D 
##D (fm1 <- scrDED(y.ded, X, ~1, ~1, rasters=elev,
##D #               start=c(-1, -1, 1),
##D                method="BFGS",
##D                control=list(trace=TRUE, REPORT=1, maxit=50)))
##D 
##D exp(fm1$par[1:2])            # 0.8, 0.1
##D exp(fm1$par[3])+nrow(y.ded)  # 50
##D 
##D 
##D (fm2 <- scrDED(y.ded, X, ~elev, ~1, rasters=elev,
##D #               start=c(log(0.8), log(0.1), log(10), 2),
##D                method="BFGS",
##D                control=list(trace=TRUE, REPORT=1, maxit=500)))
##D 
##D exp(fm2$par[1:2])               # 0.8, 0.1
##D exp(fm2$par[3])+nrow(y.ded)     # 50
##D 
##D 
##D 
##D (fm3 <- scrDED(y.ded, X, ~elev, ~elev, rasters=elev,
##D #               start=c(log(0.8), log(0.1), log(10), 2, 0),
##D                method="BFGS",
##D                control=list(trace=TRUE, REPORT=1, maxit=500)))
##D 
##D exp(fm3$par[1:2])                       # 0.8, 0.1
##D c(N=exp(fm3$par[3])+nrow(y.ded))        # 50
##D fm3$par[4:5]                            # 2, 1
##D 
##D 
## End(Not run)




cleanEx()
nameEx("scrIPP")
### * scrIPP

flush(stderr()); flush(stdout())

### Name: scrIPP
### Title: A Metropolis-within-Gibbs sampler to obtain posterior
###   distributions for parameters of an SCR model with an inhomogeneous
###   binomial point process.
### Aliases: scrIPP

### ** Examples




# Spatial covariate (with mean 0)
elev.fn <- function(s) {
    s <- matrix(s, ncol=2)        # Force s to be a matrix
    (s[,1] + s[,2] - 100) / 40.8  # Returns (standardized) "elevation"
}

mu <- function(s, beta0, beta1) exp(beta0 + beta1*elev.fn(s=s))

library(R2Cuba)
xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(100,100), beta0=0, beta1=2)

xx <- cuhre(2, 1, mu, lower=c(0,0), upper=c(1,1), beta0=0, beta1=2,
            flags=list("verbose"=0))


# Simulate PP using rejection sampling
set.seed(31025)
beta0 <- -6 # intercept of intensity function
beta1 <- 1  # effect of elevation on intensity
# Next line computes integral, which is expected value of N
EN <- cuhre(2, 1, mu, beta0=beta0, beta1=beta1,
            lower=c(0,0), upper=c(100,100))$value
EN
N <- rpois(1, EN) # Realized N
s <- matrix(NA, N, 2) # This matrix will hold the coordinates
elev.min <- elev.fn(c(0,0))
elev.max <- elev.fn(c(100, 100))
Q <- max(c(exp(beta0 + beta1*elev.min),
           exp(beta0 + beta1*elev.max)))
counter <- 1
while(counter <= N) {
  x.c <- runif(1, 0, 100); y.c <- runif(1, 0, 100)
  s.cand <- c(x.c,y.c)
  pr <- mu(s.cand, beta0, beta1) #/ EN
  if(runif(1) < pr/Q) {
    s[counter,] <- s.cand
    counter <- counter+1
    }
  }

plot(s)


xsp <- seq(20, 80, by=10); len <- length(xsp)
X <- cbind(rep(xsp, each=len), rep(xsp, times=len)) # traps
ntraps <- nrow(X); noccasions <- 5
y <- array(NA, c(N, ntraps, noccasions)) # capture data
sigma <- 5  # scale parameter
lam0 <- 1   # basal encounter rate
lam <- matrix(NA, N, ntraps)
set.seed(5588)
for(i in 1:N) {
    for(j in 1:ntraps) {
        # The object "s" was simulated in previous section
        distSq <- (s[i,1]-X[j,1])^2 + (s[i,2] - X[j,2])^2
        lam[i,j] <- exp(-distSq/(2*sigma^2)) * lam0
        y[i,j,] <- rpois(noccasions, lam[i,j])
    }
}
# data augmentation
nz <- 80
M <- nz+nrow(y)
yz <- array(0, c(M, ntraps, noccasions))
yz[1:nrow(y),,] <- y # Fill data augmentation array





# Fit the model using MCMC
# Sample the parameters: "sigma", "lam0", "beta0", "beta1", "N", "EN"

set.seed(3434)
system.time({
fm1 <- scrIPP(yz, X, M, 10000, xlims=c(0,100), ylims=c(0,100),
              space.cov=elev.fn,
              tune=c(0.4, 0.2, 0.3, 0.3, 7))
}) # 328s



library(coda)

mc1 <- mcmc(fm1$out)
plot(mc1)
summary(mc1)
rejectionRate(mc1)






cleanEx()
nameEx("scrUN")
### * scrUN

flush(stderr()); flush(stdout())

### Name: scrUN
### Title: Fit SCR model to data collected on unmarked animals.
### Aliases: scrUN
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(coda)

# Analysis of the Northern Parula data
data(nopa)
fm1 <- scrUN(nopaDat$y, nopaDat$X, M=100, obsmod="pois", niters=1000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.5,0.2,1))
plot(mcmc(fm1))

fm2 <- scrUN(nopaDat$y, nopaDat$X, M=80, obsmod="bern", niters=12000,
    xlims=c(0,24), ylims=c(0, 16), tune=c(0.4, 0.05, 15))
mc2 <- mcmc(fm2)
plot(mc2)
rejectionRate(mc2)


# Simulate data using dimensions of parula data

set.seed(34123)
N <- 20
X <- nopaDat$X
J <- nrow(X)
K <- 10 # ncol(nopaDat$y)
z <- array(NA, c(N, J, K))
#xlims <- c(0, 24)
#ylims <- c(0, 16)
xlims <- c(2, 22)
ylims <- c(2, 14)
s <- cbind(runif(N, xlims[1], xlims[2]),
           runif(N, ylims[1], ylims[2]))
p0 <- 0.8
sigma <- 1
for(j in 1:J) {
    dist <- sqrt((s[,1]-X[j,1])^2 + (s[,2]-X[j,2])^2)
    p <- p0 * exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        z[,j,k] <- rbinom(N, 1, p) # latent encounter histories
#        z[,j,k] <- rpois(N, p) # latent encounter histories
    }
}
n <- colSums(z) # observed counts = sum z over individuals
table(n)

plot(s, pch=16, col=4, asp=1, xlim=xlims)
rect(xlims[1], ylims[1], xlims[2], ylims[2])
points(X, pch="+")

fm3 <- scrUN(n, X, M=100, obsmod="pois", niters=5000,
    xlims=xlims, ylim=ylims, tune=c(0.1, 0.1, 5))

mc3 <- mcmc(fm3)
plot(mc3)
rejectionRate(mc3)
summary(window(mc3, start=10001))


# Use M-H to update z
fm4 <- scrUN(n, X, M=40, obsmod="bern", niters=50000,
    xlims=xlims, ylim=ylims, tune=c(0.03, 0.04, 2), zGibbs=FALSE)

mc4 <- mcmc(fm4)
plot(mc4)
rejectionRate(mc4)
summary(window(mc4, start=101001))


# Use Gibbs to update z (Why doesn't this work?????)
fm5 <- scrUN(n, X, M=40, obsmod="bern", niters=150000,
    xlims=xlims, ylim=ylims, tune=c(0.03, 0.05, 0.6), zGibbs=TRUE)

mc5 <- mcmc(fm5)
plot(mc5)
rejectionRate(mc5)
summary(window(mc5, start=3001))






cleanEx()
nameEx("scrbook-package")
### * scrbook-package

flush(stderr()); flush(stdout())

### Name: scrbook-package
### Title: Companion to the book "Spatial Capture Recapture"
### Aliases: scrbook-package scrbook
### Keywords: package

### ** Examples

~~ simple examples of the most important functions ~~



cleanEx()
nameEx("secr.bear")
### * secr.bear

flush(stderr()); flush(stdout())

### Name: secr.bear
### Title: analysis of the Ft. Drum black bear data using secr
### Aliases: secr.bear
### Keywords: secr AIC

### ** Examples


bear.ou<-secr.bear()
secr.bear$AIC.tab

## returns:
                                       model    detectfn npar    logLik
bear.b                     D~1 g0~bk sigma~1  halfnormal    4 -641.7215
bear.h2           D~1 g0~h2 sigma~h2 pmix~h2  halfnormal    6 -653.8382
bear.0exp                   D~1 g0~1 sigma~1 exponential    3 -663.9152
bear.B                      D~1 g0~b sigma~1  halfnormal    4 -677.6175
bear.Bt                 D~1 g0~b + t sigma~1  halfnormal   11 -668.3044
bear.sex  D~session g0~session sigma~session  halfnormal    6 -677.7151
bear.t                      D~1 g0~t sigma~1  halfnormal   10 -674.4134
bear.0                      D~1 g0~1 sigma~1  halfnormal    3 -686.2455
               AIC     AICc  dAICc AICwt
bear.b    1291.443 1292.395  0.000     1
bear.h2   1319.676 1321.776 29.381     0
bear.0exp 1333.830 1334.389 41.994     0
bear.B    1363.235 1364.187 71.792     0
bear.Bt   1358.609 1366.152 73.757     0
bear.sex  1367.430 1369.530 77.135     0
bear.t    1368.827 1374.938 82.543     0
bear.0    1378.491 1379.049 86.654     0




cleanEx()
nameEx("secr_wolverine")
### * secr_wolverine

flush(stderr()); flush(stdout())

### Name: secr_wolverine
### Title: analysis of the wolverine data using secr
### Aliases: secr_wolverine
### Keywords: ~kwd1 ~kwd2

### ** Examples

##
## The function executes a set of commands to conduct a
## likelihood analysis of the wolverine camera trapping data
## the function returns the results for the 2 x 2 habitat mask
##

## The function is currently defined by the following set of commands:


library("secr")
library("scrbook")
data("wolverine")
traps<-as.matrix(wolverine$wtraps)   #[,1:3]
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))

traps1<-as.data.frame(traps[,1:3])
traps1$x<-as.numeric(as.character(traps1$x))
traps1$y<-as.numeric(as.character(traps1$y))

# This seems to ignore the trap operation information
trapfile1<-read.traps(data=traps1,detector="proximity")

trapfile2<-scr2secr(scrtraps=traps,type="proximity")

wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)
wolvcapt1<-make.capthist(wolv.dat,trapfile1,fmt="trapID",noccasions=165)
wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)

gr<-(as.matrix(wolverine$grid2))
dimnames(gr)<-list(NULL,c("x","y"))
gr2<-read.mask(data=gr)

gr<-(as.matrix(wolverine$grid4))
dimnames(gr)<-list(NULL,c("x","y"))
gr4<-read.mask(data=gr)

gr<-(as.matrix(wolverine$grid8))
dimnames(gr)<-list(NULL,c("x","y"))
gr8<-read.mask(data=gr)

# run model without trap operation information
wolv.secr<-secr.fit(wolvcapt1,model=list(D~1, g0~1, sigma~1), buffer=20000)

# now use wolvcapt2 which has trap operation information
wolv.secr0<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000)
wolv.secr2<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr2)
wolv.secr4<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr4)
wolv.secr8<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr8)
 

# reported in the book chapter:
wolv.secr2




cleanEx()
nameEx("sim.pID.data")
### * sim.pID.data

flush(stderr()); flush(stdout())

### Name: sim.pID.data 
### Title: Simulation function for spatial mark-resight data
### Aliases: 'sim.pID.data '

### ** Examples

####example from Royle et al. 2013; Chapt. 15, sec. XXX
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



cleanEx()
nameEx("simMnSCR")
### * simMnSCR

flush(stderr()); flush(stdout())

### Name: simMnSCR
### Title: bleen
### Aliases: simMnSCR
### Keywords: ~kwd1 ~kwd2

### ** Examples


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




cleanEx()
nameEx("simPoissonSCR.fn")
### * simPoissonSCR.fn

flush(stderr()); flush(stdout())

### Name: simPoissonSCR
### Title: Simulatesion SCR data with Poisson observation model
### Aliases: simPoissonSCR.fn
### Keywords: ~kwd1 ~kwd2

### ** Examples






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





cleanEx()
nameEx("simSCR0")
### * simSCR0

flush(stderr()); flush(stdout())

### Name: simSCR0
### Title: Simulate some SCR data
### Aliases: simSCR0
### Keywords: ~kwd1 ~kwd2

### ** Examples


library("scrbook")
data<-simSCR0(discard0=TRUE,rnd=2013)


##out1<-SCR0bayes(data,M=200,engine="jags",ni=2000,nb=1000)
# Fit the model using WinBUGS:
out2<-SCR0bayes(data,M=200,engine="winbugs",ni=2000,nb=1000)









cleanEx()
nameEx("simScSCR.fn")
### * simScSCR.fn

flush(stderr()); flush(stdout())

### Name: simScSCR.fn
### Title: bleen
### Aliases: simScSCR.fn
### Keywords: ~kwd1 ~kwd2

### ** Examples


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




cleanEx()
nameEx("spatial.plot")
### * spatial.plot

flush(stderr()); flush(stdout())

### Name: spatial.plot
### Title: Makes a primative spatial plot using a grid of points shaded by
###   the attribute of interest.
### Aliases: spatial.plot

### ** Examples


spatial.plot <-
function(x,y,add=TRUE,cx=1){
 nc<-as.numeric(cut(y,20))
if(!add) plot(x,pch=" ")
 points(x,pch=20,col=terrain.colors(20)[nc],cex=cx)
image.scale(y,col=terrain.colors(20))

}




cleanEx()
nameEx("spcov")
### * spcov

flush(stderr()); flush(stdout())

### Name: spcov
### Title: Functions to generate spatial covariates using Kriging
### Aliases: spcov
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(B=1, v=20) {
    R <- data.frame(x=rep(seq(0, B, length=v), each=v),
                    y=rep(seq(0, B, length=v), times=v))
    elev <- function(x, y) (x-.5)+(y-.5) # Elevation is a function of x,y
    D<-e2dist1(R,R)
    V<-exp(-D/2)
    Vi<-solve(V)

    cov1<-t(chol(V))%*%rnorm(nrow(R))
    image(matrix(cov1,20,20))
    R$elev <- apply(R[,1:2], 1, function(x) elev(x[1], x[2]))
    R$cov1<-cov1-mean(cov1)
    cov1.fn<-function(newpt,cov1,cov1.coords=cbind(R$x,R$y),Vi){
        newpt<-matrix(newpt,ncol=2)
        k<- exp(-e2dist1(newpt,cov1.coords)/2)
        pred<-k%*%Vi%*%cov1
        as.numeric(pred)
    }
    return(cov1.fn)
  }



cleanEx()
nameEx("spiderplot")
### * spiderplot

flush(stderr()); flush(stdout())

### Name: spiderplot
### Title: Makes a spider plot of some SCR data.
### Aliases: spiderplot

### ** Examples


library("scrbook")

# make spiderplot for the beardata the hard way:
data(beardata)
X<-as.matrix(cbind(id=1:nrow(beardata$trapmat),beardata$trapmat))
opps<-matrix(1,nrow=nrow(X),ncol=8)
dimnames(opps)<-list(NULL,1:8)
X<-cbind(X,opps)
a<-SCR23darray(beardata$flat,X)
toad<-spiderplot(a,beardata$trapmat)
# now grab the distance from centroid variable
xcent<-toad$xcent
# see Chapter 3 of the book

# the easy way
spiderplot(beardata$bearArray,beardata$trapmat)

# for the wolverine data:
y3d <- SCR23darray(wolverine$wcaps,wolverine$wtraps)
spiderplot(y3d,wolverine$wtraps[,2:3])








cleanEx()
nameEx("tortoise")
### * tortoise

flush(stderr()); flush(stdout())

### Name: tortoise
### Title: desert tortoise data
### Aliases: tortoise
### Keywords: datasets

### ** Examples

data(tortoise)

data.mn<- formatDistData(tortoise,"Dist","Transect",dist.breaks=c(0:32))

library("unmarked")
 umf <- unmarkedFrameDS(y=as.matrix(data.mn), siteCovs=NULL, survey="line",
 dist.breaks=c(0:32), tlength=rep(1000, 120),
 unitsIn="m")
 summary(umf)

png("tortoise.png",width=7,height=7, units="in", res=400)
 hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)
dev.off()

m0<-distsamp(~1~1,umf,keyfun="halfnorm",output="density",unitsOut="ha")

x<-tortoise[,"Dist"]
nind<-sum(!is.na(x))
y<-rep(1,nind)
nz<-700
y<-c(y,rep(0,nz))
x<-c(x,rep(NA,nz))
z<-y

cat("
model{
alpha1~dunif(0,10)
sigma<- sqrt(1/(2*alpha1))
psi~dunif(0,1)

for(i in 1:(nind+nz)){
   z[i]~dbern(psi)    # DA Variables
   x[i]~dunif(0,B)    # B=strip width
   p[i]<-exp(logp[i])   # DETECTION MODEL
   logp[i]<-   -alpha1*(x[i]*x[i])
   mu[i]<-z[i]*p[i]
   y[i]~dbern(mu[i])  # OBSERVATION MODEL
 }

N<-sum(z[1:(nind+nz)])
D<- N/striparea  # area of transects
}
",file="dsamp.txt")
library("R2WinBUGS")
# density to be in units of ind/ha
data<-list(y=y,x=x,nz=nz,nind=nind,B=40,striparea=(120*1000*40*2/10000))
params<-list('alpha1','sigma','N','D','psi')
inits =  function() {list(z=z, psi=runif(1), alpha1=runif(1,0,.02) )}
fit = bugs(data, inits, params, model.file="dsamp.txt",working.directory=getwd(),    
       debug=T, n.chains=3, n.iter=3000, n.burnin=1000, n.thin=2)









cleanEx()
nameEx("wolvESA")
### * wolvESA

flush(stderr()); flush(stdout())

### Name: wolvESA
### Title: computes effective sample area plot for wolverine data
### Aliases: wolvESA
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (noargs = TRUE) 
{
    library("scrbook")
    data(wolverine)
    traps <- wolverine$wtraps
    y3d <- SCR23darray(wolverine$wcaps, wolverine$wtraps)
    g <- as.matrix(wolverine$grid2)
    SCRsmy(y3d)
    Sgrid <- g
    traplocs <- as.matrix(traps[, 2:3])
    mingridx <- min(traplocs[, 1])
    mingridy <- min(traplocs[, 2])
    traplocs[, 1] <- traplocs[, 1] - min(traplocs[, 1])
    traplocs[, 2] <- traplocs[, 2] - min(traplocs[, 2])
    traplocs <- traplocs/10000
    ntraps <- nrow(traplocs)
    MASK <- traps[, 4:ncol(traps)]
    ndays <- apply(MASK, 1, sum)
    Sgrid[, 1] <- Sgrid[, 1] - mingridx
    Sgrid[, 2] <- Sgrid[, 2] - mingridy
    Sgrid <- Sgrid/10000
    Dmat <- e2dist(traplocs, Sgrid)
    probs <- rep(1/nrow(Sgrid), nrow(Sgrid))
    p0 <- 0.05
    sigma <- 0.62
    D <- e2dist(traplocs, Sgrid)
    netp <- rep(NA, nrow(Sgrid))
    for (i in 1:nrow(Sgrid)) {
        d2 <- (traplocs[, 1] - Sgrid[i, 1])^2 + (traplocs[, 2] - 
            Sgrid[i, 2])^2
        pvec <- p0 * exp(-(1/(2 * sigma * sigma)) * d2)
        netp[i] <- 1 - prod((1 - pvec)^ndays)
    }
    y <- netp
    x <- Sgrid
    par(mar = c(3, 3, 3, 6))
    plot(x, pch = " ")
    nc <- as.numeric(cut(y, 10))
    cc <- topo.colors(10)
    points(x, pch = 20, col = cc[nc], cex = 1)
    image.scale(y, col = cc)
  }



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("wolvSCR0")
### * wolvSCR0

flush(stderr()); flush(stdout())

### Name: wolvSCR0
### Title: fit SCR0 to Wolverine data
### Aliases: wolvSCR0
### Keywords: ~kwd1 ~kwd2

### ** Examples


library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
toad<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100)




cleanEx()
nameEx("wolvSCR0Dss")
### * wolvSCR0Dss

flush(stderr()); flush(stdout())

### Name: wolvSCR0Dss
### Title: fit SCR0 to wolverine data using discrete state-space
### Aliases: wolvSCR0Dss
### Keywords: ~kwd1 ~kwd2

### ** Examples


library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
g<-as.matrix(wolverine$grid8)
bln<-wolvSCR0Dss(y3d,traps,nb=1000,ni=2000,M=100,Sgrid=g,engine="jags",area=8)




cleanEx()
nameEx("wolvSCR0Dssv2")
### * wolvSCR0Dssv2

flush(stderr()); flush(stdout())

### Name: wolvSCR0Dssv2
### Title: SCR0 with discrete state-space for wolverine data
### Aliases: wolvSCR0Dssv2
### Keywords: ~kwd1 ~kwd2

### ** Examples




library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
g<-as.matrix(wolverine$grid8)
bln<-wolvSCR0Dss(y3d,traps,nb=1000,ni=2000,M=160,Sgrid=g,engine="bugs",area=8)
# time trial with/without distance matrix calculation in BUGS//JAGS compare result
# run regular grid with mask?
bln2<-wolvSCR0Dss(y3d,traps,nb=1000,ni=2000,M=160,Sgrid=g,engine="jags",area=8)







cleanEx()
nameEx("wolvSCR0ms")
### * wolvSCR0ms

flush(stderr()); flush(stdout())

### Name: wolvSCR0ms
### Title: fits multiple models to wolverine data
### Aliases: wolvSCR0ms.fn
### Keywords: ~kwd1 ~kwd2

### ** Examples


# library("R2jags")
# library("scrbook")
# data(wolverine)
#
# These functions load libraries and data -- please read through them
# Default to run with jags()

## Run here with too FEW iterations -- increase to at least 10000 total
##
## This takes several hours to run
## 
 toad0<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=0)
 toad1<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=1)
 toad2<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=2)
 toad3<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=3)
 toad4<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=4)
 toad5<-wolvSCR0ms(nb=1000,ni=2000,buffer=2,M=200,model=5)

# summary function to extract DIC info
dev.smy<-function(a){
b<-print(a)
dev<-mean(a$BUGSoutput$sims.list$deviance)
c(dev,b$pD,b$DIC,a$BUGSoutput$DICbyR)
}

# make a table
rbind(dev.smy(toad0),dev.smy(toad1),dev.smy(toad2),dev.smy(toad3),dev.smy(toad4))

# summary function for parameters
parms.smy<-function(a){
b<-print(a)
b<- b$summary[,1:2]
nam<-c("D","N","alpha0","alpha.sex","beta","beta.sex","sigma","psi","psi.sex","deviance")
round(b[nam,],2)
}

# make a table
cbind(parms.smy(toad0),parms.smy(toad1),parms.smy(toad2),parms.smy(toad3),parms.smy(toad4))










cleanEx()
nameEx("wolvSCR0ms2")
### * wolvSCR0ms2

flush(stderr()); flush(stdout())

### Name: wolvSCR0ms2
### Title: model selection of encounter models
### Aliases: wolvSCR0ms2
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


toad<-wolvSCR0ms2(nb=500,ni=1000,buffer=2,M=121,engine="winbugs")

mprobs<- table(toad$sims.list$mod)
mprobs<-mprobs/sum(mprobs)





cleanEx()
nameEx("wolvSCR0pois.fn")
### * wolvSCR0pois.fn

flush(stderr()); flush(stdout())

### Name: wolvSCR0pois.fn
### Title: wolverine SCR model with Poisson observation model
### Aliases: wolvSCR0pois.fn
### Keywords: ~kwd1 ~kwd2

### ** Examples

library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray.fn(wolverine$wcaps,wolverine$wtraps)
toad<-wolvSCR0.fn(y3d,traps,nb=1000,ni=6000,delta=2,M=200)
toad2<-wolvSCR0pois.fn(y3d,traps,nb=1000,ni=6000,delta=2,M=200)



cleanEx()
nameEx("wolverine")
### * wolverine

flush(stderr()); flush(stdout())

### Name: wolverine
### Title: wolverine data from Audrey Magoun
### Aliases: wolverine
### Keywords: datasets

### ** Examples

data(wolverine)
## maybe str(wolverine) ; plot(wolverine) ...



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
