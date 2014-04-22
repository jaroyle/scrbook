\name{intlik3rsfD}
\alias{intlik3rsfD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Model for combined SCR and telemetry data on resource selection.
  Allows for inclusion of a covariate on density. 
}
\description{
   Computes the likelihood for combined SCR and telemetry data with a resource selection
  covariate. Uses only a single resource selection covariate (this is
  required input in the present version of the function) and the
  Gaussian hazard encounter probability model:

  p(x,s) = 1-exp(-lambda(x,s))

  log(lambda(x,s)) = alpha0 - (1/(2*sigma*sigma))*D(x,s)^2  + alpha2*covariate(x)

  Covariates on density are modeled as:

  log(D(s)) = beta*covariate(s)

  (note: no intercept in the density model because the model is
  described conditional on total population size N for the state-space).
  
}
\usage{
intlik3rsfD(start = NULL, y = y, K = NULL, X = traplocs, ztrap, G, ntel = NULL, zall = NULL, stel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
     Starting values for the parameter vector. Has length = 5:
    alpha0 <- start[1]    
    sigma <- exp(start[2])
    alpha2 <- start[3]
    n0 <- exp(start[4]) # n0 = number of uncaptured individuals
    beta<- start[5]     # Density covariate
  }
   \item{y}{ nind x ntraps matrix of encounter frequencies (out of K)
%%     ~~Describe \code{y} here~~
}
  \item{K}{
Number of sample occasions%%     ~~Describe \code{K} here~~
}
  \item{X}{
ntraps x 2 matrix of trap coordinates %%     ~~Describe \code{X} here~~
}
  \item{ztrap}{
Covariate values at the trap locations, ntraps x 1 vector %%     ~~Describe \code{ztrap} here~~
}
  \item{G}{
nG x 2 matrix of state-space grid cells -- should be the raster coordinates (center points) %%     ~~Describe \code{G} here~~
}
  \item{ntel}{
nTel x nG matrix of telemetry encounter frequencies in each grid
cell. nTel = number of telemetered individuals, nG = number of raster
grid cells. %%     ~~Describe \code{ntel} here~~
}
  \item{zall}{
Raster covariate, nG x 1 vector %%     ~~Describe \code{zall} here~~
}
  \item{stel}{ optional nTel x 2 matrix of location centroids for each
    individual. If provided this is used as an estimate of s[i] for the
    telemetered individuals and the likelihood for telemetered
    individuals is not marginalized over s[i]
}
 }
\details{
This function illustrates calculation of the combined likelihood for
spatial capture-recapture data with some telemetered individuals.
The model assumes independence of the two samples and so does not
accommodate the observation of telemetered individuals in some
traps (which probably induces a very slight loss of efficiency).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
   Chapter 13 of the SCR book.
  
  Royle, J. A., Chandler, R. B., Sun, C. C., and Fuller, A. K. (2013).
  Integrating resource selection information with spatial
  capture-recapture. Methods in Ecology and Evolution, 4(6), 520-530.
}
\author{
Andy Royle, aroyle@usgs.gov
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
intlik3rsf 
}
\examples{


library("scrbook")

## the following block of code makes up a covariate as a spatially correlated
## noise field, with an exponential spatial correlation function
set.seed(1234)
gr<-expand.grid(1:40,1:40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))\%*\%rnorm(nrow(gr))
par(mar=c(3,3,3,6))
spatial.plot(gr,z,cx=3,col=terrain.colors(20))


###
### Set some parameter values
###
alpha0 <- -2
sigma<- 2
alpha2<- 1
Ntel<-4      # number of individuals with telemeters
nsim<-100
Nfixes<-20   # number of telemetry fixes per individual
N<- 100      # population size


# simulate activity centers of all N individuals
Sid<- sample(1:1600,N,replace=TRUE)
# and coordinates
S<-gr[Sid,]
# now simulate activity centers of telemetered individuals
# have to draw telemetry guys interior or else make up more landscape --
# can't have truncated telemetry obs
#
poss.tel<- S[,1]>5 & S[,1]<35 & S[,2]>5 & S[,2]<35
tel.guys<-sample(Sid[poss.tel],Ntel)
sid<-tel.guys
stel<-gr[sid,]

# make up matrix to store RSF data
n<-matrix(NA,nrow=Ntel,ncol=1600)

# for each telemetered guy simulate a number of fixes.
# note that n = 0 for most of the landscape
par(mfrow=c(3,3))
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
   d<- Dmat[sid[i],]
   lam<- exp(1 - (1/(2*sigma*sigma))*d*d + alpha2* z)
   n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
   par(mar=c(3,3,3,6))
   lammat[i,]<-lam
   img<- matrix(lam,nrow=40,ncol=40,byrow=FALSE)
   image(1:40,1:40,rot(img),col=terrain.colors(10))
}

## now lets simulate some SCR data on a bunch of guys:

# make a trap array
X<-  cbind(  sort(rep( seq(5,35,5),7)), rep( seq(5,35,5),7))
ntraps<-nrow(X)
raster.point<-rep(NA,nrow(X))
for(j in 1:nrow(X)){  # which point in the raster is the trap? must be raster points
 raster.point[j]<- (1:1600)[ (X[j,1]==gr[,1]) & (X[j,2] == gr[,2])]
}
points(X,pch=20,cex=2)

D<- e2dist(S,X)  ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<-   alpha0  -(1/(2*sigma*sigma))*D*D + alpha2*Zmat
p<- 1-exp(-exp(loglam))

## Now simulate SCR data

K<- 10
y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:N){
   y[i,]<- rbinom(ntraps,K,p[i,])
}

cap<-apply(y,1,sum)>0

y<-y[cap,]
gr<-as.matrix(gr)
sbar<- (n\%*\%gr)/as.vector(n\%*\%rep(1,nrow(gr)))

###
### ORDER OF STARTING VALUES
### (intercept for cloglog-scale detection, log(sigma), coefficient on z(x), log(n0))
### where N = n_observed + n0  i.e., n0 = number of uncaptured individuals
### to transform back take , e.g., nrow(y)+exp(tmp1$estimate[4]) = Nhat
###
###


### Fit basic SCR model with RSF covariate at trap locations.

tmp1<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr)

## Note: this will generate some non-critical numerical warnings. You
##   can use optim() instead which won't generate those warnings. 

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

# use mean telemetry location as estimate of "s" instead of marginalizing it out
tmp3<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z),stel=sbar)

# no SCR data, s is random. Here there are 2 extra parameters that are not estimated: start[1] and start[4]
tmp4<-nlm(intlik3rsf,c(-3,log(3),1,0),y=NULL,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

# Fits SCR model with isotropic Gaussian encounter model
tmp5<- nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=rep(0,ntraps),G=gr)

####
####  NOW FIT SOME MODELS WITH COVARIATE ON LOG-DENSITY
####
####

# Fits SCR model with isotropic Gaussian encounter model log(D(x)) = covariate
#  AND with cloglog(p(x)) = covariate
# note one of the starting values is a dead parameter so exclude it from nlm() output
tmp6<- nlm(intlik3rsfD,c(-3,log(3),.17,4,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,
        zall=as.vector(z),hessian=TRUE)

# use telemetry data and have covariate on log(D(x))
tmp7<-nlm(intlik3rsfD,tmp6$estimate,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,
          zall=as.vector(z),hessian=TRUE)

# Fit a model with log(D(x)) = z(x) but NO COVARIATE on p(x)
tmp8<- nlm(intlik3rsfD,tmp7$estimate,y=y,K=K,X=X,
        ztrap=rep(0,length(raster.point)),G=gr,zall=as.vector(z),hessian=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
