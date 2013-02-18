\name{intlik3rsf}
\alias{intlik3rsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SCR+RSF integrated likelihood%%  ~~function to do ... ~~
}
\description{
  computes the likelihood from royle et al. 2013 mee paper
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlik3rsf(start = NULL, y = y, K = NULL, X = traplocs, ztrap, G, ntel = NULL, zall = NULL, stel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{ztrap}{
%%     ~~Describe \code{ztrap} here~~
}
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{ntel}{
%%     ~~Describe \code{ntel} here~~
}
  \item{zall}{
%%     ~~Describe \code{zall} here~~
}
  \item{stel}{
%%     ~~Describe \code{stel} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as


### before running this code, put the functions at the end of this script
### into your R workspace
###
library("scrbook")

## the following block of code makes up a covariate as a spatially correlated
## noise field, with an exponential spatial correlation function
set.seed(1234)
gr<-expand.grid(1:40,1:40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))%*%rnorm(1600)
spatial.plot(gr,z)


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
sbar<- (n%*%gr)/as.vector(n%*%rep(1,nrow(gr)))

###
### ORDER OF STARTING VALUES
### (intercept for cloglog-scale detection, log(sigma), coefficient on z(x), log(n0))
### where N = n_observed + n0  i.e., n0 = number of uncaptured individuals
### to transform back take , e.g., nrow(y)+exp(tmp1$estimate[4]) = Nhat
###
###
# Basic SCR model with RSF covariate at trap locations.
tmp1<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

# use mean "s" instead of estimating it
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
