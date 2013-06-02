
### before running this code, put the functions at the end of this script
### into your R workspace
###


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
tmp1<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

# use mean "s" instead of estimating it
tmp3<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z),stel=sbar)

# no SCR data, s is random. Here there are 2 extra parameters that are not estimated: start[1] and start[4]
tmp4<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=NULL,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

# Fits SCR model with isotropic Gaussian encounter model
tmp5<- nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=rep(0,ntraps),G=gr)

####
####  NOW WE FIT SOME MODELS WITH COVARIATE ON LOG-DENSITY
####
####

# Fits SCR model with isotropic Gaussian encounter model log(D(x)) = covariate
#  AND with cloglog(p(x)) = covariate
# note one of the starting values is a dead parameter so exclude it from nlm() output
tmp6<- nlm(intlik3rsf.v3,c(-3,log(3),.17,4,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,
zall=as.vector(z),hessian=TRUE)

# use telemetry data and have covariate on log(D(x))
tmp7<-nlm(intlik3rsf.v3,tmp6$estimate,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,
zall=as.vector(z),hessian=TRUE)

# Fit a model with log(D(x)) = z(x) but NO COVARIATE on p(x)
tmp8<- nlm(intlik3rsf.v3,tmp7$estimate,y=y,K=K,X=X,
ztrap=rep(0,length(raster.point)),G=gr,zall=as.vector(z),hessian=TRUE)



###
###
### put all the functions below this line into your R workspace
###
###

e2dist<-
function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

"rot" <-
function (m)
{
# takes a matrix m and rotates it such that image(rot(m)) is plotted "as you look at it"
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

 spatial.plot<-
function(x,y){
 nc<-as.numeric(cut(y,20))
 plot(x,pch=" ")
 points(x,pch=20,col=topo.colors(20)[nc],cex=2)
 ###image.scale(y,col=topo.colors(20))
}

### This is the likelihood function
### It computes several versions of the likelihood depending on the arguments specified
### see the 5 examples above

intlik3rsf.v2 <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL,stel=NULL){
# start = vector of length 4 = starting values
###
### ORDER OF STARTING VALUES
### (intercept for cloglog-scale detection, log(sigma), coefficient on z(x), log(n0))
### where N = n_observed + n0  i.e., n0 = number of uncaptured individuals
### to transform back take , e.g., nrow(y)+exp(tmp1$estimate[4]) = Nhat
###
###
# y = nind x ntraps encounter matrix
# K = how many samples?
# X = trap locations
# ztrap = covariate value at trap locations
# zall = all covariate values for all nG pixels
# ntel = nguys x nG matrix of telemetry fixes in each nG pixels
# stel = home range center of telemetered individuals, IF you wish to estimate it. Not necessary

nG<-nrow(G)
D<- e2dist(X,G)

alpha0<-start[1]
sigma<- exp(start[2])
alpha2<- start[3]
n0<-    exp(start[4])
a0<- 1

if(!is.null(y)){
loglam<-   alpha0  -(1/(2*sigma*sigma))*D*D + alpha2*ztrap  # ztrap recycled over nG
probcap<- 1-exp(-exp(loglam))
#probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
ymat<-y
ymat<-rbind(y,rep(0,ncol(y)))
lik.marg<-rep(NA,nrow(ymat))
for(i in 1:nrow(ymat)){
Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
lik.cond<- exp(colSums(Pm))
lik.marg[i]<- sum( lik.cond*(1/nG) )
}
nv<-c(rep(1,length(lik.marg)-1),n0)
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
part2<- sum(nv*log(lik.marg))
out<-  -1*(part1+ part2)
}
else{
out<-0
}

if(!is.null(ntel) & !is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(stel,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- lam/denom  # each column is the probs for a guy at column [j]

tel.loglik<-  -1*sum(  ntel*log(probs) )

out<- out  + tel.loglik
}

if(!is.null(ntel) & is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(G,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- t(lam/denom)  # each column is the probs for a guy at column [j]
marg<- as.vector(rowSums(exp(ntel%*%log(probs))/nG ))

tel.loglik<- -1*sum(log(marg))

out<- out  + tel.loglik
}

out
}







intlik3rsf.v3 <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL,stel=NULL){
#
# this version of the code handles a covariate on log(Density). This is starting value 5
#
# start = vector of length 5 = starting values
# y = nind x ntraps encounter matrix
# K = how many samples?
# X = trap locations
# ztrap = covariate value at trap locations
# zall = all covariate values for all nG pixels
# ntel = nguys x nG matrix of telemetry fixes in each nG pixels
# stel = home range center of telemetered individuals, IF you wish to estimate it. Not necessary

nG<-nrow(G)
D<- e2dist(X,G)

alpha0<-start[1]
sigma<- exp(start[2])
alpha2<- start[3]
n0<-    exp(start[4])
beta<- start[5]
a0<- 1
if(!is.null(zall)){
 psi<- exp(beta*zall)
 psi<-psi/sum(psi)
}
else{
psi<-rep(1/nG,nG)
}
if(!is.null(y)){
loglam<-   alpha0  -(1/(2*sigma*sigma))*D*D + alpha2*ztrap  # ztrap recycled over nG


probcap<- 1-exp(-exp(loglam))
#probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
ymat<-y
ymat<-rbind(y,rep(0,ncol(y)))
lik.marg<-rep(NA,nrow(ymat))
for(i in 1:nrow(ymat)){
Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
lik.cond<- exp(colSums(Pm))
lik.marg[i]<- sum( lik.cond*psi )
}
nv<-c(rep(1,length(lik.marg)-1),n0)
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
part2<- sum(nv*log(lik.marg))
out<-  -1*(part1+ part2)
}
else{
out<-0
}

if(!is.null(ntel) & !is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(stel,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- lam/denom  # each column is the probs for a guy at column [j]

tel.loglik<-  -1*sum(  ntel*log(probs) )

out<- out  + tel.loglik
}

if(!is.null(ntel) & is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(G,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- t(lam/denom)  # each column is the probs for a guy at column [j]
temp<-exp(ntel%*%log(probs))  # Ntel x nG matrix

marg<- as.vector(rowSums(  temp*psi ))


tel.loglik<- -1*sum(log(marg))

out<- out  + tel.loglik
}





out
}


####
####
####
#### black bear analysis
####
####
####

grid<-read.csv("AvgElevation1000by1000.csv")
z<-grid[,4]
z<- (z-mean(z))/sqrt(var(z))
gr<-grid[,7:8]/10000
spatial.plot(gr,z)

##### trap locations J= 122
traplocs<-read.csv("traplocs_1.csv")
head(traplocs)
J<-nrow(traplocs)
attach(traplocs)
plot(easting,northing)
X<-traplocs[,c("easting","northing")]/10000

##### nonspatial encounter histories for n=30
nonSpat<-read.csv("enchist_nonspat_1.csv")
head(nonSpat)
n<-nrow(nonSpat);n
K<-5
enchist_nonspat<-nonSpat[,3:7]
head(enchist_nonspat)

      ##############################
      #### n=30 ; J=122  ; K=5 #####
      ##############################

#### 3d encounter histories for n=30

enchist_spat<-array(0,dim=c(n,J,K))

enchists<-read.csv("enchist_1.csv")[,-2]
for (i in 1:nrow(enchists)){
  enchist_spat[enchists[i,1],enchists[i,2],enchists[i,3]]<-1
}

 apply(enchist_spat,c(1,2),sum)->ff

 dim(ff)
[1]  30 122
 table(ff)
ff
   0    1    2 
3628   30    2 

2 guys caught twice in 2 traps
> apply(ff,1,sum)
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 4 1 1 1 2 1 1 1 1
 
one guy caught 4 times total
one guy caught 2 times total
 == NO informatinon about model parameters! 

y<-ff
K<-5   ## (weekly periods?)

spatial.plot(grid[,c("X_UTM","Y_UTM")],z)

teldata<- read.csv("Summer2011_locations_3bears_Andy_wUTM.csv")
teldata<-teldata[seq(1,nrow(teldata),10),]
tel.id<-as.numeric(factor(teldata[,1]))
teldata<-teldata[,7:8]/10000
points(teldata,pch=20)
## THIS CAN CRASH if raster is huge
d<-e2dist(teldata,gr)
raster.point<-rep(NA,nrow(d))
for(i in 1:nrow(d)){
raster.point[i]<-(1:ncol(d))[d[i,]==min(d[i,])]
}
ntel<-matrix(0,nrow=3,ncol=nrow(grid))
for(i in 1:3){
x<-table(raster.point[tel.id==i])
ntel[i,as.numeric(names(x))]<-x
}




d<-e2dist(X,gr)
raster.point.trap<-rep(NA,nrow(X))
for(i in 1:nrow(d)){
raster.point.trap[i]<-(1:ncol(d))[d[i,]==min(d[i,])]
}

gr<-as.matrix(gr)
tel.sum<-ntel%*%gr
tel.n<-apply(ntel,1,sum)
sbar<- tel.sum/matrix(tel.n,ncol=2,nrow=nrow(tel.sum),byrow=TRUE)

# Basic SCR model with RSF covariate at trap locations. 
tmp1<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,gradtol=.01)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2<-nlm(intlik3rsf.v2,c(-4.7,log(1),-.3,5),y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z))

# use mean "s" instead of estimating it 
tmp3<-nlm(intlik3rsf.v2,c(-7, 7.2, -.2, 4),y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z),stel=sbar,gradtol=.1)

# no SCR data, s is random. Here there are 2 extra parameters that are not estimated: start[1] and start[4]
tmp4<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=NULL,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z))

# Fits SCR model with isotropic Gaussian encounter model
tmp5<- nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=rep(0,length(raster.point.trap)),G=gr,gradtol=.01)


o1<-nlm(intlik3rsf.v2,c(-6.7,-.15,4.9),y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z),sigma=1)

o05<-nlm(intlik3rsf.v2,o1$estimate,y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z),sigma=.5)
o15<-nlm(intlik3rsf.v2,o1$estimate,y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z),sigma=1.5)

o20<-nlm(intlik3rsf.v2,o15$estimate,y=y,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z),sigma=2)


library("rgeos")
png("studyarea2.png",width=7,height=7, units="in", res=400)

buff<- 1
p1<-Polygon(rbind(X,X[1,]))
p2<-Polygons(list(p1=p1),ID=1)
p3<-SpatialPolygons(list(p2=p2))
p1ch<-gConvexHull(p3)
 bp1<-(gBuffer(p1ch, width=buff))
 plot(bp1, col='gray')
 plot(p1ch, border='black', lwd=2, add=TRUE)
 gArea(bp1)
 #points(gr,pch=".")
 points(X,pch=20)


library(maptools)
library(sp)
library(scrbook)

Pcoord<-SpatialPoints(gr)
PinPoly<-over(Pcoord,bp1)  ### determine if each point is in polygon
mask<-as.numeric(!is.na(PinPoly))  ## convert to binary 0/1
gr2<-gr[mask==1,]
z2<-z[mask==1]
#points(gr2,pch=".")
dev.off()

teldata<- read.csv("Summer2011_locations_3bears_Andy_wUTM.csv")
teldata<-teldata[seq(1,nrow(teldata),10),]
tel.id<-as.numeric(factor(teldata[,1]))
teldata<-teldata[,7:8]/10000
points(teldata,pch=20)
## THIS CAN CRASH if raster is huge
d<-e2dist(teldata,gr2)
raster.point2<-rep(NA,nrow(d))
for(i in 1:nrow(d)){
raster.point2[i]<-(1:ncol(d))[d[i,]==min(d[i,])]
}
ntel2<-matrix(0,nrow=3,ncol=nrow(gr2))
for(i in 1:3){
x<-table(raster.point2[tel.id==i])
ntel2[i,as.numeric(names(x))]<-x
}


d<-e2dist(X,gr2)
raster.point.trap2<-rep(NA,nrow(X))
for(i in 1:nrow(d)){
raster.point.trap2[i]<-(1:ncol(d))[d[i,]==min(d[i,])]
}





# Basic SCR model with RSF covariate at trap locations. 
tmp1c.1<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z2[raster.point.trap2],G=gr2,hessian=TRUE)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2c.1<-nlm(intlik3rsf.v2,c(-4.7,log(1),-.3,5),y=y,K=K,X=X,
ztrap=z2[raster.point.trap2],G=gr2,ntel=ntel2,zall=as.vector(z2),hessian=TRUE)

# no SCR data, s is random. Here there are 2 extra parameters that are not estimated: start[1] and start[4]
#tmp4<-nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=NULL,K=K,X=X,ztrap=z[raster.point.trap],G=gr,ntel=ntel,zall=as.vector(z))

# Fits SCR model with isotropic Gaussian encounter model
tmp5c.1<- nlm(intlik3rsf.v2,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=rep(0,length(raster.point.trap2)),G=gr2,hessian=TRUE)


# Fits SCR model with isotropic Gaussian encounter model
tmp6c.1<- nlm(intlik3rsf.v3,c(-3,log(3),.17,4,0),y=y,K=K,X=X,ztrap=z2[raster.point.trap2],G=gr2,zall=as.vector(z2),hessian=TRUE)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2d.1<-nlm(intlik3rsf.v3,c(-4.7,-.8,-.3,5,0),y=y,K=K,X=X,ztrap=z2[raster.point.trap2],G=gr2,ntel=ntel2,
zall=as.vector(z2),hessian=TRUE)
tmp5d.1<- nlm(intlik3rsf.v3,c(-3.8,-1.2,0,5.4,0),y=y,K=K,X=X,
ztrap=rep(0,length(raster.point.trap2)),G=gr2,zall=as.vector(z2),hessian=TRUE)



n=33 [season 1 1st half?]

                 alpha0   log(sigma)    alpha2      log(n0)  beta        -loglik
SCR+p(x)      -2.8561676 -1.1174638  0.1747187   4.1395461                122.738
   SE          0.3899063  0.1389833  0.2477921   0.3656961
SCR           -2.729194   -1.122389   1.000000   4.109886                 122.990   
   SE          0.3453705   0.1403783             0.3618065
SCR+D(x)      -2.715320  -1.133076   0.000000   4.113903   1.247256       118.007
   SE          0.3526155 0.1394352              0.3575286 0.4083330       
SCR+p(x)+D(x) -2.4838347 -1.1567458 -0.3842881  4.2547317  1.5710664      117.075
               0.3910420  0.1421062  0.2760693  0.3767896  0.4630096

SCR+RSF       -3.0676938 -0.8141204 -0.2810946   3.8841581               1271.739
   SE          0.27218129  0.03640307 0.11759909 0.36255170
SCR+RSF+D(x)  -3.0701403 -0.8100523 -0.3706265  4.0284284  1.272629      1266.700
   SE          0.27199799 0.03683849 0.12387969 0.36606116 0.411030



use<- exp(-.371*z2)
odds<- use/exp(-.371*mean(z2) )
# odds ratio of use of x relative to average pixel  at same distance from s

area<- gArea(bp1)*100
area.per.pixel<- area/nrow(gr2)

Nhat<-exp(4.028)+nrow(y)

Dhat<- exp(1.27*z2)
Dhat<-(Dhat/sum(Dhat))*Nhat   # individuals per pixel
Dhat<- (Dhat/area.per.pixel)*100

png("spaceusage.png",width=7,height=7,units="in",res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr2,odds)
dev.off()

png("density.png",width=7,height=7, units="in", res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr2,Dhat)
dev.off()

png("elev_captures.png",width=7,height=7,units="in",res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr2,z2)
tmp<-X[col(y)[y>0],]
tmp<-tmp + rnorm(prod(dim(tmp)),0,.1)
points(tmp,pch=20)  # traps where captures happened.
dev.off()





###
###
###  Simulation
###
###

## the following block of code makes up a covariate as a spatially correlated
## noise field, with an exponential spatial correlation function
set.seed(1234)
gr<-as.matrix(expand.grid(1:40,1:40))
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))%*%rnorm(1600)
spatial.plot(gr,z)

###
### Set some parameter values
###
alpha0 <- -2
sigma<- 2
beta<- 1
Ntel<-4      # number of individuals with telemeters
nsim<-10
Nfixes<-20   # number of telemetry fixes per individual
N<- 100      # population size 


## create some matrices to store output. Here the naming convention is:  4 = number of telemetered
## individuals and 20 = number of fixes. 
## simout0 = basic SCR model with no telemetry guys
simout0<-simout4.20<-simout4.20b<-matrix(NA,nrow=nsim,ncol=6)


for(sim in 1:nsim){
cat("sim: ",sim,fill=TRUE)

# simulate activity centers of all N individuals
Sid<- sample(1:1600,N,replace=TRUE)
# and coordinates
S<-gr[Sid,]
# now draw centers of telemetered individuals
# have to draw telemetry guys interior or else make up more landscape -- 
# can't have truncated telemetry obs

poss.tel<- S[,1]>5 & S[,1]<35 & S[,2]>5 & S[,2]<35
tel.guys<-sample(Sid[poss.tel],Ntel)
sid<-tel.guys 
stel<-gr[sid,]
if(1==2){
stest<-source("stest.R")$value
Ntel<-8
Nfixes<-120
sid<-stest$id
stel<-stest$s
}

# make up matrix to store RSF data
n<-matrix(NA,nrow=Ntel,ncol=1600)

i<-1  # do this for each telemetered guy to simulate a number of fixes.
      # note that n = 0 for most of the landscape
par(mfrow=c(3,3))
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
   d<- Dmat[sid[i],]
   lam<- exp(1 - (1/(2*sigma*sigma))*d*d + beta* z) 
   n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
#  n[i,]<- rpois(1600,lam)
   par(mar=c(3,3,3,6))
   lammat[i,]<-lam
   img<- matrix(lam,nrow=40,ncol=40,byrow=FALSE)
   image(1:40,1:40,rot(img),col=terrain.colors(10))
}

if(1==2){

png("habitat.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1),mar=c(3,3,6,6))
image(1:40,1:40,rot(matrix(z,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
image.scale(z,col=terrain.colors(10))
points(stel,pch=20,cex=2)
dev.off()


png("homeranges8.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1))
tot<- apply(lammat,1,sum)
lammat<-lammat/tot
lamtot<-apply(lammat,2,sum)
image(1:40,1:40,rot(matrix(lamtot,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
points(stel,pch=20,cex=2)
dev.off()
}

## now lets simulate some SCR data on a bunch of guys:

# make a trap array
X<-  cbind(  sort(rep( seq(5,35,5),7)), rep( seq(5,35,5),7))
ntraps<-nrow(X)
raster.point<-rep(NA,nrow(X))
for(j in 1:nrow(X)){  # which piont in the raster is the trap? must be raster points
 raster.point[j]<- (1:1600)[ (X[j,1]==gr[,1]) & (X[j,2] == gr[,2])]
}

points(X,pch=20,cex=2)

D<- e2dist(S,X)  ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<-   alpha0  -(1/(2*sigma*sigma))*D*D + beta*Zmat
p<- 1-exp(-exp(loglam))

## Now simulate SCR data

K<- 10
y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:N){
y[i,]<- rbinom(ntraps,K,p[i,])
}

cap<-apply(y,1,sum)>0
y<-y[cap,]

sbar<- (n%*%gr)/as.vector(n%*%rep(1,nrow(gr)))

tmp1<-optim(c(-3,log(3),1,0,1),intlik3rsf.v2,y=y,K=K,X=X,ztrap=z[raster.point],G=gr)
simout0[sim,]<-c(nrow(y),tmp1$estimate)

# use telemetry s random
tmp2<-optim(c(-3,log(3),1,0,1),intlik3rsf.v2,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))
simout4.20[sim,]<-c(nrow(y),tmp2$estimate)

# use mean "s" instead of estimating it 
tmp3<-optim(c(-3,log(3),1,0,1),intlik3rsf.v2,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z),stel=sbar)
simout4.20b[sim,]<-c(nrow(y),tmp3$estimate)

# no SCR data, s is random
tmp4<-optim(c(-3,log(3),1,0,1),intlik3rsf.v2,y=NULL,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))
simout4.20c[sim,]<-c(NA,tmp4$estimate)
}













