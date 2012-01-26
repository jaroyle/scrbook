
###
### Kimmy's SCR stuff starts here
####

library("rgeos")
library("shapefiles")
library("gdistance")
library("maptools")
library("spatstat")

setwd("C:\\Documents and Settings\\Kgazenski\\My Documents\\Andy\\streambuffering")

streamnet<- readShapeLines("clip_pwrcstreamnet")
plot(streamnet, col="blue")

streambuf<- gBuffer(streamnet, width=0.002)
plot(streambuf)


pspstreamnet<- as.psp(streamnet)
streampoints<- pointsOnLines(pspstreamnet)
plot(streampoints)

plot(streambuf)
plot(streampoints, pch=".", col="blue", add=TRUE)



pg<-gBuffer(streamnet, width=0.002)
pg.coords<- pg@polygons[[1]]@Polygons[[1]]@coords
 # note: can you believe this shit?
#pg.coords<-pg.coords*3



xg<-seq(-76.858,-76.68,,30)
yg<-seq(38.97,39.16,,30)


delta<-mean(diff(xg))
pts<- cbind(sort(rep(xg,30)),rep(yg,30))
points(pts,pch=20)

in.pts<-point.in.polygon(pts[,1],pts[,2],pg.coords[,1],pg.coords[,2])
points(pts[in.pts==1,],pch=20,col="red")

cost<-rep(NA,nrow(pts))
cost[in.pts==1]<-1   # low cost to move among pixels but not 0
cost[in.pts!=1]<-10000   # high cost 


library("raster")
 r<-raster(nrows=30,ncols=30)
 projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
 extent(r)<-c(-76.858-delta/2,-76.68+delta/2,38.97-delta/2,39.16+delta/2)
 #extent(r)<-c(.5,30.5,.5,30.5)


#values(r)<-rot(rot(rot(matrix(cost,30,30,byrow=TRUE))))
values(r)<-matrix(cost,30,30,byrow=FALSE)
par(mfrow=c(1,1))
rx<-flip(r,direction="y")
plot(rx)
points(pts,pch=20,cex=.4)

points(streampoints, pch=".", col="blue")

rx<-1/rx

## use max = doesn't count moving through boundary pixel
tr1<-transition(rx,transitionFunction=max,directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
costs1<-costDistance(tr1CorrC,pts)

sp<-rasterToPoints(rx,spatial=TRUE)

outD<-as.matrix(costs1)

points(matrix(pts[116,],ncol=2),col="red")
points(matrix(pts[111,],ncol=2),col="blue")

outD[116,100:120]

plot(pts,pch=".")
points(pts[in.pts==1,],pch=20,col="red")


get.traplocs<-function(ntraps){
loc<-NULL
locid<-NULL
for(i in 1:ntraps){
x<-locator(1)
d<- sqrt((x$x - pts[,1])^2 + (x$y - pts[,2])^2 )
kp<- d==min(d)

tmp<-pts[kp,]
 loc<-rbind(loc,tmp)
locid<-c(locid,(1:length(kp))[kp])
}
list(loc=loc,locid=locid)
}


#stopped here for now 


if(!exists("traps")){
traps<-get.traplocs(10)
 }
traplocs<-traps$loc
trap.id<-traps$locid

ntraps<-nrow(traplocs)
set.seed(2013)

N<-20
S.possible<- (1:nrow(pts))[in.pts==1]
S.id<-sample(S.possible,N,replace=TRUE)
S<- pts[S.id,]

 D<- outD[S.id,trap.id]

 Dtraps<-outD[trap.id,]

plot(pts)
points(traplocs,pch=20,cex=2,col="red")
points(S+.01,pch=20,col="blue")


alpha0<- -1.5
sigma<- 1
beta<- 1/(2*sigma*sigma)
K<-10

probcap<-plogis(alpha0)*exp(-beta*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}

Y<-Y[apply(Y,1,sum)>0,]



##### 
##### RUN SCR MODEL
#####
####



y<-Y
traplocs<-traplocs
nind<-nrow(y)
X<-traplocs
J<-nrow(X)
K<-K

Sg<-pts[in.pts==1,]
DD<- outD[in.pts==1,in.pts==1]


## To fit the model:
### 1. read in the function "intlik3" whic his given BELOW
### 2. Execute this function where Dtraps = ecological distance between each trap and each
####   state-space grid point
### 3. Change D=NULL if you want to fit Euclidean distance


frog<-nlm(intlik3,c(-2.5,2,log(4)),hessian=TRUE,y=y,K=K,X=traplocs,S=pts,D=Dtraps,inpoly=in.pts)

## Note Dtraps = ngridcellsinSS x ntraps  = Ecological distance
## inpoly = ngridcellsinSS x 1  = 1 if in polygon = 0 o/w

##### READ THIS IN FIRST

intlik3<-function(start=NULL,y=y,K=NULL,X=traplocs,S=NULL,D,inpoly){
if(is.null(K)) return("need sample size")
if(is.null(S)) return("beat it")
G<-S
nG<-nrow(G)

# note this assumes that input S is already subsetted

#inpoly<- rep(1,nrow(G))

###

# below this computes this integrated likelihood
#G<-G[inpoly==1,]
#nG<-nrow(G)
#PrS<-inpoly  # weight 0 or 1

G<-G[inpoly==1,]
nG<-nrow(G)
if(!is.null(D))
 D<-D[,inpoly==1]

if(is.null(D))
D<- e2dist(X,G)      # this is Dtraps if input


if(is.null(start)) start<-c(0,0,0)
alpha0<-start[1]
alpha1<-start[2]
n0<-exp(start[3])


probcap<- expit(alpha0)*exp(-alpha1*D*D)
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

out
}




























