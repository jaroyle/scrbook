

### Estimating pixel value related to a covariate


library("shapefiles")
library("gdistance")
library("raster")

###png("raster_krige.png",width=5,height=5, units="in", res=400)

par(mfrow=c(1,1))
set.seed(12)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
delta<- (4.5-.5)/20
gx<- seq(.5 + delta/2, 4.5-delta/2,,20)
gy<- rev(gx)
gx<-sort(rep(gx,20))
gy<-rep(gy,20)
grid<-cbind(gx,gy)
Dmat<-as.matrix(dist(grid))
V<-exp(-Dmat/.5)
z<-t(chol(V))%*%rnorm(400)
z<- (z-mean(z))/sqrt(var(as.vector(z)))
#values(r)<-rot(rot(rot(matrix(z,20,20,byrow=TRUE))))
#values(r)<-(matrix(z,20,20,byrow=FALSE))
values(r)<-matrix(z,20,20,byrow=FALSE)
spatial.plot(grid,z)
plot(r)
####dev.off()
cost<-r
covariate.patchy<- z*sqrt(1.68) + .168



# systematic covariate
cost<-matrix(NA,nrow=20,ncol=20)
cost<-row(cost)+col(cost)
covariate.trend<- (cost-20)/10




###
### generate home range sizes...D = dist from individual to trap
###  D<- dist(grid)  
#
sigma<-.25
theta<- 1/(2*sigma*sigma)
gx<- seq(.5 + delta/2, 4.5-delta/2,,20)
gy<- rev(gx)
gx<-sort(rep(gx,20))
gy<-rep(gy,20)
grid2<-cbind(gx,gy)
D<-e2dist(grid,grid2)
probcap<-plogis(alpha0)*exp(-theta*D*D)
spatial.plot(grid2,probcap[1,])

cost<- exp(beta*covariate.trend)
values(r)<-matrix(cost,20,20,byrow=FALSE)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
costs1<-costDistance(tr1CorrC,grid,grid2)
outD<-as.matrix(costs1)
probcap<-plogis(alpha0)*exp(-theta*outD*outD)

###png("home_ranges.png",width=5,height=8, units="in", res=400)
par(mfrow=c(3,2),mar=c(2,2,2,2))
#spatial.plot(grid2,probcap[63,])
#spatial.plot(grid2,probcap[76,])
values(r)<-matrix(probcap[63,],20,20,byrow=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[63,],nrow=1),pch=20,cex=2)
values(r)<-matrix(probcap[76,],20,20,byrow=FALSE)
#plot(r,axes=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[76,],nrow=1),pch=20,cex=2)


#spatial.plot(grid2,probcap[183,])
#spatial.plot(grid2,probcap[196,])
#spatial.plot(grid2,probcap[363,])
#spatial.plot(grid2,probcap[376,])
values(r)<-matrix(probcap[183,],20,20,byrow=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[183,],nrow=1),pch=20,cex=2)

values(r)<-matrix(probcap[196,],20,20,byrow=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[196,],nrow=1),pch=20,cex=2)

values(r)<-matrix(probcap[363,],20,20,byrow=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[363,],nrow=1),pch=20,cex=2)

values(r)<-matrix(probcap[376,],20,20,byrow=FALSE)
image(r,col=terrain.colors(10),axes=FALSE)
points(matrix(grid[376,],nrow=1),pch=20,cex=2)

#####dev.off()

 
smy.fn <-
function(x){
c(mean(x),sqrt(var(x)),quantile(x,c(0.025,0.50,0.975)))
}

###
### R commands to carry-out the simulations. MUST SOURCE likelihood function first -- SEE BELOW
###

nsims<-2

simout.low.N100<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=5,covariate=covariate.trend)
simout.low.N200<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=5,covariate=covariate.trend)
simout.reallylow.N100<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=3,covariate=covariate.trend)
simout.reallylow.N200<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=3,covariate=covariate.trend)
simout.high.N100<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=10,covariate=covariate.trend)
simout.high.N200<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=10,covariate=covariate.trend)

###
### R commands to summarize the simulation output
###
mat<-matrix(NA,nrow=9,ncol=10)
for(i in 1:3){
mat[i,1:5]<- smy.fn(exp(simout.reallylow.N100[[i]][,3]) + simout.reallylow.N100[[i]][,5])
mat[i,6:10]<- smy.fn(exp(simout.reallylow.N200[[i]][,3]) + simout.reallylow.N200[[i]][,5])
mat[3+i,1:5]<- smy.fn(exp(simout.low.N100[[i]][,3]) + simout.low.N100[[i]][,5])
mat[3+i,6:10]<- smy.fn(exp(simout.low.N200[[i]][,3]) + simout.low.N200[[i]][,5])
mat[6+i,1:5]<- smy.fn(exp(simout.high.N100[[i]][,3]) + simout.high.N100[[i]][,5])
mat[6+i,6:10]<- smy.fn(exp(simout.high.N200[[i]][,3]) + simout.high.N200[[i]][,5])
}

###
### 
### Using the "patchy" cost function execute these commands
###
###

simout.low.N100.k<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=5,covariate=covariate.patchy)
simout.low.N200.k<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=5,covariate=covariate.patchy)
simout.reallylow.N100.k<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=3,covariate=covariate.patchy)
simout.reallylow.N200.k<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=3,covariate=covariate.patchy)
simout.high.N100.k<-sim.fn(N=100,nsim=nsims,alpha0=-2,sigma=.5,K=10,covariate=covariate.patchy)
simout.high.N200.k<-sim.fn(N=200,nsim=nsims,alpha0= -2, sigma=.5,K=10,covariate=covariate.patchy)

###
### R commands to summarize the simulation output
###
mat.k<-matrix(NA,nrow=9,ncol=10)
for(i in 1:3){
mat.k[i,1:5]<- smy.fn(exp(simout.reallylow.N100.k[[i]][,3]) + simout.reallylow.N100.k[[i]][,5])
mat.k[i,6:10]<- smy.fn(exp(simout.reallylow.N200.k[[i]][,3]) + simout.reallylow.N200.k[[i]][,5])
mat.k[3+i,1:5]<- smy.fn(exp(simout.low.N100.k[[i]][,3]) + simout.low.N100.k[[i]][,5])
mat.k[3+i,6:10]<- smy.fn(exp(simout.low.N200.k[[i]][,3]) + simout.low.N200.k[[i]][,5])
mat.k[6+i,1:5]<- smy.fn(exp(simout.high.N100.k[[i]][,3]) + simout.high.N100.k[[i]][,5])
mat.k[6+i,6:10]<- smy.fn(exp(simout.high.N200.k[[i]][,3]) + simout.high.N200.k[[i]][,5])
}




sim.fn<-function(N=200,nsim=100,alpha0= -2, sigma=.5, K=5,covariate){
cl<-match.call()
set.seed(2013)
simout2<-simout1<-simout3<-matrix(NA,nrow=nsim,ncol=5)

r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
beta<-1
cost<- exp(beta*covariate)
values(r)<-matrix(cost,20,20,byrow=FALSE)
par(mfrow=c(1,1))
plot(r)


## use max = doesn't count moving through boundary pixel
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

xg<-seq(1,4,1)
yg<-4:1
pts<-cbind( sort(rep(xg,4)),rep(yg,4))
#costs1<-costDistance(tr1CorrC,pts)
#outD<-as.matrix(costs1)
traplocs<-pts
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)


for(sim in 1:nsim){

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))
D<-costDistance(tr1CorrC,S,traplocs)
theta<- 1/(2*sigma*sigma)
probcap<-plogis(alpha0)*exp(-theta*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}
Y<-Y[apply(Y,1,sum)>0,]

# raster has to be defined for state-space and ssbuffer = 0 only
n0<- N-nrow(Y)
frog<-nlm(intlik3ed,c(alpha0,beta,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,ssbuffer=0.5,distmet="euclid",covariate=covariate,beta=1)
simout1[sim,]<-c(frog$estimate,NA,nrow(Y))

frog<-nlm(intlik3ed,c(alpha0,beta,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,ssbuffer=0.5,distmet="ecol",covariate=covariate,beta=1)
simout2[sim,]<-c(frog$estimate,NA,nrow(Y))

frog<-nlm(intlik3ed,c(alpha0,beta,log(n0),-.3),hessian=TRUE,y=Y,K=K,X=traplocs,ssbuffer=0.5,distmet="ecol",covariate=covariate,beta=NA)
simout3[sim,]<-c(frog$estimate,nrow(Y))
}
list(simout1=simout1,simout2=simout2,simout3=simout3,call=cl)

}

## 
## the object of class "scr" will need to have y, X, covariate information
##

intlik3ed<-function(start=NULL,y=y,K=NULL,delta=.1,X=traplocs,ssbuffer=2,distmet="ecol",covariate,beta=NA){

Xl<-min(X[,1]) -ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])- ssbuffer
SSarea<- (Xu-Xl)*(Yu-Yl)
if(is.null(K)) return("need sample size")
#delta<- (Xu-Xl)/npix
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
area<- (Xu-Xl)*(Yu-Yl)/((npix.x)*(npix.y))
G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))
nG<-nrow(G)
if(distmet=="euclid")
D<- e2dist(X,G)  
if(distmet=="ecol"){
if(is.na(beta))
beta<-exp(start[4])

r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
cost<- exp(beta*covariate)
values(r)<-matrix(cost,20,20,byrow=FALSE)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D<-costDistance(tr1CorrC,X,G)
}

if(is.null(start)) start<-c(0,0,0,0)
alpha0<-start[1]
alpha1<-start[2]
n0<-exp(start[3])


probcap<- (exp(alpha0)/(1+exp(alpha0)))*exp(-alpha1*D*D)
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
attr(out,"SSarea")<- SSarea
out
}



































