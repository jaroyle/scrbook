### Source all of the utility functions at the end of this script before
### running any of the code. 
###
###

### We require 3 R libraries
library("shapefiles")
library("gdistance")
library("raster")

###
### following block of code creates a "patchy" looking covariate to use
### in our cost function. It uses a standard method for generating a correlated
### multivariate normal vector of length, in this case, 400 (one value for each
### pixel). One can use any correlation function here but we chose a standard
### exponential model with range parameter 0.5
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
par(mfrow=c(2,1))
plot(r)
####dev.off()
covariate.patchy<- z*sqrt(1.68) + .168   #approx. same scale as systematic covariate below
values(r)<-matrix(covariate.patchy,20,20,byrow=FALSE)
covariate.patchy<- r
class(covariate.patchy)

##
## build systematic covariate
## defined here as a trend from NW to SE
##
cost<-matrix(NA,nrow=20,ncol=20)
cost<-row(cost)+col(cost)
covariate.trend<- (cost-20)/10
values(r)<-matrix(covariate.trend,20,20,byrow=FALSE)
covariate.trend<-r
class(covariate.trend)




###
### generate Pr(encounter) for all pixels, for hypothetical traps in every other pixel
###  
### sigma = a bivariate normal home range parameter -- use different values to
### see how this looks
### theta0 = intercept of detection prob. model
### theta1 = coefficient on "distance" -- derived from "sigma"
### theta2 = coefficient on covariate used in cost function
###
theta2<- 1
theta0<- -2   
sigma<-.25
theta1<- 1/(2*sigma*sigma)
gx<- seq(.5 + delta/2, 4.5-delta/2,,20)
gy<- rev(gx)
gx<-sort(rep(gx,20))
gy<-rep(gy,20)
grid2<-cbind(gx,gy)
D<-e2dist(grid,grid2)
## Here is the stationary and isotropic detection model:
probcap<-plogis(theta0)*exp(-theta1*D*D)
## spatial.plot is a utility function (end of this script)
spatial.plot(grid2,probcap[1,])

###
### Set one of the cost functions -- either covariate.trend or covariate.patchy
###
cost<- exp(theta2*covariate.patchy)
#values(r)<-matrix(cost,20,20,byrow=FALSE)
r<-cost
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
costs1<-costDistance(tr1CorrC,grid,grid2)
outD<-as.matrix(costs1)
probcap<-plogis(theta0)*exp(-theta1*outD*outD)

###
### The following block of code will plot the "space usage" by 6 individuals with
### activity centers selected arbitrarily. You can plot this for any hypothetical
### individual if you wish
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




###
### R commands to carry-out the simulations. MUST SOURCE likelihood function first -- SEE BELOW
### This takes 1-2 days to run for nsims=100
###
###

nsims<-50

simout.low.N100<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=5,covariate=covariate.trend)
simout.low.N200<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=5,covariate=covariate.trend)
simout.reallylow.N100<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=3,covariate=covariate.trend)
simout.reallylow.N200<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=3,covariate=covariate.trend)
simout.high.N100<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=10,covariate=covariate.trend)
simout.high.N200<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=10,covariate=covariate.trend)

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

simout.low.N100.k<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=5,covariate=covariate.patchy)
simout.low.N200.k<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=5,covariate=covariate.patchy)
simout.reallylow.N100.k<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=3,covariate=covariate.patchy)
simout.reallylow.N200.k<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=3,covariate=covariate.patchy)
simout.high.N100.k<-sim.fn(N=100,nsim=nsims,theta0=-2,sigma=.5,K=10,covariate=covariate.patchy)
simout.high.N200.k<-sim.fn(N=200,nsim=nsims,theta0= -2, sigma=.5,K=10,covariate=covariate.patchy)

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


#
# to do a single run you have to define a few things and then execute the code within the "sim.fn" below
#
covariate<-covariate.patchy
N<-200
theta0<- -2
sigma<- .5
K<- 5

sim.fn<-function(N=200,nsim=100,theta0= -2, sigma=.5, K=5,covariate){
# input covariate as a RasterLayer
#
cl<-match.call()
set.seed(2013)
simout2<-simout1<-simout3<-matrix(NA,nrow=nsim,ncol=5)
theta1<- 1/(2*sigma*sigma)

r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
theta2<-1
cost<- exp(theta2*covariate)
#values(r)<-matrix(cost,20,20,byrow=FALSE)
#par(mfrow=c(1,1))
#plot(r)
r<-cost    

## use max = doesn't count moving through boundary pixel
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

xg<-seq(1,4,1)
yg<-4:1
pts<-cbind( sort(rep(xg,4)),rep(yg,4))
#costs1<-costDistance(tr1CorrC,pts)
#D<-as.matrix(costs1)

traplocs<-pts
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)


for(sim in 1:nsim){

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))
D<-costDistance(tr1CorrC,S,traplocs)
probcap<-plogis(theta0)*exp(-theta1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}
Y<-Y[apply(Y,1,sum)>0,]

# raster has to be defined for state-space and ssbuffer = 0 only
n0<- N-nrow(Y)
frog<-nlm(intlik3ed,c(theta0,theta1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,distmet="euclid",covariate=covariate,theta2=1)
simout1[sim,]<-c(frog$estimate,NA,nrow(Y))

frog<-nlm(intlik3ed,c(theta0,theta1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,distmet="ecol",covariate=covariate,theta2=1)
simout2[sim,]<-c(frog$estimate,NA,nrow(Y))

frog<-nlm(intlik3ed,c(theta0,theta1,log(n0),-.3),hessian=TRUE,y=Y,K=K,X=traplocs,distmet="ecol",covariate=covariate,theta2=NA)
simout3[sim,]<-c(frog$estimate,nrow(Y))
}
list(simout1=simout1,simout2=simout2,simout3=simout3,call=cl)

}

##
##
## UTILITY FUNCTIONS
## 
## PUT ALL OF THE OBJECTS BELOW INTO YOUR WORKSPACE
##
	
intlik3ed<-function(start=NULL,y=y,K=NULL,delta=.2,X=traplocs,distmet="ecol",covariate,theta2=NA){
if(is.null(K)) return("need sample size")
if(class(covariate)!="RasterLayer") {
 cat("make a raster out of this",fill=TRUE)
 return(NULL)

}
# do a check here that trap locations exist in same space as raster. 
# forthcoming

# build integration grid. This derives from the covariate raster
# i.e., potential values of s are the mid-point of each raster pixel
nc<-covariate@ncols
nr<-covariate@nrows
Xl<-covariate@extent@xmin
Xu<-covariate@extent@xmax
Yl<-covariate@extent@ymin
Yu<-covariate@extent@ymax
SSarea<- (Xu-Xl)*(Yu-Yl)
### ASSUMES SQUARE RASTER -- NEED TO GENERALIZE THIS
delta<- (Xu-Xl)/nc
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
if(is.na(theta2))
theta2<-exp(start[4])

cost<- exp(theta2*covariate)
tr1<-transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
D<-costDistance(tr1CorrC,X,G)
}

if(is.null(start)) start<-c(0,0,0,0)
theta0<-start[1]
theta1<-start[2]
n0<-exp(start[3])


probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
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















image.scale <-
function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", 
    "ranges"))
{
    # sort out the location
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05*dx/2	# default x to right of image
    else if (is.list(x)) {
        if (length(x$x) == 2) 
          size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
          size <- 0.618*dy/n	# default size, golden ratio
          y <- my + 0.618*dy/2	# default y to give centred scale
        } else size <- (y-my)*2/n
    if (length(size)==1)
        size <- rep(size, 2)	# default square boxes
    if (is.null(y))
        y <- my + n*size[2]/2
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
        col = rev(col), xpd = TRUE)
    # sort out the labels
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
        ifelse(size[1]>0, 0, 1), xpd = TRUE) 
}



spatial.plot <-
function(x,y,add=TRUE,cx=1){
 nc<-as.numeric(cut(y,20))
if(!add) plot(x,pch=" ")
 points(x,pch=20,col=terrain.colors(20)[nc],cex=cx)
image.scale(y,col=terrain.colors(20))

}




e2dist <- function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}






smy.fn <-
function(x){
c(mean(x),sqrt(var(x)),quantile(x,c(0.025,0.50,0.975)))
}












