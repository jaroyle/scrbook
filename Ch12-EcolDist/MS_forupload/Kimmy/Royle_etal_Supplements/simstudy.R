### We require 3 R libraries
library("shapefiles")
library("gdistance")
library("raster")

###Before we get started, we'll run some utility functions 
## UTILITY FUNCTIONS 
## 
## PUT ALL OF THE OBJECTS BELOW INTO YOUR WORKSPACE
##

####This function creates a color scale used with the function image
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

####This function makes a nice 3D plot
spatial.plot <-
function(x,y,add=TRUE,cx=1){
 nc<-as.numeric(cut(y,20))
if(!add) plot(x,pch=" ")
 points(x,pch=20,col=terrain.colors(20)[nc],cex=cx)
image.scale(y,col=terrain.colors(20))
}

####This function calculates geographic distances among two sets of locations 
####e.g., will work for traplocs and integration grid instead of just among traplocs 
e2dist <- function (x, y) {
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
    }
####This function calculates summary statistics for an input x
smy.fn <- function(x){
    c(mean(x),sqrt(var(x)),quantile(x,c(0.025,0.50,0.975)))
    }

##########################################################################
####CREATE COVARIATES
### Now we create the covariates that we will use in the simulation.
### The following block of code creates a "patchy" looking covariate to use
### in our cost function. It uses a standard method for generating a correlated
### multivariate normal vector of length, in this case, 400 (one value for each
### pixel). One can use any correlation function here but we chose a standard
### exponential model with range parameter 0.5.
par(mfrow=c(1,1))
#setting the seed ensures that your results will match ours 
#delete this line for stochastic results
set.seed(12)  

n.pix= 20 #number of pixels in raster
r<-raster(nrows=n.pix,ncols=n.pix) #create an empty raster
projection(r)<- "+proj=utm +zone=12 +datum=WGS84" #set the projection to be UTMs
xmin<- ymin<- .5
xmax<- ymax<- 4.5
extent(r)<-c(xmin,xmax,ymin,ymax) #set the extent of the raster
delta<- (xmax-xmin)/n.pix  #find the distance between points 
gx<- seq(xmin + delta/2, xmax-delta/2,,n.pix) #create a sequence for the x locations
gy<- rev(gx) #create a sequence for the y locations in reverse order
gx<-sort(rep(gx,20)) #sort and repeat the sequence of x locations for all columns
gy<-rep(gy,20) #repeat the y locations for all rows
grid<-cbind(gx,gy) #bind the x and y components of the locations
Dmat<-as.matrix(dist(grid)) #calculate the Euclidean distances among the potential activity centers 
V<-exp(-Dmat/.5)  #Create the correlation function: standard exponential with range paramater 0.5
z<-t(chol(V))%*%rnorm(400) #Create the correlated variable z, with some random noise
z<- (z-mean(z))/sqrt(var(as.vector(z))) #center the correlated variable
values(r)<-matrix(z,20,20,byrow=FALSE)  #assign z to the raster r
#plot the covariate
par(mfrow=c(2,1)) 
hist(z)
plot(r)
points(grid)
#dev.off() #removes figure
covariate.patchy<- z*sqrt(1.68) + .168   #approx. same scale as systematic covariate below
values(r)<-matrix(covariate.patchy,20,20,byrow=FALSE) #assign the scaled covariate values to r
covariate.patchy<- r  #create a new raster called covariate.patchy
#class(covariate.patchy) #check to confirm that covariate.patchy is a RasterLayer

## SYSTEMATIC COVARIATE
## It's much simpler to build a systematic (or trend) covariate.
## Here we defined as a trend from NW to SE
## 
cost<-matrix(NA,nrow=n.pix,ncol=n.pix)
cost<-row(cost)+col(cost)
covariate.trend<- (cost-20)/10 #define cost values
values(r)<-matrix(covariate.trend,n.pix,n.pix,byrow=FALSE) #assign costs to raster r
covariate.trend<-r
class(covariate.trend)
plot(covariate.trend)


####SIMULATION FUNCTION
####Next we define a few conditions and then the simulation function, 
####where we simulate activity centers and trap locations, calculate
####the probability of capture, simulate captures, and then fit the 
####simulated data.  To do this a single time, first run the utility
####functions below, run the following conditions, and then run only 
####the code inside of the function.

N<-200  #the number of individuals=activity centers
alpha0<- -2  
sigma<- .5  #scale parameter for the distance distribution
K<- 5
nsim=1
covariate<-covariate.patchy

sim.fn<-function(N=200,nsim=100,alpha0= -2, sigma=.5, K=5,covariate){
#N= true number of individuals
#nsim = the number of simulations you'd like to run
#alpha0 is the resistance covariate
#sigma is the scale parameter
#K is the number of occasions
#covariate refers to the RasterLayer object that is the covariate

cl<-match.call() #stores all the arguments in the function call for later reference
#Create matrices to store the output
simout2<-simout1<-simout3<-matrix(NA,nrow=nsim,ncol=5)
alpha1<- 1/(2*sigma*sigma) #calculate alpha based on the sigma (scale parameter)

#Create a raster object
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84" #assign the projection
extent(r)<-c(.5,4.5,.5,4.5)  #assign the extent of the raster
alpha2<-1 #assign the resistance coefficient
cost<- exp(alpha2*covariate) #calculate the resistance surface 
#(i.e. the cost per pixel for each pixel of the covariate

#par(mfrow=c(1,1))
#plot(r)
r<-cost #rename the cost raster r
#identify neighbors and calculate distances as conductances (1/average distance between 2 pixels)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)
#correct the diagonal distances 
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE) 

#Create the trap locations composed of locations xg, yg
xg<-seq(1,4,1)
yg<-4:1
traplocs<-cbind( sort(rep(xg,4)),rep(yg,4))
#graph the traplocations
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)


for(sim in 1:nsim){

#Create the activity centers, distributed uniformly in the state space= extent
#Note that the raster defines the state space because least cost distance 
#can not be calculated where there are no covariate values.

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))

#calculate the least cost distances between activity centers and traplocations 

D<- costDistance(tr1CorrC,S,traplocs)  #N x ntraps matrix

#calculate the probability of capture based on the weighted distance distribution

probcap<-plogis(alpha0)*exp(-alpha1*D*D) #N X ntraps matrix

# now generate the number of encounters of every individual in every trap

Y<-matrix(NA,nrow=N,ncol=ntraps) #N X ntraps matrix
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,]) #simulate # of encounters from K occasions
}

#select only those individuals that were counted at least once for the capture history

Y<-Y[apply(Y,1,sum)>0,] #matrix of N actually captured X ntraps
n0<- N-nrow(Y)  #number of zero-encounter histories

#calculate the likelihood based solely on Euclidean distance
frog<-optim(c(alpha0,alpha1,log(n0)),intlik3ed,hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="euclid",covariate=covariate,alpha2=1)
# if using nlm and sometimes with optim, warnings can be produced. This
# is due to parameters outside of parameter space or small number 
# arithmetic. 
# ignore these warnings.

simout1[sim,]<-c(frog$par,NA,nrow(Y))

#calculate the likelihood based on ecological distance with cost coefficient fixed
frog<-optim(c(alpha0,alpha1,log(n0)),intlik3ed,hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="ecol",covariate=covariate,alpha2=1)
simout2[sim,]<-c(frog$par,NA,nrow(Y))

#calculate the likelihood based on ecological distance with cost coefficient estimated
frog<-optim(c(alpha0,alpha1,log(n0),-.3),intlik3ed,hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="ecol",covariate=covariate,alpha2=NA)
simout3[sim,]<-c(frog$par,nrow(Y))
}
#output results of simulations, plus the inputs to the function (cl defined above) 
list(simout1=simout1,simout2=simout2,simout3=simout3,call=cl)

}  #end of the simulation function



#####The next groups of code were used to specify and summarize the simulations reported in the paper.
###
### R commands to carry-out the simulations. MUST SOURCE likelihood function first -- SEE BELOW
### This takes 1-2 days to run for nsims=100
###
###

nsims<-50
###
####Systematic covariate 
###
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
#### Patchy covariate
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
