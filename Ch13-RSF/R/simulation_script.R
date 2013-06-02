RSFsim<-function(){

library(scrbook)

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
simout0<-simout1<-simout2<-simout3<-matrix(NA,nrow=nsim,ncol=6)


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

tmp1<-optim(c(-3,log(3),1,0,1),intlik3rsf,y=y,K=K,X=X,ztrap=z[raster.point],G=gr)
simout0[sim,]<-c(nrow(y),tmp1$estimate)

# use telemetry s random
tmp2<-optim(c(-3,log(3),1,0,1),intlik3rsf,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))
simout1[sim,]<-c(nrow(y),tmp2$estimate)

# use mean "s" instead of estimating it 
tmp3<-optim(c(-3,log(3),1,0,1),intlik3rsf,y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z),stel=sbar)
simout2[sim,]<-c(nrow(y),tmp3$estimate)

# no SCR data, s is random
tmp4<-optim(c(-3,log(3),1,0,1),intlik3rsf,y=NULL,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))
simout3[sim,]<-c(NA,tmp4$estimate)
}


return(simout0,simout1,simout2,simout3)

}

