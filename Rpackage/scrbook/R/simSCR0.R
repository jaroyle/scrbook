simSCR0 <-
function(N=100,K=20,alpha0=-2.5,sigma=.5, discard0=TRUE,array3d=FALSE,rnd=2013){
set.seed(rnd)

# make trapping grid. Normally you would provide a 2-dimensional matrix
#   of trap coordinates and read it in like this:
#    tralocs<-read.csv("traplocs.csv") or similar
traplocs<- cbind(sort(rep(1:5,5)),rep(1:5,5))
Dmat<-e2dist(traplocs,traplocs)
ntraps<-nrow(traplocs)
plot(traplocs)

# define state-space of point process. (i.e., where animals live).
# Here "delta" just adds
# a fixed buffer to the outer extent of the traps.
buffer<-2
Xl<-min(traplocs[,1] - buffer)
Xu<-max(traplocs[,1] + buffer)
Yl<-min(traplocs[,2] - buffer)
Yu<-max(traplocs[,2] + buffer)

sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
S<-cbind(sx,sy)

# how far is each individual from each trap?
D<- e2dist(S,traplocs)

#alpha0<- -2.5
#sigma<- 0.5
alpha1<- 1/(2*sigma*sigma)

#cloglog.probcap<- alpha0  - alpha1*D*D
#probcap<- 1-exp(-exp(cloglog.probcap))
# this is logit model here:
#probcap<- expit(-2.5 - alpha1*D)

probcap<-plogis(alpha0)*exp(-alpha1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:nrow(Y)){
 Y[i,]<-rbinom(ntraps,K,probcap[i,])
}

## NOTE
## Y is a matrix of encounter frequencies of EACH individual in EACH trap
## As simulated here it includes the "all 0" observations.  We want
## to delete those to mimic real data.
if(discard0){
  totalcaps<-apply(Y,1,sum)
  Y<-Y[totalcaps>0,]
}

dimnames(Y)<-list(1:nrow(Y),paste("trap",1:ncol(Y),sep=""))

if(array3d){

## Here we demonstrate how to simulate the full 3-dimensional
## encounter history
# now generate the encounters of every individual in every trap

Y<-array(NA,dim=c(N,ntraps, K))
for(i in 1:nrow(Y)){
for(j in 1:ntraps){
 Y[i,j,1:K]<-rbinom(K,1,probcap[i,j])
}
}

if(discard0){
## create nind x ntraps array
  Y2d<- apply(Y,c(1,2),sum)
## which individuals were captured?
  ncaps<-apply(Y2d,1,sum)
## keep those ones that were captured
  Y<-Y[ncaps>0,,]
}

}


list(Y=Y,traplocs=traplocs,xlim=c(Xl,Xu),ylim=c(Yl,Yu),N=N,alpha0=alpha0,alpha1=alpha1,sigma=sigma,K=K)
}
