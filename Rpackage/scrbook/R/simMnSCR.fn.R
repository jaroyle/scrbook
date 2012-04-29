simMnSCR.fn <-
function(parms,K=7,ssbuff=2){

traplocs<- cbind(sort(rep(1:5,5)),rep(1:5,5))
Dmat<-e2dist(traplocs,traplocs)
ntraps<-nrow(traplocs)
plot(traplocs)

delta<-ssbuff
Xl<-min(traplocs[,1] - delta)
Xu<-max(traplocs[,1] + delta)
Yl<-min(traplocs[,2] - delta)
Yu<-max(traplocs[,2] + delta)
N<-parms$N
sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
S<-cbind(sx,sy) 

# how far is each individual from each trap?
D<- e2dist(S,traplocs)

sigma<-parms$sigma
alpha0<-parms$alpha0
###alpha1<-parms$alpha1
alpha1<- 1/(2*sigma*sigma)
alpha2<-parms$alpha2

Ycat<-matrix(NA,nrow=N,ncol=K)
Xlag<-matrix(0,nrow=N,ncol=K+1)
Xlag[,1]<-rep(0,N)
for(i in 1:N){
for(k in 1:K){
lp<- alpha0 + alpha2*Xlag[i,k] - alpha1*D[i,]*D[i,]
cp<- exp(c(lp,0))
cp<- cp/sum(cp)
Ycat[i,k]<- sample(1:(ntraps+1),1,prob=cp)
if(Ycat[i,k] <=ntraps)
 Xlag[i,(k+1):ncol(Xlag)]<-1
}
}
captured<-apply(Ycat<=ntraps,1,sum)
captured<-captured>0
###
### INPUTS
####
Ycat<-Ycat[captured,]
Xlag<-Xlag[captured,]
reencounter=Xlag[,1:K]
X<-traplocs
K<-ncol(reencounter)
S1<-S[captured,]
S2<-S[!captured,]
S<-rbind(S1,S2)
## function returns actual activity centers for evaluation
list(Ycat=Ycat,X=traplocs,reencounter=reencounter,
X=X,ssbuff=ssbuff,S=S1,xlim=c(Xl,Xu),ylim=c(Yl,Yu),K=K)
}
