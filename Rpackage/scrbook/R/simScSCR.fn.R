simScSCR.fn <-
function(parms,K=7,ssbuff=2){
# alpha2 = behavior effect not yet implemented
N<-parms$N
alpha0<-parms$alpha0
sigma<-parms$sigma
alpha2<- parms$alpha2
alpha1<- 1/(2*sigma*sigma)

traplocs<- cbind(sort(rep(1:10,10)),rep(1:10,10))
Dmat<-e2dist(traplocs,traplocs)
ntraps<-nrow(traplocs)
#plot(traplocs)

delta<-ssbuff
Xl<-min(traplocs[,1] - delta)
Xu<-max(traplocs[,1] + delta)
Yl<-min(traplocs[,2] - delta)
Yu<-max(traplocs[,2] + delta)



sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
S<-cbind(sx,sy) 

# how far is each individual from each trap?
D<- e2dist(S,traplocs)

Ycat<-matrix(NA,nrow=N,ncol=K)
Xlag<-matrix(0,nrow=N,ncol=K+1)
Xlag[,1]<-rep(0,N)

for(k in 1:K){
random<-sample(1:N,N,replace=FALSE)
avail<-rep(1,nrow(traplocs)+1)

for(i in random){
lp<- alpha0 + alpha2*Xlag[i,k] - alpha1*D[i,]*D[i,]

cp<- exp(c(lp,0))*avail
cp<- cp/sum(cp)
Ycat[i,k]<- sample(1:(ntraps+1),1,prob=cp)
if(Ycat[i,k] <=ntraps){
 Xlag[i,(k+1):ncol(Xlag)]<-1
 avail[Ycat[i,k]]<-0   # remove that trap
}
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

S1<-S[captured,]
S2<-S[!captured,]
S<-rbind(S1,S2)
list(Ycat=Ycat,reencounter=reencounter,S=S,X=X,ssbuff=ssbuff,K=K,xlim=c(Xl,Xu),ylim=c(Yl,Yu))
}
