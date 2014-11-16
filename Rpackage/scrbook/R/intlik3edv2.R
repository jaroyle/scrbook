intlik3edv2 <-
function(start=NULL,y=y,K=NULL,X=traplocs,S=NULL,D,inpoly){
if(is.null(K)) return("need sample size")
if(is.null(S)) return("beat it")
G<-S
nG<-nrow(G)
nind<-nrow(y)
J<-nrow(X)

 if(length(K)==1) K<- rep(K,nrow(X))
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


probcap<- plogis(alpha0)*exp(-alpha1*D*D)
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
