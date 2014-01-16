intlik1 <-
function(parm,y=y,delta=.2,X=traplocs,ssbuffer=2){

Xl<-min(X[,1])-ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])-ssbuffer

## These commands set up the integration grid
xg<-seq(Xl+delta/2,Xu-delta/2,by=delta)
yg<-seq(Yl+delta/2,Yu-delta/2,by=delta)
npix<-length(xg)
G<-cbind(rep(xg,npix),sort(rep(yg,npix)))
nG<-nrow(G)

D<- e2dist(X,G)

alpha0<-parm[1]
## alpha1 should probably be restricted to be positive here using exp()
alpha1<-parm[2]
probcap<- plogis(alpha0)*exp(-alpha1*D*D)


Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
# n0 here is the number of individuals with all-zero encounter histories
# the next 3 lines replace ALL of the "all-0" encounter histories with a
# single one, so that the calculation is only done once.
n0<-sum(apply(y,1,sum)==0)
ymat<-y[apply(y,1,sum)>0,]
ymat<-rbind(ymat,rep(0,ncol(ymat)))
lik.marg<-rep(NA,nrow(ymat))
# for each encounter history, compute the likelihood for each value of s
# then average over all possible values
for(i in 1:nrow(ymat)){
   Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),K,probcap[1:length(Pm)],log=TRUE))
   lik.cond<- exp(colSums(Pm))
   lik.marg[i]<- sum( lik.cond*(1/nG))
}
# nv is a vector of the frequency of each encounter history. Here it is
#   a vector of 1's with n0 being the number of all-0 histories.
nv<-c(rep(1,length(lik.marg)-1),n0)
-1*( sum(nv*log(lik.marg)) )

}
