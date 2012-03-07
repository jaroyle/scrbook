intlik1 <-
function(parm,y=y,delta=.2,X=traplocs,ssbuffer=2){

Xl<-min(X[,1])-ssbuffer
Xu<-max(X[,1])+ ssbuffer
Yu<-max(X[,2])+ ssbuffer
Yl<-min(X[,2])-ssbuffer

xg<-seq(Xl+delta/2,Xu-delta/2,by=delta)
yg<-seq(Yl+delta/2,Yu-delta/2,by=delta)
npix<-length(xg)
###area<- (Xu-Xl)*(Yu-Yl)/((npix)*(npix))

G<-cbind(rep(xg,npix),sort(rep(yg,npix)))
nG<-nrow(G)
D<- e2dist(X,G)

alpha0<-parm[1]
alpha1<-parm[2]
probcap<- plogis(alpha0)*exp(-alpha1*D*D)
Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
n0<-sum(apply(y,1,sum)==0)
ymat<-y[apply(y,1,sum)>0,]
ymat<-rbind(ymat,rep(0,ncol(ymat)))
lik.marg<-rep(NA,nrow(ymat))
for(i in 1:nrow(ymat)){
Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),K,probcap[1:length(Pm)],log=TRUE))
lik.cond<- exp(colSums(Pm))
lik.marg[i]<- sum( lik.cond*(1/nG))
}
nv<-c(rep(1,length(lik.marg)-1),n0)
####-1*( lgamma( nrow(y)+1) -lgamma(n0+1)  +sum(nv*log(lik.marg)) )
-1*( sum(nv*log(lik.marg)) )

}
