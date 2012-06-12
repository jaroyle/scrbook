intlik3Poisson <-
function(start=NULL,y=y,K=NULL,delta=.3,X=traplocs,ssbuffer=2){

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
D<- e2dist(X,G)  

if(is.null(start)) start<-c(0,0,0)
alpha0<-start[1]
alpha1<-start[2]
Dens<-exp(start[3])


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
nv<-c(rep(1,length(lik.marg)-1),1)
atheta<- 1-lik.marg[nrow(ymat)]
nind<-nrow(ymat)-1
part1<- nind*log(Dens) - Dens*atheta
part2<- sum(nv[1:nind]*log(lik.marg[1:nind]))
out<-  -1*(part1+ part2)
attr(out,"SSarea")<- SSarea
out
}
