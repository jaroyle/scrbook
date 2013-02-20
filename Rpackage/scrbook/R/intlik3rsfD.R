intlik3rsfD <-
function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL,stel=NULL){
#
# this version of the code handles a covariate on log(Density). This is starting value 5
#
# start = vector of length 5 = starting values
# y = nind x ntraps encounter matrix
# K = how many samples?
# X = trap locations
# ztrap = covariate value at trap locations
# zall = all covariate values for all nG pixels
# ntel = nguys x nG matrix of telemetry fixes in each nG pixels
# stel = home range center of telemetered individuals, IF you wish to estimate it. Not necessary

nG<-nrow(G)
D<- e2dist(X,G)

alpha0<-start[1]
sigma<- exp(start[2])
alpha2<- start[3]
n0<-    exp(start[4])
beta<- start[5]
a0<- 1
if(!is.null(zall)){
 psi<- exp(beta*zall)
 psi<-psi/sum(psi)
}
else{
psi<-rep(1/nG,nG)
}
if(!is.null(y)){
loglam<-   alpha0  -(1/(2*sigma*sigma))*D*D + alpha2*ztrap  # ztrap recycled over nG


probcap<- 1-exp(-exp(loglam))
#probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
ymat<-y
ymat<-rbind(y,rep(0,ncol(y)))
lik.marg<-rep(NA,nrow(ymat))
for(i in 1:nrow(ymat)){
Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
lik.cond<- exp(colSums(Pm))
lik.marg[i]<- sum( lik.cond*psi )
}
nv<-c(rep(1,length(lik.marg)-1),n0)
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
part2<- sum(nv*log(lik.marg))
out<-  -1*(part1+ part2)
}
else{
out<-0
}

if(!is.null(ntel) & !is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(stel,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- lam/denom  # each column is the probs for a guy at column [j]

tel.loglik<-  -1*sum(  ntel*log(probs) )

out<- out  + tel.loglik
}

if(!is.null(ntel) & is.null(stel) ){

# this is a tough calculation here
D2<-  e2dist(G,G)^2
# lam is now nG x nG!
lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall))  # recycle zall over all ntel guys
denom<-rowSums(lam)
probs<- t(lam/denom)  # each column is the probs for a guy at column [j]
temp<-exp(ntel%*%log(probs))  # Ntel x nG matrix

marg<- as.vector(rowSums(  temp*psi ))


tel.loglik<- -1*sum(log(marg))

out<- out  + tel.loglik
}

out
}
