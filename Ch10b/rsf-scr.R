set.seed(1234)

gr<-make.grid(minx=1,maxx=40,miny=1,maxy=40,nx=40,ny=40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))%*%rnorm(1600)
spatial.plot(gr,z)
Ntel<-10
sid<-sample(1:1600,Ntel,replace=TRUE)
s<-gr[sid,]

choose.s<-function(n){
x<-locator(n)
u<-x$x
v<-x$y
dd<- e2dist(cbind(u,v),gr)
id<-rep(NA,n)
for(i in 1:nrow(dd)){
xx<-dd[i,]
kp<-xx==min(xx)
id[i]<- (1:1600)[ kp ]
}

list("id"=id,"s"=gr[id,])

}
stest<-choose.s(8)
s<-stest$s
sid<-stest$id
Ntel<-nrow(s)

nsim<-100

simout0<-matrix(NA,nrow=nsim,ncol=5)

for(sim in 1:nsim){

# we combine distance and a covariate into a model for space usage
# individuals use each pixel in proportion to exp( distance + covariate)
# as in the following expression. We simulate counts here, as if we were
# recording telemetry location frequencies:
sigma<- 2
beta<- 1
n<-matrix(NA,nrow=Ntel,ncol=1600)
Nfixes<-50

i<-1  # do this for each telemetered guy to simulate a number of fixes.
      # note that n = 0 for most of the landscape
par(mfrow=c(3,3))
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
   d<- Dmat[sid[i],]
   lam<- exp(1 - (1/(2*sigma*sigma))*d*d + beta* z) 
n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
#   n[i,]<- rpois(1600,lam)
   par(mar=c(3,3,3,6))
   #spatial.plot(gr,lam)
lammat[i,]<-lam
   img<- matrix(lam,nrow=40,ncol=40,byrow=FALSE)
   image(1:40,1:40,rot(img),col=terrain.colors(10))
}
if(1==2){
png("habitat.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1),mar=c(3,3,3,6))
image(1:40,1:40,rot(matrix(z,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
image.scale(z,col=terrain.colors(10))
points(s,pch=20)
dev.off()


png("homeranges8.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1))
tot<- apply(lammat,1,sum)
lammat<-lammat/tot
lamtot<-apply(lammat,2,sum)
image(1:40,1:40,rot(matrix(lamtot,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
points(s,pch=20)
dev.off()
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
N<- 100
Sid<- sample(1:1600,N,replace=TRUE)
S<-gr[Sid,]

D<- e2dist(S,X)  ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<-   -2  -(1/(2*sigma*sigma))*D*D + beta*Zmat
p<- 1-exp(-exp(loglam))

## Now simulate SCR data

K<- 10
y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:N){
y[i,]<- rbinom(ntraps,K,p[i,])
}

cap<-apply(y,1,sum)>0

y<-y[cap,]

#####nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr)
tmp<-nlm(intlik3rsf.v2,c(-3,log(3),1,0,1),y=y,K=K,X=X,ztrap=z[raster.point],G=gr)
#tmp<-nlm(intlik3rsf.v2,c(-3,log(3),1,0,1),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))


simout0[sim,]<-tmp$estimate

}




# use estimated activity centers, with multinomial likelihood
nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[raster.point],G=gr,ntel=n,zall=as.vector(z))

intlik3rsf <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL){
if(is.null(K)) return("need sample size")
# Z = trap specific covariate
# zall = all covariate values for all nG pixels
# ntel = nguys x nG matrix of telemetry fixes in each nG pixels

nG<-nrow(G)
D<- e2dist(X,G)

if(is.null(start)) start<-c(0,0,0,0)
alpha0<-start[1]
sigma<- exp(start[2])
alpha2<- start[3]
n0<-    exp(start[4])

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
lik.marg[i]<- sum( lik.cond*(1/nG) )
}
nv<-c(rep(1,length(lik.marg)-1),n0)
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
part2<- sum(nv*log(lik.marg))
out<-  -1*(part1+ part2)

if(!is.null(ntel)){
sbar<- (ntel%*%G)/as.vector(ntel%*%rep(1,nG))
D2<-  e2dist(G,sbar)^2
lam<- exp(0 - (1/(2*sigma*sigma))*D2+ alpha2*zall)  # recycle zall over all ntel guys 
# lam is nG x ntel
denom<- as.vector(t(lam)%*%rep(1,nG))
probs<- t(lam)/denom  # now this is ntel x nG
tel.loglik<- sum(ntel*log(probs))

out<- out - tel.loglik
}





out
}




intlik3rsf.v2 <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL){
if(is.null(K)) return("need sample size")

# this version computes the marginal for the telemetry data too

# Z = trap specific covariate
# zall = all covariate values for all nG pixels
# ntel = nguys x nG matrix of telemetry fixes in each nG pixels

nG<-nrow(G)
D<- e2dist(X,G)

alpha0<-start[1]
sigma<- exp(start[2])
alpha2<- start[3]
n0<-    exp(start[4])
a0<- start[5]

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
lik.marg[i]<- sum( lik.cond*(1/nG) )
}
nv<-c(rep(1,length(lik.marg)-1),n0)
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
part2<- sum(nv*log(lik.marg))
out<-  -1*(part1+ part2)

if(!is.null(ntel)){

# this is a tough calculation here
D2<-  e2dist(G,G)^2
# lam is now nG x nG!
## should estimate the nuuisance parameter
lam<- exp(a0 - (1/(2*sigma*sigma))*D2+ alpha2*zall)  # recycle zall over all ntel guys 
#denom<- colSums(lam)
denom<-rep(1,nG)
probs<- t(t(lam)/denom)  # each column is the probs for a guy at column [j]

marg<-rep(NA,length=nrow(ntel))
for(j in 1:nrow(ntel)){
 cond<- exp( colSums(ntel[j,]*log(probs) -probs)) 
 marg[j]<- sum(cond*(1/nG) )
}
tel.loglik<- -1*sum(log(marg))

out<- out  + tel.loglik
}





out
}







