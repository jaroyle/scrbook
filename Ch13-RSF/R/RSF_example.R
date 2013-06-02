RSF_example <-
function(){
set.seed(1234)

gr<-make.statespace(minx=1,maxx=40,miny=1,maxy=40,nx=40,ny=40)
### Note: expand.grid puts the points in a different order.
### if you use expand.grid then you have to change the call to image() below
##gr<-as.matrix(expand.grid(1:40,1:40))

Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
z<-t(chol(V))%*%rnorm(1600)
spatial.plot(gr,z)

s<- matrix(c(
  9, 31,
  8, 18,
  9,  6,
 20, 24,
 19, 12,
 31, 32,
 31, 17,
 32,  7),ncol=2,byrow=TRUE)

# we combine distance and a covariate into a model for space usage
# individuals use each pixel in proportion to exp( distance + covariate)
# as in the following expression. We simulate counts here, as if we were
# recording telemetry location frequencies:
alpha0 <- -2
sigma<- 2
alpha2<- 1
Ntel<-8
nsim<-100
Nfixes<-20
N<- 100

#Sid<- sample(1:1600,N,replace=TRUE)
#S<-gr[Sid,]

## To simulate a new set of telemetered individuals do this:
#tel.guys<-sample(Sid,Ntel)
#sid<-tel.guys ####sample(1:1600,Ntel,replace=TRUE)
sid<- c( 330,  303,  355,  777,  749, 1209, 1224, 1274)
s<-gr[sid,]

n<-matrix(NA,nrow=Ntel,ncol=1600)

i<-1  # do this for each telemetered guy to simulate a number of fixes.
      # note that n = 0 for most of the landscape
par(mfrow=c(3,3))
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
   d<- Dmat[sid[i],]
   lam<- exp(1 - (1/(2*sigma*sigma))*d*d + alpha2* z) 
   n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
   par(mar=c(3,3,3,6))
   #spatial.plot(gr,lam)
   lammat[i,]<-lam
   img<- matrix(lam,nrow=40,ncol=40,byrow=FALSE)
   image(1:40,1:40,rot(img),col=terrain.colors(10))
}


#png("habitat.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1),mar=c(3,3,3,6))
image(1:40,1:40,rot(matrix(z,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
image.scale(z,col=terrain.colors(10))
points(s,pch=20)
#dev.off()


#png("homeranges8.png",width=7,height=7, units="in", res=400)
par(mfrow=c(1,1))
tot<- apply(lammat,1,sum)
lammat<-lammat/tot
lamtot<-apply(lammat,2,sum)
image(1:40,1:40,rot(matrix(lamtot,40,40,byrow=FALSE)),col=terrain.colors(10),xlab=" ",ylab=" ")
points(s,pch=20)
#dev.off()

}
