SCRgof <-
function(out,nx=6,ny=6,traplocs=NULL,buffer=2,Xl=NULL,Xu=NULL,Yl=NULL,Yu=NULL){
### need to have a "Discard" or shrink option, or provide trap locs?

## works with SCRf.fn output
## NEEDS "S" and state-space grid G --
## OR OR OR -- Just S if continuous space.
### needs "z" too
S<-out$s
Sxout<-S[,,1]
Syout<-S[,,2]
z<-out$w
Xl<- min(traplocs[,1]) - buffer
Xu<- max(traplocs[,1]) + buffer
Yl<- min(traplocs[,2]) - buffer
Yu<- max(traplocs[,2]) + buffer

niter<-nrow(z)

xg<-seq(Xl,Xu,,nx)
yg<-seq(Yl,Yu,,ny)
area<- (Xu-Xl)*(Yu-Yl)

Sxout2<-cut(Sxout[z==1],breaks=xg)
Syout2<-cut(Syout[z==1],breaks=yg)

Dn<-table(Sxout2,Syout2)/niter

image(xg,yg,Dn,col=terrain.colors(10))
image.scale(Dn,col=terrain.colors(10))

stat2<-statsim2<-stat<-statsim<-rep(NA,niter)


for(i in 1:niter){

N<- sum(z[i,])
D<- N/area
E<- N/( (nx-1)*(ny-1) )

inside  <- (Sxout[i,] < Xu & Sxout[i,] > Xl) & (Syout[i,] < Yu & Syout[i,] < Yl)

Dn<- table(cut(Sxout[i,][z[i,]==1 &inside],breaks=xg),cut(Syout[i,][z[i,]==1 &inside],breaks=yg))
Dnv<-Dn[1:length(Dn)]

E<-mean(Dnv)
stat[i]<-  (var(Dnv)/mean(Dnv))
stat2[i]<-   sum(   (sqrt(Dnv) - sqrt(E))^2 )

Sxsim<-runif(sum(z[i,][inside]),Xl,Xu)
Sysim<-runif(sum(z[i,][inside]),Yl,Yu)

Dnsim<- table(cut(Sxsim,breaks=xg),cut(Sysim,breaks=yg))
Dnsimv<-Dnsim[1:length(Dnsim)]
statsim[i]<- (var(Dnsimv)/mean(Dnsimv))
statsim2[i]<-  sum(   (sqrt(Dnsimv) - sqrt(E))^2)


}
out<-cbind(data=stat2,newdata=statsim2)
cat("Cluster index observed: ",mean(stat),fill=TRUE)
cat("Cluster index simulated: ",mean(statsim),fill=TRUE)
cat("P-value  index of dispersion: ",mean(statsim>stat),fill=TRUE)
cat("P-value2 freeman-tukey: ",mean(statsim2>stat2),fill=TRUE)
invisible(out)

}
