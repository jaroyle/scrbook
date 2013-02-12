SCRdensity <-
function(obj,nx=30,ny=30,Xl=NULL,Xu=NULL,Yl=NULL,Yu=NULL,scalein=100,scaleout=100,ncolors=10){

Sxout<-obj$Sx
Syout<-obj$Sy
z<-obj$z
niter<-nrow(z)

if(is.null(Xl)){
Xl<-min(Sxout)*.999
Xu<-max(Sxout)*1.001
Yl<-min(Syout)*.999
Yu<-max(Syout)*1.001
}
xg<-seq(Xl,Xu,,nx)
yg<-seq(Yl,Yu,,ny)

Sxout<-cut(Sxout[z==1],breaks=xg)
Syout<-cut(Syout[z==1],breaks=yg)
Dn<-table(Sxout,Syout)/niter  # Dn = avg # guy (posterior)
area<-  (yg[2]-yg[1])*(xg[2]-xg[1])*scalein # this is in sq km now

Dn<- (Dn/area)*scaleout   # now wolverines per 100 km
cat("mean: ",mean(Dn),fill=TRUE)
par(mar=c(3,3,3,6))
image(xg,yg,Dn,col=terrain.colors(ncolors))
image.scale(Dn,col=terrain.colors(ncolors))
box()
return(list(grid=cbind(xg,yg),Dn=Dn))
}
