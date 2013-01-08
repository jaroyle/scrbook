library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
# this takes 341 seconds on a standard CPU circa 2011
unix.time(out<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,keepz=TRUE))

Sx<-out$sims.list$s[,,1]
Sy<-out$sims.list$s[,,2]
w<- out$sims.list$z
obj<-list(Sx=Sx,Sy=Sy,z=w)
SCRdensity(obj)

allx<-Sx[1:length(Sx)]
ally<-Sy[1:length(Sy)]

plot(allx,ally,pch=" ")

points(Sx[,1],Sy[,1],pch=".")



 SCRdensity<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, 
    Yu = NULL, scalein = 100, scaleout = 100, col="gray",ncolors = 10,whichguy=NULL) 
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
guy<-col(Sxout)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
if(is.null(whichguy)){
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
}
else{
    Dn<-table(Sxout[guy==whichguy],Syout[guy==whichguy] )/niter
}

    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
 if (col == "gray") {
        cc <- seq(3, 17, , 10)/20
        cc <- rev(gray(cc))
    }
    else cc <- terrain.colors(ncolors)

    image(xg, yg, Dn, col = cc)
    image.scale(Dn, col = cc)
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}







png("wolv_post_s.png",width=5.25,height=5.25, units="in", res=400)

par(mfrow=c(1,1),mar=c(2,2,2,6))
X<-wolverine$wtraps[,2:3]
X[,1]<-(X[,1]-min(X[,1]))/10000
X[,2]<-(X[,2]-min(X[,2]))/10000

SCRdensity(obj,whichguy=1)
points(X,pch=20,cex=1,col="white")
points(X[30,],pch=1,cex=2,col="black")
dev.off()







par(mfrow=c(2,1),mar=c(2,2,2,6))
X<-wolverine$wtraps[,2:3]
X[,1]<-(X[,1]-min(X[,1]))/10000
X[,2]<-(X[,2]-min(X[,2]))/10000

SCRdensity(obj,whichguy=22)
points(X,pch=20,cex=1,col="white")
#points(X[30,],pch=1,cex=2,col="black")







