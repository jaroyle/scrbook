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


