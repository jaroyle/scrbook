
# Make a cool figure

library(gdistance)

len <- 100
sp <- seq(0,1,len=len)
u <- cbind(rep(sp,len), rep(sp,each=len))


fu <- function(x, b=5) {
    exp(0.5 + b*x[1] + -b*x[1]^2)
}
elev <- 1/apply(u, 1, fu)
elev <- matrix(elev, len, len)
plot(raster(elev))
range(elev)

sigma <- 0.05

p <- apply(u, 1, function(x) {
    d <- sqrt((x[1]-0.5)^2 + (x[2]-0.5)^2)
    p <- 0.5*exp(-d^2/(2*sigma))
    p
})
p <- matrix(p, len)
image(sp, sp, p)


elev.z <- elev-min(elev)
elev.z <- elev.z/max(elev.z)
cost <- exp(20*raster(matrix(elev.z, len)))
tr <- transition(cost, function(x) mean(1/x), directions=8)
tr1 <- geoCorrection(tr)
plot(raster(tr1))

lcd <- costDistance(tr1, c(0.5,0.5), u)
lcdr <- raster(matrix(lcd, len, len))
#plot(lcdr)

p2 <- 0.5*exp(-lcd^2/(2*sigma))
p2 <- matrix(p2, len)
#image(sp, sp, p2)


png("../figs/distort.png", width=6, height=1.8, units="in", res=400)
op <- par(mfrow=c(1,3), mai=c(0.1,0.1,0.1,0.1))
pr <- seq(0, 1, by=0.1)
contour(sp, sp, t(p), levels=pr, axes=FALSE, frame=TRUE, xlim=c(0,1),asp=1)
plot(raster(elev.z), axes=FALSE)
contour(sp, sp, p2, levels=pr, axes=FALSE, frame=TRUE, xlim=c(0,1),asp=1)
par(op)
dev.off()

