



npts <- 500
s1 <- c(0,-2)
s2 <- c(0,2)
sigma <- 0.8

set.seed(23)
u1 <- cbind(rnorm(npts, s1[1], sigma), rnorm(npts, s1[2], sigma))
u2 <- cbind(rnorm(npts, s2[1], sigma), rnorm(npts, s2[2], sigma))

xl <- 5
xn <- 10
xr <- seq(-xl, xl, length=xn)
x <- cbind(rep(xr, each=xn), rep(xr, times=xn))





# Figures


png("../figs/heuristic.png", width=5, height=5, res=400, units="in")
par(mai=rep(0.1,4))
plot(0, type="n", xlim=c(-xl-1, xl+1), asp=1, ann=FALSE,
#     xlab="Easting", ylab="Northing",
     axes=FALSE, frame=TRUE)
points(x, pch="+")
points(u1, pch=16, col=rgb(0,0,0,0.3))
points(u2, pch=16, col=rgb(0,0,0,0.3))
dev.off()

