



npts <- 500
s1 <- c(0,-2)
s2 <- c(0,2)
sigma <- 0.8

# Movement outcomes
set.seed(23)
u1 <- cbind(rnorm(npts, s1[1], sigma), rnorm(npts, s1[2], sigma))
u2 <- cbind(rnorm(npts, s2[1], sigma), rnorm(npts, s2[2], sigma))

# Trap locations
xl <- 5
xn <- 10
xr <- seq(-xl, xl, length=xn)
x <- cbind(rep(xr, each=xn), rep(xr, times=xn))


# Simulate count data (T=1)
library(scrbook)
dist <- e2dist(rbind(s1,s2), x)
lam0 <- 10
lambda <- lam0*exp(-dist^2/(2*sigma^2))
n <- rpois(nrow(x), colSums(lambda))



# Figures

pdf("../figs/heuristic.pdf", width=6, height=6)
par(mai=rep(0.1,4))
plot(0, type="n", xlim=c(-xl-1, xl+1), asp=1, ann=FALSE,
     axes=FALSE, frame=TRUE)
points(rbind(s1, s2), pch=16, cex=1.2)
points(x, pch="+", cex=0.9, col=gray(0.5))
text(x, labels=ifelse(y>0, y, ""), cex=2)
dev.off()
system("open ../figs/heuristic.pdf")


png("../figs/heuristic2.png", width=6, height=3, res=400, units="in")
par(mfrow=c(1,2), mai=rep(0.1,4))
plot(0, type="n", xlim=c(-xl-1, xl+1), asp=1, ann=FALSE,
#     xlab="Easting", ylab="Northing",
     axes=FALSE, frame=TRUE)
points(x, pch="+", cex=0.3)
points(u1, pch=16, col=rgb(0,0,0,0.3))
points(u2, pch=16, col=rgb(0,0,0,0.3))
plot(0, type="n", xlim=c(-xl-1, xl+1), asp=1, ann=FALSE,
     axes=FALSE, frame=TRUE)
points(x, pch="+", cex=0.3)
text(x, labels=ifelse(y>0, y, ""))
dev.off()



