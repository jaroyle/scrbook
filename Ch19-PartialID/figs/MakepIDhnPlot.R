### make figure partialID.fig.pIDhn
N <- 200  # Guys in S

# Dimensions of S 
xlimS <- c(0, 4)
ylimS <- c(0, 4)
cen <- c(mean(xlimS), mean(ylimS))

# Activity centers in S
set.seed(43)
s <- cbind(runif(N, xlimS[1], xlimS[2]),
           runif(N, ylimS[1], ylimS[2]))

dc <- apply(s, 1, function(x) sqrt((x[1]-cen[1])^2 + (x[2]-cen[2])^2))
tau <- 0.85
w0 <- 1
PrMark <- w0*exp(-dc^2/(2*tau^2))
w <- rbinom(N, 1, PrMark)

co <- seq(0.8, 3.2, len=6)
X <- cbind(rep(co, each=length(co)), rep(co, times=length(co)))

plot(xlimS, ylimS, xlab="Easting", ylab="Northing", cex=1.3, cex.lab=1.5, cex.axis=1.5, pch="")
points(s[w==0,], col="grey", pch=16, cex=1.3)
points(s[w==1,], col="black", pch=16, cex=0.8)
points(X, pch="+", cex=1.6)
points(x=2,y=2, pch="*", cex=2)

