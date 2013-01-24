


set.seed(666)
(n <- rbinom(3, size=10, prob=0.5))


dbinom(5, 10, 0.5)
N <- 10
n <- 5
p <- 0.5
factorial(N)/(factorial(n)*factorial(N-n))*p^n*(1-p)^(N-n)
exp(lgamma(N+1) - (lgamma(n+1) + lgamma(N-n+1)))*p^n*(1-p)^(N-n)
choose(N, n)*p^n*(1-p)^(N-n)




plot(0:15, dbinom(0:15, 10, 0.5), type="h", lwd=10, lend="butt",
     col=gray(0.5), xlab="n", ylab="Pr(x=n|N,p)")


set.seed(334)


plot(0:15, dbinom(0:15, 10, 0.5), type="h", lwd=10, lend="butt",
     col=gray(0.5), xlab="n", ylab="Pr(x=n|N,p)")
hist(rbinom(10000, 10, 0.5), freq=FALSE, add=TRUE)



N <- 20
p <- 0.5
x <- rbinom(100000, 20, 0.5)
p.hat <- x/20
hist(p.hat) # sampling distribution of p.hat
mean(p.hat) # expected value of the sampling distribution



y <- 3
(p^y * (1-p)^(N-y))*120
dbinom(y, N, p)


plot(0:N, dbinom(0:N, N, p), type="h", lwd=3, #ylim=c(0, 0.5),
     xlab="Number of captures", ylab="Probability")
abline(h=0, col=gray(0.5))




integrate(dnorm, -Inf, Inf, mean=0, sd=1)$value
sum(dbinom(0:5, size=5, p=0.1))



dbinom(7, 20, 0.35)

plot(0:20, dbinom(0:20, 20, 0.35), type="h", ylab="Probability",
     xlab="Number of shad caught (X)")

plot(table(rbinom(1000, 20, 0.35))/1000)

png("../figs/bin.png", width=7, height=7, units="in", res=400)
plot(0:20, dbinom(0:20, 20, 0.35), type="h", ylab="Probability",
     xlab=expression(paste("Number of shad caught (", italic(x), ") after 20 casts", sep="")), lwd=3, lend="butt",
     cex.lab=1.3)
abline(h=0, col=gray(0.8))
dev.off()

system("open ../figs/bin.png")





set.seed(30394)
nFish <- 50
mean <- 60
SD <- sqrt(5)
nSamples <- 30
sample.data <- matrix(NA, nFish, nSamples)
sample.means <- sample.SDs <- rep(NA, nSamples)
curve(dnorm(x, mean, SD), 45, 75, ylim=c(0, 0.4), xlim=c(50, 70),
     xlab="Fish length", ylab="Probability density", lwd=2)
for(i in 1:nSamples) {
  sample.data[,i] <- rnorm(nFish, mean, SD)
  sample.means[i] <- mean(sample.data[,i])
  sample.SDs[i] <- sd(sample.data[,i])
  curve(dnorm(x, sample.means[i], sample.SDs[i]), 45, 75, add=TRUE,
        col=i, lty=i)
}
curve(dnorm(x, mean, SD), 45, 75, add=TRUE, lwd=2)
#lines(density(sample.means))
dev.off()





# Monte Carlo average
set.seed(0343)
mean(rbinom(10000, 9, 0.35))




# Expected value
integrate(function(x) x*dnorm(x, 3, 1), -Inf, Inf)


# Expected value
sum(dbinom(0:100, 20, 0.35)*0:100)


# Variance
set.seed(340)
20*0.35*(1-0.35)             # Population variance
x <- rbinom(100000, 20, 0.35)
mean((x-mean(x))^2)          # Monte Carlo approximation






# Multinomial

set.seed(2321)
caphist.probs <- c("11"=0.09, "10"=0.21, "01"=0.21, "00"=0.49)
drop(rmultinom(1, 10, caphist.probs))





# Uniform on plane


set.seed(3656)
plot(runif(100), runif(100))




# Bivariate normal

library(mvtnorm)
set.seed(3)
mu <- c(0,0)
Sigma <- matrix(c(1, .9, .9, 1), 2, 2)
X1 <- cbind(rnorm(50, mu[1], Sigma[1,1]), # No correlation (rho=0)
            rnorm(50, mu[2], Sigma[2,2]))
X2 <- rmvnorm(50, mu, Sigma)              # rho=0.9



png("../figs/bvn.png", width=12, height=6, units="in", res=400)
par(mfrow=c(1,2), mai=c(0.1, 0.1, 0.1, 0.1))
plot(X1, axes=FALSE, frame=TRUE, xlim=c(-3, 3), ylim=c(-3, 3))
points(0, 0, pch=16, col="gray", cex=2)
plot(X2, axes=FALSE, frame=TRUE, xlim=c(-3, 3), ylim=c(-3, 3))
points(0, 0, pch=16, col="gray", cex=2)
dev.off()
system("open ../figs/bvn.png")


# Maximum likelihood


nSites1 <- 10
nSites2 <- 100
replicates <- 5000 # Substitute for infinity

lambda <- 4 # warblers/site
lambda.hat1 <- rep(NA, replicates)
lambda.hat2 <- rep(NA, replicates)

for(i in 1:replicates) {
    counts1 <- rpois(nSites1, lambda)
    counts2 <- rpois(nSites2, lambda)
    lambda.hat1[i] <- mean(counts1)
    lambda.hat2[i] <- mean(counts2)
}

plot(density(lambda.hat1), lty=1, ylim=c(0, 2))
lines(density(lambda.hat2), lty=2)






set.seed(3440)
lambda <- 3
y1 <- rpois(100, lambda)
negLogLike1 <- function(par) -sum(dpois(y1, par, log=TRUE))
starting.value <- c('lambda'=1)
optim(starting.value, negLogLike1)$par


set.seed(540)
nsites <- 100
elevation <- rnorm(100)
veght <- rnorm(100)
beta0 <- 1
beta1 <- -1
beta2 <- 0
lambda <- exp(beta0 + beta1*elevation + beta2*veght)
y2 <- rpois(nsites, lambda)
negLogLike2 <- function(pars) {
    beta0 <- pars[1]
    beta1 <- pars[2]
    beta2 <- pars[3]
    lambda <- exp(beta0 + beta1*elevation + beta2*veght)
    -sum(dpois(y2, lambda, log=TRUE))
}
starting.values <- c('beta0'=0, 'beta1'=0, 'beta2'=0)
optim(starting.values, negLogLike2)$par





# Joint, marginal, conditional distributions

X <- 0:20 # All possible values of X
Y <- 0:10  # All possible values of Y
lambda <- 0.6
p <- plogis(-0.62 + -2*Y) # p as function of Y
round(p,2)
joint <- matrix(NA, length(X), length(Y))
rownames(joint) <- paste("X=", X, sep="")
colnames(joint) <- paste("Y=", Y, sep="")

# Joint distribution [X,Y]
for(i in 1:length(Y)) {
    joint[,i] <- dbinom(X, 20, p[i]) * dpois(Y[i], lambda)
}
round(joint,2)
sum(joint)  # As dictated by law of total probability

# Marginal distributions [X] and [Y]
margX <- rowSums(joint)
round(margX, 2)

margY <- colSums(joint)
round(margY, 2)
all(margY==dpois(Y,lambda)) # Y is independent of X, but not vice versa

# Conditional distributions [X|Y] and [Y|X]
XgivenY <- joint/matrix(margY, nrow(joint), ncol(joint), byrow=TRUE)
round(XgivenY, 2)
YgivenX <- joint/matrix(margX, nrow(joint), ncol(joint))
round(YgivenX, 2)

colSums(XgivenY)
rowSums(YgivenX)









# Basic simulation


set.seed(36372)
Area <- 1                               # area of state-space (unit square)
x <- cbind(rep(seq(.1,.9,.2), each=5),  # trap locations
           rep(seq(.1,.9,.2), times=5))
p0 <- 0.3                               # baseline capture probability
sigma <- 0.05                           # Gaussian scale parameter
mu <- 50                                # population density
N <- rpois(1, mu*Area)                  # population size
s <- cbind(runif(N, 0, 1),              # activity centers in unit square
           runif(N, 0, 1))
K <- 5
y <- matrix(NA, N, nrow(x))             # capture data
for(i in 1:N) {
  d.ij <- sqrt((x[,1] - s[i,1])^2 +     # distance between x and s[i]
               (x[,2] - s[i,2])^2)
  p.ij <- p0*exp(-d.ij^2 / (2*sigma^2)) # capture probability
  y[i,] <- rbinom(nrow(x), K, p.ij)     # capture history for animal i
}


png("../figs/SCR0.png", width=7, height=7, units="in", res=400)
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(x, xlim=c(0,1), ylim=c(0,1), pch="+", cex=1.5, axes=FALSE, frame=TRUE)
points(s)
for(i in 1:N) {
    t1 <- y[i,]==1
    nj <- sum(t1)
    if(nj == 0)
       next
    points(s[i,,drop=FALSE], pch=16, col="gray")
    segments(rep(s[i,1],nj), rep(s[i,2],nj), x[t1,1], x[t1,2])
}
dev.off()
system("open ../figs/SCR0.png")



set.seed(36372)       # so that results can be reproduced
N <- 10               # population size
                      # create trap coordinates:
x <- cbind(rep(seq(0.1,0.9,0.2), each=5), rep(seq(0.1,0.9,0.2), times=5))
                      # generate individual home range centroids
s <- cbind(runif(N), runif(N))
                      # create nice graphic:
plot(x, pch= "+", xlim=c(0,1), ylim=c(0,1), xlab="Easting", ylab="Northing")
points(s, pch=16, col="blue")
for(t in 1:5) {
  points(cbind(rnorm(N, s[,1], 0.05), rnorm(N, s[,2], 0.05)), col="green",pch=20)
}


















# DAGS



# The full Monte, minus hyper
png("../figs/DAGr0.png", width=6, height=6.5, units="in", res=400)
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(0, 0, xlim=c(-0.5, 1.4), ylim=c(-.85, 0.9),
     asp=1, type="n", axes=FALSE,
     ann=FALSE, frame=TRUE)
symbols(0, 0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.8, expression(bold(s)))
text(0.5, 0.8, "Activity center location", pos=4)
arrows(0, 0.7, 0, 0.5, length=0.1, lwd=2)
symbols(0, 0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.4, expression(bold(u)))
text(0.5, 0.4, "Location of animal", pos=4)
arrows(0, 0.3, 0, 0.1, length=0.1)
symbols(0, 0, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0, "d")
text(0.5, 0, "Distance between trap and animal", pos=4)
arrows(0, -0.1, 0, -0.3, length=0.1)
symbols(0, -0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.4, "p")
text(0.5, -0.4, "Capture/detection probability", pos=4)
arrows(0, -0.5, 0, -0.7, length=0.1, lwd=2)
symbols(0, -0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.8, "y")
text(0.5, -0.8, "Capture/detection data", pos=4)
dev.off()
system("open ../figs/DAGr0.png")






# The full Monte
png("../figs/DAGr1.png", width=6, height=6.5, units="in", res=400)
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(0, 0, xlim=c(-0.5, 1.4), ylim=c(-.85, 0.9),
     asp=1, type="n", axes=FALSE,
     ann=FALSE, frame=TRUE)
symbols(0, 0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.8, expression(bold(s)))
text(0.5, 0.8, "Activity center location", pos=4)
arrows(0, 0.7, 0, 0.5, length=0.1, lwd=2)
symbols(-.4, 1, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, 1, expression(mu))
arrows(-.36, 0.965, -.09, 0.84, length=0.1, lwd=2)
symbols(.4, 1, rectangles=matrix(c(0.1, 0.1), 1), inches=FALSE, add=TRUE)
text(.4, 1, expression(italic(S)))
arrows(.35, 0.965, .09, 0.84, length=0.1, lwd=2)
symbols(0, 0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.4, expression(bold(u)))
text(0.5, 0.4, "Animal location", pos=4)
arrows(0, 0.3, 0, 0.1, length=0.1)
symbols(-.4, 0.6, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, 0.6, expression(tau))
arrows(-.36, 0.565, -.09, 0.45, length=0.1, lwd=2)
symbols(0, 0, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0, "d")
text(0.5, 0, "Distance between trap and animal", pos=4)
arrows(0, -0.1, 0, -0.3, length=0.1)
symbols(-.4, 0.2, rectangles=matrix(c(.1, .1), 1), inches=FALSE, add=TRUE)
text(-.4, 0.2, expression(bold(X)))
arrows(-.35, 0.18, -.09, 0.05, length=0.1)
symbols(0, -0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.4, "p")
text(0.5, -0.4, "Capture/detection probability", pos=4)
arrows(0, -0.5, 0, -0.7, length=0.1, lwd=2)
symbols(-.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, -0.2, expression(sigma))
arrows(-.36, -0.226, -.09, -0.35, length=0.1)
symbols(.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(.4, -0.2, expression(italic(p)[0]))
arrows(.36, -0.226, .09, -0.35, length=0.1)
symbols(0, -0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.8, "y")
text(0.5, -0.8, "Capture/detection data", pos=4)
dev.off()
system("open ../figs/DAGr1.png")





# Typical SCR
png("../figs/DAGr2.png", width=6, height=6.5, units="in", res=400)
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(0, 0, xlim=c(-0.5, 1.4), ylim=c(-.85, 0.9),
     asp=1, type="n", axes=FALSE,
     ann=FALSE, frame=TRUE)
symbols(0, 0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.8, expression(bold(s)))
text(0.5, 0.8, "Activity center location", pos=4)
arrows(0, 0.7, 0, 0.1, length=0.1)
symbols(-.4, 1, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, 1, expression(mu))
arrows(-.36, 0.965, -.09, 0.84, length=0.1)
symbols(.4, 1, rectangles=matrix(c(0.1, 0.1), 1), inches=FALSE, add=TRUE)
text(.4, 1, expression(italic(S)))
arrows(.35, 0.965, .09, 0.84, length=0.1)
#symbols(0, 0.4, circles=0.1, inches=FALSE, add=TRUE)
#text(0, 0.4, expression(bold(u)))
#text(0.5, 0.4, "Animal location", pos=4)
#arrows(0, 0.3, 0, 0.1, length=0.1)
#symbols(-.4, 0.6, circles=0.05, inches=FALSE, add=TRUE)
#text(-.4, 0.6, expression(tau))
#arrows(-.36, 0.565, -.09, 0.45, length=0.1)
symbols(0, 0, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0, "d")
text(0.5, 0, "Distance between trap and animal", pos=4)
arrows(0, -0.1, 0, -0.3, length=0.1)
symbols(-.4, 0.2, rectangles=matrix(c(.1, .1), 1), inches=FALSE, add=TRUE)
text(-.4, 0.2, expression(bold(X)))
arrows(-.35, 0.18, -.09, 0.05, length=0.1)
symbols(0, -0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.4, "p")
text(0.5, -0.4, "Capture/detection probability", pos=4)
arrows(0, -0.5, 0, -0.7, length=0.1)
symbols(-.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, -0.2, expression(sigma))
arrows(-.36, -0.226, -.09, -0.35, length=0.1)
symbols(.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(.4, -0.2, expression(italic(p)[0]))
arrows(.36, -0.226, .09, -0.35, length=0.1)
symbols(0, -0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.8, "y")
text(0.5, -0.8, "Capture/detection data", pos=4)
dev.off()
system("open ../figs/DAGr2.png")








# Typical distance sampling
png("../figs/DAGr3.png", width=6, height=6.5, units="in", res=400)
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(0, 0, xlim=c(-0.5, 1.4), ylim=c(-.85, 0.9),
     asp=1, type="n", axes=FALSE,
     ann=FALSE, frame=TRUE)
#symbols(0, 0.8, circles=0.1, inches=FALSE, add=TRUE)
#text(0, 0.8, expression(bold(s)))
#text(0.5, 0.8, "Activity center location", pos=4)
#arrows(0, 0.7, 0, 0.5, length=0.1)
#symbols(-.4, 1, circles=0.05, inches=FALSE, add=TRUE)
#text(-.4, 1, expression(mu))
#arrows(-.36, 0.965, -.09, 0.84, length=0.1)
#symbols(.4, 1, rectangles=matrix(c(0.1, 0.1), 1), inches=FALSE, add=TRUE)
#text(.4, 1, expression(italic(S)))
#arrows(.35, 0.965, .09, 0.84, length=0.1)
symbols(0, 0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0.4, expression(bold(u)))
text(0.5, 0.4, "Animal location", pos=4)
arrows(0, 0.3, 0, 0.1, length=0.1)
#symbols(-.4, 0.6, circles=0.05, inches=FALSE, add=TRUE)
#text(-.4, 0.6, expression(tau))
#arrows(-.36, 0.565, -.09, 0.45, length=0.1)
symbols(0, 0, circles=0.1, inches=FALSE, add=TRUE)
text(0, 0, "d")
text(0.5, 0, "Distance between trap and animal", pos=4)
arrows(0, -0.1, 0, -0.3, length=0.1)
symbols(-.4, 0.2, rectangles=matrix(c(.1, .1), 1), inches=FALSE, add=TRUE)
text(-.4, 0.2, expression(bold(X)))
arrows(-.35, 0.18, -.09, 0.05, length=0.1)
symbols(0, -0.4, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.4, "p")
text(0.5, -0.4, "Capture/detection probability", pos=4)
arrows(0, -0.5, 0, -0.7, length=0.1)
symbols(-.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(-.4, -0.2, expression(sigma))
arrows(-.36, -0.226, -.09, -0.35, length=0.1)
symbols(.4, -0.2, circles=0.05, inches=FALSE, add=TRUE)
text(.4, -0.2, expression(italic(p)[0]))
arrows(.36, -0.226, .09, -0.35, length=0.1)
symbols(0, -0.8, circles=0.1, inches=FALSE, add=TRUE)
text(0, -0.8, "y")
text(0.5, -0.8, "Capture/detection data", pos=4)
dev.off()
system("open ../figs/DAGr3.png")









