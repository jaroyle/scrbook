


(x <- rbinom(3, size=10, prob=0.5))







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







png("../figs/sampleDists.png", width=7, height=7, units="in", res=400)
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


