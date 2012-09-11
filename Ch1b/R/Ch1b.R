


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
     xlab="Number of shad caught (x)")

plot(table(rbinom(1000, 20, 0.35))/1000)

png("../figs/bin.png", width=7, height=7, units="in", res=400)
plot(0:20, dbinom(0:20, 20, 0.35), type="h", ylab="Probability",
     xlab="Number of shad caught (x) after 20 casts", lwd=3, lend="butt",
     cex.lab=1.3)
abline(h=0, col=gray(0.8))
dev.off()





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
set.seed(344)
integrate(function(x) x*dnorm(x, 3, 1), -Inf, Inf)


# Expected value
set.seed(344)
integrate(function(x) x*dnorm(x, 3, 1), -Inf, Inf)


