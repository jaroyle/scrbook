
tr <- seq(15, 85, length=10)
X <- cbind(rep(tr, each=length(tr)),
           rep(tr, times=length(tr)))    # trap coords
set.seed(10)
xlim <- c(0, 100); ylim <- c(0, 100)     # S is [0,100]x[0,100] square
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/1e4 # area of S
mu <- 50                                 # density (animals/unit area)
(N <- rpois(1, mu*A))                    # Generate N=75 as Poisson deviate
s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
#plot(X, xlim=xlim, ylim=ylim, pch="+")
#points(s, col=gray(0.5), pch=16)

sigma <- 5
lam0 <- 0.4
J <- nrow(X)
K <- 5
y <- array(NA, c(N, J, K))
for(j in 1:J) {
    dist <- sqrt((X[j,1]-s[,1])^2 + (X[j,2] - s[,2])^2)
    lambda <- lam0*exp(-dist^2/(2*sigma^2))
    for(k in 1:K) {
        y[,j,k] <- rpois(N, lambda)
    }
}
table(y)

table(apply(apply(y, c(1,2), sum)>0, 1, sum))


n <- apply(y, c(2,3), sum)
dimnames(n) <- list(paste("trap", 1:J, sep=""),
                    paste("night", 1:K, sep=""))
n[1:4,]

table(n)

plot(X, cex=rowSums(n), asp=1, xlim=xlim)
points(X, pch="+", cex=.5)
points(s, col=gray(0.5), pch=16)





# Analyze using scrUN()

library(scrbook)
library(coda)

# source("../../Rpackage/scrbook/R/scrUN.R")

set.seed(4569)
system.time({
fm1 <- scrUN(n=n, X=X, M=300, niter=250000, xlims=xlim, ylims=ylim,
             inits=list(lam0=0.4, sigma=5),
             updateY=TRUE,
#             priors=list(sigma=list("dgamma",
#                                    list(shape=0.001, rate=0.001)),
#                         lam0=list("dgamma",
#                                   list(shape=0.001, rate=0.001)),
#                         psi=list("dbeta",
#                                  list(shape1=1, shape2=1))),
             tune=c(0.27, 0.08, 10))
}) # 24000 it/hr

mc1 <- mcmc(fm1$sims)

plot(mc1)
plot(window(mc1, start=1001))

summary(mc1)
summary(window(mc1, start=1001))

rejectionRate(mc1)
rejectionRate(window(mc1, start=1001))



# debugonce(scrUN)





save(mc1, file="scrUNmc1.gzip")

#save(mc1, file="scrUNmc1prior.gzip")



ss1 <- summary(mc1)
out1 <- cbind(ss1$stat[,1:2], ss1$quant[,c(1,3,5)])

write.table(format(out1, digits=2), quote=FALSE, sep=" & ", eol="\\\\\n")










ls()
# load("scrUNmc1.gzip")
ls()


# No y updates
set.seed(4569)
system.time({
fm2 <- scrUN(n=n, X=X, M=300, niter=250000, xlims=xlim, ylims=ylim,
             inits=list(lam0=0.4, sigma=5),
             updateY=FALSE,
             priors=list(sigma=list("dgamma",
                                    list(shape=0.001, rate=0.001)),
                         lam0=list("dgamma",
                                   list(shape=0.001, rate=0.001)),
                         psi=list("dbeta",
                                  list(shape1=0.001, shape2=1))),
             tune=c(0.3, 0.6, 5))
}) # 40463 it/hr


mc2 <- mcmc(fm2$sims)
plot(mc2)
plot(window(mc2, start=1001))

summary(mc2)
summary(window(mc2, start=1001))

rejectionRate(mc2)
rejectionRate(window(mc2, start=1001))


save(mc2, file="scrUNmc2.gzip")

mc2.prior <- mc2
# save(mc2.prior, file="scrUNmc2prior.gzip")



ls()
# load("scrUNmc2.gzip")
ls()





ss2 <- summary(window(mc2, start=5001, end=25000))
out2 <- cbind(ss2$stat[,1:2], ss2$quant[,c(1,3,5)])

write.table(format(out2, digits=2), quote=FALSE, sep=" & ", eol="\\\\\n")


ss2.prior <- summary(window(mc2.prior, start=5001, end=25000))
out2.prior <- cbind(ss2.prior$stat[,1:2], ss2.prior$quant[,c(1,3,5)])

write.table(format(out2.prior, digits=2),
            quote=FALSE, sep=" & ", eol="\\\\\n")


plot(mc2.prior)














# Yes y, but no priors

set.seed(4569)
system.time({
fm3 <- scrUN(n=n, X=X, M=300, niter=6000, xlims=xlim, ylims=ylim,
             inits=list(lam0=0.3, sigma=0.01),
             updateY=TRUE,
             tune=c(0.003, 0.06, 0.15))
}) # 39700 it/hr

mc3 <- mcmc(fm3$sims)
plot(mc3)
plot(window(mc3, start=4001))
summary(mc3)
summary(window(mc3, start=4001))

rejectionRate(mc3)
rejectionRate(window(mc3, start=4001))



save(mc3, file="scrUNmc3.gzip")





# No y, no priors


# No y updates
set.seed(4569)
system.time({
fm4 <- scrUN(n=n, X=X, M=300, niter=6000, xlims=xlim, ylims=ylim,
             inits=list(lam0=0.3, sigma=0.01),
             updateY=FALSE,
             tune=c(0.003, 0.06, 0.15))
}) # 40463 it/hr


mc4 <- mcmc(fm4$sims)
plot(mc4)
plot(window(mc4, start=1001))
summary(mc4)
summary(window(mc4, start=1001))

rejectionRate(mc4)
rejectionRate(window(mc4, start=1001))


save(mc4, file="scrUNmc4.gzip")












# Analyze in JAGS


library(rjags)
dat1 <- list(n=n, X=X, J=J, K=K, M=150, xlim=xlim, ylim=ylim)
init1 <- function() {
    yi <- array(0, c(dat1$M, dat1$J, dat1$K))
    for(j in 1:dat1$J) {
        for(k in 1:dat1$K) {
            yi[sample(1:dat1$M, dat1$n[j,k]),j,k] <- 1
        }
    }
    list(sigma=runif(1, 1, 2), lam0=runif(1),
         y=yi, z=rep(1, dat1$M))
}
pars1 <- c("lam0", "sigma", "N", "mu")

system.time({
jm <- jags.model("SCmod1.jag", data=dat1, inits=init1, n.chain=1,
                 n.adapt=1000)
jc1.1 <- coda.samples(jm, pars1, n.iter=16000)
}) # 1032it/hr


plot(jc1.1)
summary(jc1.1)


save(jc1.1, file="scrUNjc1.1.gzip")



jc1.2 <- coda.samples(jm, pars1, n.iter=5000)








#










library(rjags)
dat2 <- list(n=n, X=X, J=J, K=K, M=300, xlim=xlim, ylim=ylim)
init2 <- function() {
    list(sigma=runif(1, 1, 2), lam0=runif(1),
         z=rep(1, dat2$M))
}
pars2 <- c("lam0", "sigma", "N", "mu")

system.time({
jm2 <- jags.model("SCmod2.jag", data=dat2, inits=init2, n.chain=1,
                 n.adapt=1000)
jc2.1 <- coda.samples(jm2, pars2, n.iter=10000)
}) # 14000 it/hr

plot(jc2.1)
summary(jc2.1)


save(jc2.1, file="scrUNjc2.1.gzip")


ls()
# load("scrUNjc2.1.gzip")
ls()





jc2.2 <- coda.samples(jm2, pars2, n.iter=5000)
plot(jc2.2)

summary(jc2.2)


jc2.3 <- coda.samples(jm2, pars2, n.iter=5000)

plot(jc2.3)
summary(jc2.3)


jc2.4 <- coda.samples(jm2, pars2, n.iter=5000)

plot(jc2.4)
summary(jc2.4)


jc2.5 <- coda.samples(jm2, pars2, n.iter=5000)

plot(jc2.5)
summary(jc2.5)




save.image("sim.RData")
















# Compare results

ls()
load("scrUNmc1.gzip")
load("scrUNmc2.gzip")
load("scrUNjc1.1.gzip")
load("scrUNjc2.1.gzip")
ls()


plot(mc1[,c("sigma", "lam0", "N")])
plot(jc1.1[,c("sigma", "lam0", "N")])


plot(mc2[,c("sigma", "lam0", "N")])
plot(jc2.1[,c("sigma", "lam0", "N")])


# Update y
summary(mc1[,"N"])$quant
summary(jc1.1[,"N"])$quant

# Don't update y
summary(mc2[,"N"])$quant
summary(jc2.1[,"N"])$quant




# Update y
summary(mc1)$quant   # R version
summary(jc1.1)$quant # JAGS version

# Don't update y
summary(mc2)$quant   # R version
summary(jc2.1)$quant # JAGS version





par(mfrow=c(2,2), mai=c(0.3, 0.7, 0.7, 0.2))
hist(as.vector(mc1[,"N"]), xlim=c(0, 200),
     main="R: y updated"); abline(v=50, col=4, lwd=2)
hist(as.vector(mc2[,"N"]), xlim=c(0, 200),
     main="R: y not updated"); abline(v=50, col=4, lwd=2)
hist(as.matrix(jc1.1)[,"N"], xlim=c(0, 200),
     main="JAGS: y updated"); abline(v=50, col=4, lwd=2)
hist(as.matrix(jc2.1)[,"N"], xlim=c(0, 200),
     main="JAGS: y not updated"); abline(v=50, col=4, lwd=2)







