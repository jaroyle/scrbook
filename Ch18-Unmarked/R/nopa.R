

library(scrbook)
library(rjags)

data(nopa)
str(nopa)

nj <- rowSums(nopa$n)

pdf("../figs/nopaCounts.pdf", width=6, height=4)
par(mai=c(0.9, 0.9, 0.1, 0.1))
plot(nopa$X, asp=1, pch="+", col=gray(0.5), cex=0.8,
     xlab="Easting (m)", ylab="Northing (m)", cex.lab=1.4)
text(nopa$X[nj>0,], label=nj[nj>0], cex=1.5)
dev.off()
system("open ../figs/nopaCounts.pdf")








# Fit model using custom code

set.seed(450)
fmy1 <- scrUN(n=nopa$n, X=nopa$X, M=250, updateY=TRUE, niters=250000,
              xlims=c(-600, 600), ylims=c(-400, 400),
              inits=list(sigma=rnorm(1, 100)),
              tune=c(9, 0.05, 300))

mcy1 <- mcmc(fmy1$sims)
plot(mcy1)
summary(mcy1)

rejectionRate(mcy1)
rejectionRate(window(mcy1, start=4001))






save(mcy1, file="mcy1.gzip")

ls()
# load("mcy1.gzip")
ls()




set.seed(450)
fmnoy1 <- scrUN(n=nopa$n, X=nopa$X, M=400, updateY=FALSE, niters=5000,
              xlims=c(-600, 600), ylims=c(-400, 400),
              inits=list(sigma=rnorm(1, 100)),
              tune=c(9, 0.05, 300))

mcnoy1 <- mcmc(fmnoy1$sims)
plot(mcnoy1)
summary(mcnoy1)

rejectionRate(mcnoy1)
rejectionRate(window(mcnoy1, start=1001))

summary(as.matrix(mcnoy1)[,"N"])



save(mcnoy1, file="mcnoy1.gzip")


ls()
# load("mcnoy1.gzip")
ls()


















# JAGS

library(rjags)

paste("y[", 1:200, ",j,k]", sep="", collapse=", ")


dat1 <- list(n = nopa$n, X = nopa$X, M=200, J=nrow(nopa$n), K=ncol(nopa$n),
             xlim=c(-600, 600), ylim=c(-400, 400))

init1 <- function() {
    n <- dat1$n
    J <- nrow(n)
    K <- ncol(n)
    M <- dat1$M
    y <- array(0L, c(M, J, K))
    for(j in 1:J) {
        for(k in 1:K) {
            y[sample(1:M, n[j,k]),j,k] <- 1
        }
    }
    list(y = y, sigma=rnorm(1, 100), lam0=0.5, z=rep(1, M))
}


str(dat1)
str(init1())

pars1 <- c("sigma", "lam0", "N", "ED")

library(parallel)

cl1 <- makeCluster(3) # Open 3 parallel R instances

clusterExport(cl1, c("dat1", "init1", "pars1"))

system.time({
out1 <- clusterEvalQ(cl1, {
    library(rjags)
    jm <- jags.model("nopa1.jag", dat1, init1, n.chains=1, n.adapt=50)
    jc <- coda.samples(jm, pars1, n.iter=50)
    return(as.mcmc(jc))
})
})

mc1 <- mcmc.list(out1)
plot(mc1)


jm1 <- jags.model("nopa1.jag", dat1, init1, n.chains=1, n.adapt=50)
jc1 <- coda.samples(jm1, pars1, n.iter=50)


plot(jc)
summary(jc)


stopCluster(cl1)









# Fit the model without the latent encounter histories


dat2 <- list(n = nopa$n, X = nopa$X, M=300, J=nrow(nopa$n), K=ncol(nopa$n),
             xlim=c(-600, 600), ylim=c(-400, 400))


init2 <- function() {
    list(sigma=rnorm(1, 100), lam0=0.5, z=rep(1, dat2$M))
}

cl2 <- makeCluster(3) # Open 3 parallel R instances

clusterExport(cl2, c("dat2", "init2", "pars1"))

system.time({
out2 <- clusterEvalQ(cl2, {
    library(rjags)
    jm <- jags.model("nopa2.jag", dat2, init2, n.chains=1, n.adapt=500)
    jc <- coda.samples(jm, pars1, n.iter=55500)
    return(as.mcmc(jc))
})
})


mc2 <- mcmc.list(out2)

plot(mc2)
summary(mc2)

sort(table(as.matrix(mc2)[,"N"]))



system.time({
out2.2 <- clusterEvalQ(cl2, {
    jc <- coda.samples(jm, pars1, n.iter=7000)
    return(as.mcmc(jc))
})
})


mc2.2 <- mcmc.list(out2.2)

plot(mc2.2, ask=TRUE)
summary(mc2.2)







stopCluster(cl2)












dat2 <- list(n = nopa$n, X = nopa$X, M=300, J=nrow(nopa$n), K=ncol(nopa$n),
             xlim=c(-600, 600), ylim=c(-400, 400))


init2 <- function() {
    list(sigma=rnorm(1, 100), lam0=runif(1), z=rep(1, dat2$M))
}

cl2 <- makeCluster(3) # Open 3 parallel R instances

clusterExport(cl2, c("dat2", "init2", "pars1"))

system.time({
out2 <- clusterEvalQ(cl2, {
    library(rjags)
    jm <- jags.model("nopa2i.jag", dat2, init2, n.chains=1, n.adapt=500)
    jc <- coda.samples(jm, pars1, n.iter=2500)
    return(as.mcmc(jc))
})
}) # 6242


mc2i <- mcmc.list(out2)

plot(mc2i)
summary(mc2i)
HPDinterval(mc2i)

summary(mc2)

autocorr.plot(mc2i)
crosscorr.plot(mc2i)

system.time({
out2 <- clusterEvalQ(cl2, {
    jc <- coda.samples(jm, pars1, n.iter=7500)
    return(as.mcmc(jc))
})
})




round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 0.001, 0.001), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 0.01, 0.01), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 0.1, 0.1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 10, 10), 3)


round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 10), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 0.1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 0.01), 3)


round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 10, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 0.1, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 0.01, 0.01), 3)



round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), .01, 100), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), .1, 10), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 1, 1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 10, 0.1), 3)
round(qgamma(c(0.025, 0.25, 0.5, 0.75, 0.975), 100, 0.01), 3)





stopCluster(cl2)


















# New analysis with Pr(y=1|x)=HN and x~Norm(s, tau)

datD1 <- list(n = nopa$n, X = nopa$X, M=200,
              J=nrow(nopa$n), K=ncol(nopa$n),
             xlim=c(-600, 600), ylim=c(-400, 400))


initD1 <- function() {
    n <- datD1$n
    J <- nrow(n)
    K <- ncol(n)
    M <- datD1$M
    y <- array(0L, c(M, J, K))
    s <- cbind(rnorm(M), rnorm(M))
    for(j in 1:J) {
        for(k in 1:K) {
            y[sample(1:M, n[j,k]),j,k] <- 1
        }
    }
    list(y = y, sigma=rnorm(1, 500), tau=rnorm(1, 1),
         s=s, z=rep(1, M))
}


parsD1 <- c("sigma", "tau", "N", "ED")



clD1 <- makeCluster(3) # Open 3 parallel R instances

clusterExport(clD1, c("datD1", "initD1", "parsD1"))

system.time({
outD1 <- clusterEvalQ(clD1, {
    library(rjags)
    jm <- jags.model("nopaD1.jag", datD1, initD1, n.chains=1, n.adapt=500)
    jc <- coda.samples(jm, parsD1, n.iter=2500)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD1 <- mcmc.list(outD1)

plot(mcD1)
summary(mcD1)







system.time({
outD1.2 <- clusterEvalQ(clD1, {
    jc <- coda.samples(jm, parsD1, n.iter=7000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD1.2 <- mcmc.list(outD1.2)






system.time({
outD1.3 <- clusterEvalQ(clD1, {
    jc <- coda.samples(jm, parsD1, n.iter=15000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD1.3 <- mcmc.list(outD1.3)
plot(mcD1.3)
summary(mcD1.3)






system.time({
outD1.4 <- clusterEvalQ(clD1, {
    jc <- coda.samples(jm, parsD1, n.iter=15000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD1.3 <- mcmc.list(outD1.3)
plot(mcD1.3)
summary(mcD1.3)








stopCluster(clD1)




















library(parallel)



clD2 <- makeCluster(3) # Open 3 parallel R instances

clusterExport(clD2, c("datD1", "initD1", "parsD1"))

system.time({
outD2 <- clusterEvalQ(clD2, {
    library(rjags)
    jm <- jags.model("nopaD2.jag", datD1, initD1, n.chains=1, n.adapt=500)
    jc <- coda.samples(jm, parsD1, n.iter=2500)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD2 <- mcmc.list(outD2)

plot(mcD2)
summary(mcD2)



system.time({
outD2.2 <- clusterEvalQ(clD2, {
    jc <- coda.samples(jm, parsD1, n.iter=5000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD2.2 <- mcmc.list(outD2.2)
plot(mcD2.2)



system.time({
outD2.3 <- clusterEvalQ(clD2, {
    jc <- coda.samples(jm, parsD1, n.iter=15000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD2.3 <- mcmc.list(outD2.3)
plot(mcD2.3)






system.time({
outD2.4 <- clusterEvalQ(clD2, {
    jc <- coda.samples(jm, parsD1, n.iter=100000)
    return(as.mcmc(jc))
})
}) # 1000it/hr


mcD2.4 <- mcmc.list(outD2.4)
plot(mcD2.4)









stopCluster(clD2)









edr <- function(sigma, r=150) {
    ea <- 2 * pi * integrate(function(x, sig=sigma)
                    exp(-x^2/(2*sig^2))*x, 0, r)$value
#    a <- ea/(pi*r^2)
    sqrt(ea/pi)
#    a
}

# NOPA EDR between 50 and 100, so sigma between 14 and 77

edr(14)
edr(77)

curve(dnorm(x, 75, 9), 50, 100, ylim=c(0, 0.04))

curve(dgamma(x, 50, 4), 50, 100)

curve(dgamma(x, 30, 3), 0, 100)



bigN7510 <- rnorm(10000, 75, 10)

nll <- function(pars) {
    a <- pars[1]
    b <- pars[2]
    -sum(dgamma(bigN7510, a, b, log=TRUE))
}

optim(c(5, 5), nll)


curve(dgamma(x, 55, 0.75), 50, 100)
