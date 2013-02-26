

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
})


mc2 <- mcmc.list(out2)


