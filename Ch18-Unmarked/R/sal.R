
# Load data
load("salrem11.gzip")
ls()


str(dfus11.rem)
str(pt.data)


# Do names match?
all(pt.data$SiteName %in% rownames(dfus11.rem))


# Put them in same order, get rid of 1 site
sal.n <- dfus11.rem[pt.data$SiteName,,]
str(sal.n)

all(rownames(sal.n) == pt.data$SiteName) # check


apply(sal.n, c(1,3), sum, na.rm=TRUE)

apply(sal.n, c(3), sum, na.rm=TRUE)


# plots

plot(coords.x2 ~ coords.x1, pt.data, asp=1)




# Just use data from stream 27!!!!!!!!!
stream27 <- pt.data[16:79,]
plot(coords.x2 ~ coords.x1, stream27, asp=1)

sal.n27 <- dfus11.rem[stream27$SiteName,,]


pdf("../figs/saln27.pdf", width=6, height=4)
par(mfrow=c(1,2), mai=c(0.1, 0.1, 0.2, 0.1))
plot(coords.x2 ~ coords.x1, stream27, pch="+", cex=0.5, asp=1,
     type="n",
     axes=FALSE, frame=TRUE,
     xlab="", ylab="",
     main="Occasion 1")
text(coords.x2 ~ coords.x1, stream27, cex=0.8,
     label=rowSums(sal.n27[,,1], na.rm=TRUE))
plot(coords.x2 ~ coords.x1, stream27, pch="+", cex=0.5, asp=1,
     type="n",
     axes=FALSE, frame=TRUE,
     xlab="", ylab="",
     main="Occasion 2")
text(coords.x2 ~ coords.x1, stream27, cex=0.8,
     label=rowSums(sal.n27[,,2], na.rm=TRUE))
#plot(coords.x2 ~ coords.x1, stream27, pch="+", cex=0.5, asp=1,
#     type="n",
#     axes=FALSE, frame=TRUE,
#     xlab="", ylab="",
#     main="Occasion 3")
#text(coords.x2 ~ coords.x1, stream27,
#     label=rowSums(sal.n27[,,3], na.rm=TRUE))
dev.off()
system("open ../figs/saln27.pdf")





# Coordinates of traps
# In this case, every location in the state-space has a trap, so the number
# of places in the state-space G is equal to the number of trap, i.e. G=J
X <- data.matrix(stream27[,c("coords.x1", "coords.x2")])
G <- J <- nrow(X)

# Distance between all locations in the state-space
distmat <- matrix(NA, G, G)
for(g1 in 1:G) {
    for(g2 in 1:G) {
        distmat[g1,g2] <- sqrt((X[g1,1]-X[g2,1])^2 + (X[g1,2]-X[g2,2])^2)
    }
}

summary(as.vector(distmat))

seg <- 1:G # ID of each point in state-space


# JAGS
library(rjags)

dat1 <- list(n=sal.n27, distmat=distmat, seg=seg, PrSeg=rep(1/G, G),
             J=J, G=G, K=3, M=200)
str(dat1)

n.jk <- apply(dat1$n, c(1,3), sum, na.rm=TRUE)
xxx <- rowSums(n.jk)/sum(n.jk)

set.seed(53450)
si <- sample(seg, dat1$M, replace=TRUE, prob=xxx)
ui <- matrix(NA, dat1$M, 3)
yi <- array(0L, c(dat1$M, dat1$G, 3))



for(i in 1:dat1$M) {
#    cat("guy", i, "\n")
    for(j in 1:dat1$J) {
        for(k in 1:dat1$K) {
            if((sum(yi[1:i,j,k]) < n.jk[j,k]) & (sum(yi[i,1:j,k]) < 1))
#            if((sum(yi[1:i,j,k]) < dat1$M) & (sum(yi[i,1:j,k]) < 1))
                yi[i,j,k] <- 1
        }
    }
}


all(apply(yi, c(2,3), sum) == n.jk)

(apply(yi, c(2,3), sum))


cbind(apply(yi, c(2,3), sum), n.jk)


yi[1,,1]

for(i in 1:dat1$M) {
    for(k in 1:dat1$K) {
        if(any(yi[i,,k]>0))
            ui[i,k] <- which(yi[i,,k] == 1)
        else {
            PrU <- exp(-distmat[si[i],]^2/(2*1000))
            ui[i,k] <- sample(dat1$seg, 1, prob=PrU)
        }
    }
}




init1 <- function() list(p=runif(1), psi=runif(1),
                         tau=500, #runif(1, 10, 20),
                         phi=runif(1),
                         s=si,
                         u=ui,
#                         y=yi,
                         z=matrix(1, dat1$M, dat1$K))
str(init1())

pars1 <- c("phi", "tau", "p", "Ntot")

system.time({
    jm1 <- jags.model("sal1.jag", data=dat1, inits=init1, n.chains=1,
                      n.adapt=100)
    jc1 <- coda.samples(jm1, pars1, n.iter=100)
})

plot(jc1, ask=TRUE)


system.time({
    jc2 <- coda.samples(jm1, pars1, n.iter=1000)
})


summary(jc2)
plot(jc2, ask=TRUE)

summary(window(jc2, start=901), ask=TRUE)
plot(window(jc2, start=901), ask=TRUE)



system.time({
    jc3 <- coda.samples(jm1, pars1, n.iter=1000)
})

summary(jc3)
plot(jc3)




# unmarked

library(unmarked)


sal.n27mat <- matrix(sal.n27, 64)
umf <- unmarkedFrameGMM(y=sal.n27mat,
                        yearlySiteCovs=list(visit=
                        matrix(c('1','2','3'), 64, 3,
                               byrow=TRUE)),
                        numPrimary=3, type="removal")


fm1 <- gmultmix(~1, ~1, ~1, umf)
fm1


fm2 <- gmultmix(~1, ~visit, ~1, umf)
fm2

re2 <- ranef(fm2)
plot(re2)

re2mode <- bup(re2, "mode")







# visit 1

umf1 <- unmarkedFrameMPois(y=sal.n27[,,1], type="removal")

fm1.1 <- multinomPois(~1~1, umf1)

exp(0.86)*64

