
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


# Coordinates of traps
# In this case, every location in the state-space has a trap, so the number
# of places in the state-space G is equal to the number of trap, i.e. G=J
X <- data.matrix(pt.data[,c("coords.x1", "coords.x2")])
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

dat1 <- list(n=sal.n, distmat=distmat, seg=seg, PrSeg=rep(1/G, G),
             J=J, G=G, K=3, M=500)
str(dat1)

n.jk <- apply(dat1$n, c(1,3), sum, na.rm=TRUE)
xxx <- rowSums(n.jk)/sum(n.jk)

set.seed(53450)
si <- sample(seg, dat1$M, replace=TRUE, prob=xxx)
ui <- matrix(NA, dat1$M, 3)
yi <- array(0L, c(dat1$M, dat1$G, 3))

#for(i in 1:dat1$M) {
#    PrU <- exp(-distmat[si[i],]^2/(2*800))
#    ui[i,] <- sample(dat1$seg, 3, PrU, replace=TRUE)
#    for(g in 1:dat1$G) {
#        yi[i,g,] <- as.integer(ui[i,] == g)
#    }
#}



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

ui <- apply(yi, c(1,3), function(x) {
    if(any(x > 0))
        return(which(x==1))
    else
        return(sample(1:dat1$M, 1))
})

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





#for(i in 1:dat1$M) {
#    for(k in 1:3) {
#        g <- which(yi[i,,k] == 1)
#        ui[i,k] <- which(


init1 <- function() list(p=runif(1), psi=runif(1),
                         tau=runif(1, 4000, 5000),
                         s=si,
                         u=ui,
#                         y=yi,
                         z=rep(1, dat1$M))
str(init1())


jm1 <- jags.model("sal1.jag", data=dat1, inits=init1, n.chains=1,
                  n.adapt=100)
