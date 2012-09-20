\name{SCRovenbird}
\alias{SCRovenbird}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fit SCR model to ovenbird data using JAGS%%  ~~function to do ... ~~
}
\description{
yes, does that    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRovenbird()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

    library("scrbook")
    library("secr")
    data(ovenbird)
    set.seed(2013)
    png("Ovenbird_traps.png", width = 7, height = 7, units = "in", 
        res = 400)
    par(mfrow = c(1, 3))
    plot(ovenCH[["2005"]])
    plot(ovenCH[["2007"]])
    plot(ovenCH[["2009"]])
    dev.off()
    X <- traps <- traps(ovenCH)
    xlim <- c(min(X[[1]][, 1]) - 150, max(X[[1]][, 1]) + 150)
    ylim <- c(min(X[[1]][, 2]) - 150, max(X[[1]][, 2]) + 150)
    ntraps <- nrow(traps[[1]])
    Y <- ovenCH
    K <- 10
    M <- 50
    Sst <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], 
        ylim[2]))
    Sst <- array(Sst, dim = c(M, 2, 5))
    Xb.any <- array(0, dim = c(M, K, 5))
    Xb.trap <- array(0, dim = c(M, K, ntraps, 5))
    Yarr <- array(ntraps + 1, dim = c(M, K, 5))
    for (i in 1:5) {
        tmp <- Y[[i]]
        nind <- nrow(tmp)
        nrep <- ncol(tmp)
        first <- rep(NA, nrow(tmp))
        tmp[tmp < 0] <- tmp[tmp < 0] * (-1)
        tmp[tmp == 0] <- ntraps + 1
        for (j in 1:nind) {
            first[j] <- (1:nrep)[tmp[j, ] < (ntraps + 1)][1]
            if (first[j] == nrep) 
                next
            Xb.any[j, (first[j] + 1):K, i] <- 1
            Xb.trap[j, (first[j] + 1):K, tmp[j, first[j]], i] <- 1
        }
        tmp2 <- matrix(ntraps + 1, nrow = M, ncol = 10)
        tmp2[1:nind, 1:nrep] <- tmp
        Stmp <- Sst[, , i]
        sout <- spiderplot(tmp2[1:nind, 1:nrep], as.matrix(X[[i]]))$avg.s
        browser()
        Stmp[1:nind, 1:2] <- sout
        Sst[, , i] <- Stmp
        Yarr[, , i] <- tmp2
    }


## 
## This bit of code dumps out the BUGS model file
##
cat("
model {
 # year-specific N parameterized in DA parameter 
for(t in 1:5){
psi[t] ~ dunif(0,1)
}
alpha0 ~ dnorm(0,.1)
sigma ~dunif(0,200)
alpha1<- 1/(2*sigma*sigma)
alpha2 ~ dnorm(0,.1)
alpha3~  dnorm(0,.1)

A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
for(t in 1:5){
N[t] <- sum(z[1:M,t]) 
D[t] <- (N[t]/A)*10000  # put in units of per ha

for(i in 1:M){
  z[i,t] ~ dbern(psi[t])
  S[i,1,t] ~ dunif(xlim[1],xlim[2])
  S[i,2,t] ~ dunif(ylim[1],ylim[2])
  for(j in 1:ntraps){
    #distance from capture to the center of the home range
    d2[i,j,t] <- pow(pow(S[i,1,t]-X[j,1],2) + pow(S[i,2,t]-X[j,2],2),1)
  }
  for(k in 1:K){
   for(j in 1:ntraps){
    bleen1[i,k,j,t]<- alpha2*Xb.any[i,k,t] 
    bleen2[i,k,j,t]<- alpha3*Xb.trap[i,k,j,t]
    lp[i,k,j,t] <- exp(alpha0 - alpha1*d2[i,j,t] + 0)*z[i,t]            
    cp[i,k,j,t] <- lp[i,k,j,t]/(1+sum(lp[i,k,1:ntraps,t]))
    }
    cp[i,k,ntraps+1,t] <- 1-sum(cp[i,k,1:ntraps,t])  # last cell = not captured
    Ycat[i,k,t] ~ dcat(cp[i,k,,t])
  }  
 }   
}
   
}
",file="model.txt")
###
###
###

nind <- c(nrow(Y[[1]]), nrow(Y[[2]]), nrow(Y[[3]]), nrow(Y[[4]]), 
        nrow(Y[[5]]))
    zst <- matrix(NA, nrow = M, ncol = 5)
    for (i in 1:5) {
        zst[, i] <- c(rep(1, nind[i]), rep(0, M - nind[i]))
    }
    inits <- function() {
        list(z = zst, sigma = runif(1, 50, 100), S = Sst, alpha0 = runif(1, 
            -2, -1), alpha2 = -2, alpha3 = -2)
    }
    parameters <- c("psi", "alpha0", "alpha1", "alpha2", "alpha3", 
        "sigma", "N", "D")
    data <- list(X = as.matrix(X[[1]]), K = 10, Ycat = Yarr, 
        M = M, ntraps = ntraps, ylim = ylim, xlim = xlim, Xb.any = Xb.any, 
        Xb.trap = Xb.trap)
    library("R2jags")
    out <- jags(data, inits, parameters, "model.txt", n.thin = 1, 
        n.chains = 3, n.burnin = 1000, n.iter = 2000, DIC = FALSE)
    out

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
