\name{wolvSCR0ms.fn}
\alias{wolvSCR0ms.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wolvSCR0ms.fn(y3d, traps, wsex, nb = 1000, ni = 2000, delta = 2, M = 200, model = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{
%%     ~~Describe \code{y3d} here~~
}
  \item{traps}{
%%     ~~Describe \code{traps} here~~
}
  \item{wsex}{
%%     ~~Describe \code{wsex} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
}
  \item{ni}{
%%     ~~Describe \code{ni} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
}
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

## The function is currently defined as
function (y3d, traps, wsex, nb = 1000, ni = 2000, delta = 2, 
    M = 200, model = 1) 
{
    library("R2WinBUGS")
    traplocs <- as.matrix(traps[, 2:3])
    mingridx <- min(traplocs[, 1])
    mingridy <- min(traplocs[, 2])
    traplocs[, 1] <- traplocs[, 1] - min(traplocs[, 1])
    traplocs[, 2] <- traplocs[, 2] - min(traplocs[, 2])
    traplocs <- traplocs/10000
    ntraps <- nrow(traplocs)
    Xl <- min(traplocs[, 1] - delta)
    Xu <- max(traplocs[, 1] + delta)
    Yl <- min(traplocs[, 2] - delta)
    Yu <- max(traplocs[, 2] + delta)
    area <- (Xu - Xl) * (Yu - Yl)/10
    plot(traplocs, pch = 20, xlim = c(Xl, Xu), ylim = c(Yl, Yu))
    nz <- M - dim(y3d)[1]
    wsex <- c(wsex, rep(NA, nz))
    MASK <- traps[, 4:ncol(traps)]
    Dmat <- as.matrix(dist(traplocs))
    nind <- dim(y3d)[1]
    K <- dim(y3d)[2]
    newy <- array(0, dim = c(nind + nz, K, ntraps))
    for (j in 1:nind) {
        newy[j, 1:K, 1:ntraps] <- y3d[j, 1:K, 1:ntraps]
    }
    y3d <- newy
    M <- nind + nz
    K <- apply(MASK, 1, sum)
    ncaps <- apply(y3d, c(1, 3), sum)
    sink("modelfile0.txt")
    cat("\nmodel {\nalpha.sex ~ dunif(-3,3)\nbeta.sex  ~ dunif(-3,3)\nsigma~dunif(0,50)\nlogitp0~dnorm(0,.1)\nbeta<- (1/(2*sigma*sigma) )\npsi ~ dunif(0,1)\npsi.sex  ~ dunif(0,1)\n\nfor(i in 1:M){\n wsex[i] ~ dbern(psi.sex)\n w[i]~dbern(psi)\n s[i,1]~dunif(Xl,Xu)\n s[i,2]~dunif(Yl,Yu) \n logit(p0[i])<- logitp0 +  mod[1]*alpha.sex*wsex[i]\nlog(beta.vec[i])<- log( beta ) + mod[2]*wsex[i]*beta.sex\nfor(j in 1:ntraps){\n   mu[i,j]<-w[i]*p[i,j]\nmu2[i,j]<- mu[i,j]*K[j]\n ncaps[i,j]~ dbin(mu[i,j],K[j]) \n dd[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2) \n \n p[i,j]  <-  p0[i]*exp( - beta.vec[i]*dd[i,j] )\nncapsnew[i,j]~dbin(mu[i,j],K[j])\nerr[i,j]<-  pow(pow(ncaps[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\nerrnew[i,j]<- pow(pow(ncapsnew[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\n}\n\n\nexpected[i]<- sum(mu2[i,])\nnsum[i]<- sum(ncaps[i,])\nnsumnew[i] <- sum(ncapsnew[i,])\n\nerr1[i]<- pow( pow(nsum[i],.5)  - pow(expected[i],0.5),2)\nerr1new[i]<- pow( pow(nsumnew[i],.5) - pow(expected[i],0.5),2)\n}\nfor(j in 1:ntraps){\ntraptotals[j]<- sum(ncaps[1:M,j])\ntraptotalsnew[j]<-sum(ncapsnew[1:M,j])\nexpectedtrap[j]<- sum( mu[,j]) *K[j]\nerr3[j]<- pow( pow(traptotals[j],.5) - pow(expectedtrap[j],.5),2)\nerr3new[j]<- pow( pow(traptotalsnew[j],.5) - pow(expectedtrap[j],.5),2)\n}\nX3obs<- sum(err3[])\nX3new<- sum(err3new[])\nX1obs<- sum(err1[])\nX1new<- sum(err1new[])\nXobs<-sum(err[,])\nXnew<-sum(errnew[,])\n\nN<-sum(w[1:M])\nD<-N/area\n}\n", 
        fill = TRUE)
    sink()
    sink("modelfile5.txt")
    cat("\nmodel {\nalpha.sex ~ dunif(-3,3)\nbeta.sex  ~ dunif(-3,3)\nbeta0~dunif(0,10)\n\n\nlogitp0~dnorm(0,.1)\npsi ~ dunif(0,1)\npsi.sex  ~ dunif(0,1)\nmod[1] ~ dbern(.5)\nmod[2] ~ dbern(.5)\nfor(i in 1:M){\n wsex[i] ~ dbern(psi.sex)\n w[i]~dbern(psi)\n s[i,1]~dunif(Xl,Xu)\n s[i,2]~dunif(Yl,Yu) \n logit(p0[i])<- logitp0 +  mod[1]*frog[i]\nfrog[i]<- alpha.sex*wsex[i]\n\nbeta.vec[i]<- beta0 + mod[2]*wsex[i]*beta.sex\n\nfor(j in 1:ntraps){\n   mu[i,j]<-w[i]*p[i,j]\nmu2[i,j]<- mu[i,j]*K[j]\n ncaps[i,j]~ dbin(mu[i,j],K[j]) \n dd[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2) \n \n p[i,j]  <-  p0[i]*exp( - beta.vec[i]*dd[i,j] )\nncapsnew[i,j]~dbin(mu[i,j],K[j])\nerr[i,j]<-  pow(pow(ncaps[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\nerrnew[i,j]<- pow(pow(ncapsnew[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\n}\n\n\nexpected[i]<- sum(mu2[i,])\nnsum[i]<- sum(ncaps[i,])\nnsumnew[i] <- sum(ncapsnew[i,])\n\nerr1[i]<- pow( pow(nsum[i],.5)  - pow(expected[i],0.5),2)\nerr1new[i]<- pow( pow(nsumnew[i],.5) - pow(expected[i],0.5),2)\n}\nfor(j in 1:ntraps){\ntraptotals[j]<- sum(ncaps[1:M,j])\ntraptotalsnew[j]<-sum(ncapsnew[1:M,j])\nexpectedtrap[j]<- sum( mu[,j]) *K[j]\nerr3[j]<- pow( pow(traptotals[j],.5) - pow(expectedtrap[j],.5),2)\nerr3new[j]<- pow( pow(traptotalsnew[j],.5) - pow(expectedtrap[j],.5),2)\n}\nX3obs<- sum(err3[])\nX3new<- sum(err3new[])\nX1obs<- sum(err1[])\nX1new<- sum(err1new[])\nXobs<-sum(err[,])\nXnew<-sum(errnew[,])\n\nN<-sum(w[1:M])\nD<-N/area\n}\n", 
        fill = TRUE)
    sink()
    sink("modelfile5b.txt")
    cat("\nmodel {\nalpha.sex ~ dunif(-3,3)\nbeta.sex  ~ dunif(-3,3)\nsigma0~dunif(0,10)\n\nlogitp0~dnorm(0,.1)\npsi ~ dunif(0,1)\npsi.sex  ~ dunif(0,1)\nmod[1] ~ dbern(.5)\nmod[2] ~ dbern(.5)\nfor(i in 1:M){\n wsex[i] ~ dbern(psi.sex)\n w[i]~dbern(psi)\n s[i,1]~dunif(Xl,Xu)\n s[i,2]~dunif(Yl,Yu) \n logit(p0[i])<- logitp0 +  mod[1]*frog[i]\nfrog[i]<- alpha.sex*wsex[i]\n\nlog(sigma[i])<- log(sigma0) +mod[2]*wsex[i]*beta.sex\nbeta.vec[i]<- 1/(2*sigma[i]*sigma[i])\n\n\nfor(j in 1:ntraps){\n   mu[i,j]<-w[i]*p[i,j]\nmu2[i,j]<- mu[i,j]*K[j]\n ncaps[i,j]~ dbin(mu[i,j],K[j]) \n dd[i,j]<- pow(s[i,1] - traplocs[j,1],2)  + pow(s[i,2] - traplocs[j,2],2) \n \n p[i,j]  <-  p0[i]*exp( - beta.vec[i]*dd[i,j] )\nncapsnew[i,j]~dbin(mu[i,j],K[j])\nerr[i,j]<-  pow(pow(ncaps[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\nerrnew[i,j]<- pow(pow(ncapsnew[i,j],.5) - pow(K[j]*mu[i,j],.5),2)\n}\n\n\nexpected[i]<- sum(mu2[i,])\nnsum[i]<- sum(ncaps[i,])\nnsumnew[i] <- sum(ncapsnew[i,])\n\nerr1[i]<- pow( pow(nsum[i],.5)  - pow(expected[i],0.5),2)\nerr1new[i]<- pow( pow(nsumnew[i],.5) - pow(expected[i],0.5),2)\n}\nfor(j in 1:ntraps){\ntraptotals[j]<- sum(ncaps[1:M,j])\ntraptotalsnew[j]<-sum(ncapsnew[1:M,j])\nexpectedtrap[j]<- sum( mu[,j]) *K[j]\nerr3[j]<- pow( pow(traptotals[j],.5) - pow(expectedtrap[j],.5),2)\nerr3new[j]<- pow( pow(traptotalsnew[j],.5) - pow(expectedtrap[j],.5),2)\n}\nX3obs<- sum(err3[])\nX3new<- sum(err3new[])\nX1obs<- sum(err1[])\nX1new<- sum(err1new[])\nXobs<-sum(err[,])\nXnew<-sum(errnew[,])\n\nN<-sum(w[1:M])\nD<-N/area\n}\n", 
        fill = TRUE)
    sink()
    sst <- sample(1:nrow(traplocs), M, replace = TRUE)
    wst <- c(rep(1, nind), rep(0, M - nind))
    if (model == 5) {
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex")
        inits <- function() {
            list(beta0 = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0, 
                mod = c(1, 1))
        }
        parameters <- c("mod", "beta.sex", "alpha.sex", "psi", 
            "psi.sex", "beta0", "logitp0", "N", "D", "Xobs", 
            "Xnew", "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile5.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    if (model == 6) {
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex")
        inits <- function() {
            list(sigma0 = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0, 
                mod = c(1, 1))
        }
        parameters <- c("mod", "beta.sex", "alpha.sex", "psi", 
            "psi.sex", "sigma0", "logitp0", "N", "D", "Xobs", 
            "Xnew", "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile5b.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    if (model == 1) {
        mod <- c(0, 0)
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex", "mod")
        inits <- function() {
            list(sigma = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0)
        }
        parameters <- c("beta.sex", "alpha.sex", "psi", "psi.sex", 
            "sigma", "logitp0", "N", "D", "beta", "Xobs", "Xnew", 
            "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile0.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    if (model == 2) {
        mod <- c(1, 0)
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex", "mod")
        inits <- function() {
            list(sigma = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0)
        }
        parameters <- c("beta.sex", "alpha.sex", "psi", "psi.sex", 
            "sigma", "logitp0", "N", "D", "beta", "Xobs", "Xnew", 
            "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile0.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    if (model == 3) {
        mod <- c(0, 1)
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex", "mod")
        inits <- function() {
            list(sigma = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0)
        }
        parameters <- c("beta.sex", "alpha.sex", "psi", "psi.sex", 
            "sigma", "logitp0", "N", "D", "beta", "Xobs", "Xnew", 
            "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile0.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    if (model == 4) {
        mod <- c(1, 1)
        data <- list("ncaps", "traplocs", "M", "ntraps", "K", 
            "Xl", "Xu", "Yl", "Yu", "area", "wsex", "mod")
        inits <- function() {
            list(sigma = runif(1, 0.4, 1), logitp0 = rnorm(1, 
                -2, 0.2), w = wst, alpha.sex = 0, beta.sex = 0)
        }
        parameters <- c("beta.sex", "alpha.sex", "psi", "psi.sex", 
            "sigma", "logitp0", "N", "D", "beta", "Xobs", "Xnew", 
            "X1obs", "X1new", "X3obs", "X3new")
        out <- jags(data, inits, parameters, "modelfile0.txt", 
            n.thin = 1, n.chains = 3, n.burnin = nb, n.iter = ni, 
            working.dir = getwd())
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
