\name{snakeline}
\alias{snakeline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 "Search-encounter" data from sampling along a line
}
\description{
 R script to simulate individual "search-encounter" data by sampling along a line,
  and fit a model using JAGS and the R2jags package.
}
\usage{
execute line-by-line in an R console
}
%- maybe also 'usage' for other objects documented here.
\details{
See chapter 15 %%  ~~ If necessary, more details than the description above ~~
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
Andy Royle, aroyle@usgs.gov 
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

    basex <- c(0, 0, 1, 1, 0)
    basey <- c(0, 1, 1, 0, 0)
    base <- cbind(basex, basey)
    plot(basex, basey, pch = " ")
    polygon(base)
    p1 <- cbind(basex, basey + 3)
    p2 <- cbind(basex + 1, basey + 3)
    p3 <- cbind(basex + 1, basey + 2)
    p4 <- cbind(basex + 2, basey + 2)
    p5 <- cbind(basex + 1, basey + 1)
    p6 <- cbind(basex + 2, basey + 1)
    p7 <- cbind(basex + 2, basey)
    pp <- rbind(p1, p2, p3, p4, p5, p6, p7, base)


    plot(pp, xlab = "Easting", ylab = "Northing", pch = " ", 
        cex.axis = 1.5, cex = 2, cex.lab = 1.5)
    polygon(p1)
    polygon(p2)
    polygon(p3)
    polygon(p4)
    polygon(p5)
    polygon(p6)
    polygon(p7)

line1 <-
structure(list(x = c(0.139453125, 0.315234375, 0.4734375, 0.64921875, 
0.76640625, 0.866015625, 1.006640625, 1.100390625, 1.217578125, 
1.340625, 1.46953125, 1.639453125, 1.70390625, 1.79765625, 1.903125, 
1.932421875, 1.932421875, 1.932421875, 1.95, 2.026171875, 2.1609375, 
2.272265625, 2.36015625, 2.5125, 2.600390625, 2.75859375, 2.840625, 
2.87578125, 2.79375, 2.68828125, 2.47734375, 2.36015625, 2.184375, 
2.073046875, 1.897265625, 1.75078125, 1.569140625, 1.4578125, 
1.264453125, 1.28203125, 1.440234375, 1.5984375, 1.84453125, 
1.932421875, 2.014453125, 2.073046875, 2.1375, 2.190234375, 2.325, 
2.45390625), y = c(3.16135609677158, 3.32738894732584, 3.49342179788009, 
3.68567036167975, 3.73810178817056, 3.77305607249777, 3.77305607249777, 
3.73810178817056, 3.67693179059794, 3.61576179302533, 3.5371146532891, 
3.40603608706206, 3.35360466057124, 3.26621894975322, 3.02153895946274, 
2.86424467999029, 2.75064325592685, 2.60208754753621, 2.51470183671818, 
2.46227041022736, 2.40110041265475, 2.39236184157294, 2.38362327049114, 
2.37488469940934, 2.35740755724573, 2.2612832753459, 2.07777328262804, 
1.85057043450117, 1.72823043935593, 1.64084472853791, 1.58841330204709, 
1.6233675863743, 1.6321061574561, 1.64084472853791, 1.6146290152925, 
1.58841330204709, 1.57093615988349, 1.54472044663808, 1.43985759365644, 
1.34373331175661, 1.26508617202039, 1.24760902985678, 1.22139331661138, 
1.16896189012056, 1.08157617930253, 1.02040618172991, 0.915543328748281, 
0.836896189012056, 0.714556193866818, 0.583477627639777)), .Names = c("x", 
"y"))




###    line1 <- source("line1.R")
    line1 <- cbind(line1$x, line1$y)
    line1 <- as.matrix(line1)
    lines(line1, lwd = 2)




    library(rgeos)
    library(sp)


    points <- SpatialPoints(line1)
    sLine <- Line(points)
    regpoints <- spsample(sLine, 100, type = "regular")
    plot(line1, type = "l")
    points(points, col = "grey")
    points(regpoints, col = "red", pch = 20, lwd = 2)




    perbox <- 4
    N <- 30 * perbox
    xlim <- c(-1, 4)
    ylim <- c(-1, 5)
    set.seed(2014)
    sx <- runif(N, xlim[1], xlim[2])
    sy <- runif(N, ylim[1], ylim[2])
    points(sx, sy, pch = 20, col = "red")
    sigma.move <- 0.35
    sigma <- 0.4
    alpha0 <- 0.8
    alpha1 <- 1/(2 * (sigma^2))
    X <- regpoints@coords
    J <- nrow(X)
    K <- 10
    U <- array(NA, dim = c(N, K, 2))
    y <- pmat <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N) {
        for (k in 1:K) {
            U[i, k, ] <- c(rnorm(1, sx[i], sigma.move), rnorm(1, 
                sy[i], sigma.move))
            dvec <- sqrt((U[i, k, 1] - X[, 1])^2 + (U[i, k, 2] - 
                X[, 2])^2)
            loghaz <- alpha0 - alpha1 * dvec * dvec
            H <- sum(exp(loghaz))
            pmat[i, k] <- 1 - exp(-H)
            y[i, k] <- rbinom(1, 1, pmat[i, k])
        }
    }
    Ux <- U[, , 1]
    Uy <- U[, , 2]
    Ux[y == 0] <- NA
    Uy[y == 0] <- NA
    points(Ux, Uy, pch = 20, col = "black")
    ncap <- apply(y, 1, sum)
    y <- y[ncap > 0, ]
    Ux <- Ux[ncap > 0, ]
    Uy <- Uy[ncap > 0, ]
    M <- 200
    nind <- nrow(y)
    y <- rbind(y, matrix(0, nrow = (M - nrow(y)), ncol = ncol(y)))
    Namat <- matrix(NA, nrow = (M - nind), ncol = ncol(y))
    Ux <- rbind(Ux, Namat)
    Uy <- rbind(Uy, Namat)
    S <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], 
        ylim[2]))
    for (i in 1:nind) {
        S[i, ] <- c(mean(Ux[i, ], na.rm = TRUE), mean(Uy[i, ], 
            na.rm = TRUE))
    }
    Ux.st <- Ux
    Uy.st <- Uy
    for (i in 1:M) {
        Ux.st[i, !is.na(Ux[i, ])] <- NA
        Uy.st[i, !is.na(Uy[i, ])] <- NA
        Ux.st[i, is.na(Ux[i, ])] <- S[i, 1]
        Uy.st[i, is.na(Uy[i, ])] <- S[i, 2]
    }
    cat("\nmodel {\n\n# Priors\nalpha0~dunif(-25,25)\nalpha1~dunif(0,60)\nsigma<- sqrt(1/(2*alpha1))\nlsigma~dunif(-5,5)\nsigma.move<-exp(lsigma)\ntau<-1/(sigma.move*sigma.move)\npsi~dunif(0,1)\n\n# Likelihood\nfor(i in 1:M){ # Loop over individuals\n z[i]~dbern(psi)\n s[i,1]~dunif(xlim[1],xlim[2])\n s[i,2]~dunif(ylim[1],ylim[2])\n for(k in 1:K){ # Loop over temporal replicates\n    u[i,k] ~ dnorm(s[i,1],tau) \n    v[i,k] ~ dnorm(s[i,2],tau) \n    for(j in 1:J){ # Loop over each point defining line segments\n      d[i,k,j]<-  pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)\n      h[i,k,j]<-exp(alpha0-alpha1*d[i,k,j]*d[i,k,j])\n   }\n   H[i,k]<-sum(h[i,k,1:J])\n   p[i,k]<- z[i]*(1-exp(-H[i,k]))\n   y[i,k] ~ dbern(p[i,k])\n }\n}\n\n# Derived quantity\nN<-sum(z[])\n}\n", 
        file = "model0.txt")
    data <- list(y = y, u = Ux, v = Uy, X = X, K = K, M = M, 
        J = J, xlim = xlim, ylim = ylim)
    inits <- function() {
        list(alpha0 = alpha0 - 0.3, alpha1 = alpha1 - 1.5, lsigma = log(0.5), 
            s = S, z = c(rep(1, nind), rep(0, M - nind)), u = Ux.st, 
            v = Uy.st)
    }
    parameters <- c("alpha0", "alpha1", "N", "psi", "sigma.move", 
        "sigma")
    nthin <- 1
    nc <- 3
    nb <- 500
    ni <- 3500
    library("R2jags")
    wbout2 <- jags(data, inits, parameters, "model0.txt", n.thin = nthin, 
        n.chains = nc, n.burnin = nb, n.iter = ni, working.dir = getwd())
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
