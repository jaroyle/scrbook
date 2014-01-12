\name{uniform_search}
\alias{uniform_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
R script to fit SCR model based on area search with uniform search
intensity. %%  ~~function to do ... ~~
}
\description{
R script to fit SCR model based on area search with uniform search
intensity, and a spatial autoregressive movement model.
The model is fitted in JAGS using the R2jags package. %%  ~~function to do ... ~~

}
\usage{
uniform_search()
}
%- maybe also 'usage' for other objects documented here.
\details{
  This is an R script disguised as a function. It returns the JAGS
  output but does not take any arguments. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Standard list output from running R2jags. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  SCR book, Chapt. 15

  Royle, J. A., & Young, K. V. (2008). A hierarchical model for spatial capture-recapture data. Ecology, 89(8), 2281-2289.
   (this paper uses a model of "random movement", not autoregressive). 
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



    N <- 100
    nyear <- 4
    Sx <- Sy <- matrix(NA, nrow = N, ncol = nyear)
    sigma.move <- 0.25
    Sx[, 1] <- runif(N, 0, 16)
    Sy[, 1] <- runif(N, 0, 16)
    for (t in 2:nyear) {
        Sx[, t] <- rnorm(N, Sx[, t - 1], sigma.move)
        Sy[, t] <- rnorm(N, Sy[, t - 1], sigma.move)
    }
    f <- rep(NA, N)
    Y <- matrix(0, nrow = N, ncol = nyear)
    for (i in 1:N) {
        for (t in 1:nyear) {
            if (Sx[i, t] > 3 & Sx[i, t] < 13 & Sy[i, t] > 3 & 
                Sy[i, t] < 13) 
                Y[i, t] <- rbinom(1, 1, 0.5)
        }
        if (sum(Y[i, ]) > 0) 
            f[i] <- min((1:nyear)[Y[i, ] == 1][1])
    }
    Y <- Y[!is.na(f), ]
    Sx <- Sx[!is.na(f), ]
    Sy <- Sy[!is.na(f), ]
    Sx[Y == 0] <- NA
    Sy[Y == 0] <- NA


    M <- 200
    Y <- rbind(Y, matrix(0, nrow = (M - nrow(Y)), ncol = nyear))
    Sx <- rbind(Sx, matrix(NA, nrow = (M - nrow(Sx)), ncol = nyear))
    Sy <- rbind(Sy, matrix(NA, nrow = (M - nrow(Sy)), ncol = nyear))
    G <- array(NA, dim = c(M, nyear, 2))
    G[, , 1] <- Sx
    G[, , 2] <- Sy
    sink("model1.mod")
    cat("\nmodel {\npsi ~ dunif(0,1)\ntau ~dgamma(.1,.1)\np0 ~ dunif(0,1)\nsigma.move<- sqrt(1/tau)\n\n\n# Likelihood \nfor (i in 1:M){\nz[i] ~ dbern(psi)\n\n  G[i,1,1] ~ dunif(0,16)\n  G[i,1,2] ~ dunif(0,16)\n\n   for (t in 2:n.occasions){\n   \n## See here I can only make a model for LOCATION\n      G[i,t,1] ~ dnorm(G[i,t-1,1], tau)\n      G[i,t,2] ~ dnorm(G[i,t-1,2], tau)\n      \n# Test whether the actual location is in- or outside the study area. Needs to be done for each grid cell\n    }\n   for(t in 1:n.occasions){ \n      inside[i,t] <- step(G[i,t,1]-3) * step(13-G[i,t,1]) *step(G[i,t,2]-3) * step(13-G[i,t,2])\n      Y[i,t] ~ dbern(mu2[i,t])\n      mu2[i,t] <- p0 * inside[i,t] * z[i]\n      } #t\n   } #i\nN<- sum(z[])\n}\n", 
        fill = TRUE)
    sink()

    ni <- 2000
    nt <- 1
    nb <- 1000
    nc <- 2
    jags.data <- list(Y = Y, G = G, M = M, n.occasions = ncol(Y))
    inits <- function() {
        list(psi = runif(1, 0, 1), p0 = runif(1, 0, 1), z = c(rep(1, 
            N), rep(0, M - N)), tau = 1)
    }
    library("R2jags")
    parameters <- c("p0", "tau", "N", "psi", "sigma.move")
    model1 <- jags(jags.data, inits, parameters, "model1.mod", 
        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
        working.directory = getwd())



  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
