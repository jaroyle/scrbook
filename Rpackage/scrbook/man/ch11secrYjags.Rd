\name{ch11secrYjags}
\alias{ch11secrYjags}
\title{
  Fit IPP using secr and JAGS
}
\description{
  The analysis of is of simulated data with a single state-space
  covariate. See Chapter 11 for an application using real data.
}
\details{
  See Chapter 11?
}
\seealso{
  \code{\link{scrIPP}}
}
\examples{

\dontrun{

library(secr)
library(rjags)

data(ch11simData)

ch <- ch11simData$ch.secr
msk <- ch11simData$spcov.secr


# SECR analysis

secr1 <- secr.fit(ch, model=D~canht, mask=msk)

region.N(secr1, se.N=TRUE)





# JAGS analysis

# JAGS model
sink("ippDiscrete.txt")
cat("
model{
sigma ~ dunif(0, 20)
lam0 ~ dunif(0, 5)
beta0 ~ dunif(-10, 10)
beta1 ~ dunif(-10, 10)
for(j in 1:nPix) {
  mu[j] <- exp(beta0 + beta1*CANHT[j])*pixArea
  probs[j] <- mu[j]/EN
}
EN <- sum(mu[])
psi <- EN/M
for(i in 1:M) {
  w[i] ~ dbern(psi)
  s[i] ~ dcat(probs[])
  x0g[i] <- Sgrid[s[i],1]
  y0g[i] <- Sgrid[s[i],2]
  for(j in 1:ntraps) {
    dist[i,j] <- sqrt((Sgrid[i,1]-traps[j,1])^2 +
                     (Sgrid[i,2]-traps[j,2])^2)
    lambda[i,j] <- lam0*exp(-dist[i,j]*dist[i,j]/(2*sigma*sigma)) * w[i]
    y[i,j] ~ dpois(lambda[i,j])
    }
  }
N <- sum(w[])
D <- N/1 # 1ha state-space
}
", fill=TRUE)
sink()



modfile <- "ippDiscrete.txt"

jags.data <- with(ch11simData, {
    list(y=ch.jags, CANHT=drop(spcov.jags$CANHT),
            nPix=nrow(spcov.jags),
            M=nrow(ch.jags), ntraps=nrow(traps),
            Sgrid=as.matrix(spcov.jags[,1:2]),
            pixArea=25,
            traps=traps)
    })
str(jags.data)

init <- function() {
    list(sigma=runif(1, 10, 11), lam0=runif(1),
         beta0=-8, beta1=1, #rnorm(1, 1),
         s=sample.int(jags.data$nPix, jags.data$M, replace=TRUE),
         w=rep(1, jags.data$M))
}
str(init())

pars <- c("sigma", "lam0", "beta0", "beta1", "N")

# Obtain posterior samples. This takes a few minutes
# Compile and adapt
set.seed(03453)
jm <- jags.model(modfile, jags.data, init, n.chains=2, n.adapt=1000)
# MCMC
jags1 <- coda.samples(jm, pars, n.iter=2000)

plot(jags1, ask=TRUE)
summary(window(jags1, start=1001))


unlink(modfile)


}



}
