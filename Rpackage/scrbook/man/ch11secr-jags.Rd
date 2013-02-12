\name{ch11secr-jags}
\alias{ch11secr-jags}
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
  z[i] ~ dbern(psi)
  s[i] ~ dcat(probs[])
  x0g[i] <- Sgrid[s[i],1]
  y0g[i] <- Sgrid[s[i],2]
  for(j in 1:ntraps) {
    dist[i,j] <- sqrt((x0g[i]-traps[j,1])^2 +
                     (y0g[i]-traps[j,2])^2)
    lambda[i,j] <- lam0*exp(-dist[i,j]^2/(2*sigma^2)) * z[i]
    y[i,j] ~ dpois(lambda[i,j])
    }
  }
N <- sum(z[])
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
         z=rep(1, jags.data$M))
}
str(init())

pars <- c("sigma", "lam0", "beta0", "beta1", "N", "EN")

# Obtain posterior samples. This takes a few minutes
# Compile and adapt
set.seed(03453)
jm <- jags.model(modfile, jags.data, init, n.chains=2, n.adapt=1000)
# MCMC
jags1 <- coda.samples(jm, pars, n.iter=6000)

plot(jags1, ask=TRUE)
summary(window(jags1, start=1001))


unlink(modfile)




jags.est <- summary(jags1)
jags.r <- cbind(jags.est$stat[,1:2], jags.est$quant[,c(1,5)])
jags.r

secr.est <- predict(secr1)
secr.r <- cbind(secr.est[2:3,2:5])
secr.r <- rbind(beta=as.numeric(coef(secr1)[2,]), secr.r)
secr.r <- data.matrix(rbind(region.N(secr1)[,1:4], secr.r))
secr.r

jagsVsecr <-
data.frame(Par=rep(c("$\\lambda_0$", "$\\sigma$", "$\\beta_1$", "$N$",
                   "$\\mathbb{E}[N]$"), each=2),
           Truth=rep(c(1, 10, 1, 30, 32.3), each=2),
           Software = rep(c("\\textbf{JAGS}", "\\texttt{secr}"), 5),
           rbind(jags.r[5,], secr.r[4,],
                 jags.r[6,], secr.r[5,],
                 jags.r[4,], secr.r[3,],
                 jags.r[2,], secr.r[2,],
                 jags.r[1,], secr.r[1,]))


format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE)

write.table(format(jagsVsecr, digits=2, nsmall=2, scientific=FALSE),
            file="../../../Ch11-Statespace/R/jagsVsecr2.txt",
            row.names=FALSE,
            quote=FALSE, sep=" \t& ", eol=" \\\\\n ")




}



}
