\name{PoisGLMM}
\alias{PoisGLMM}
\title{
MCMC algorithm for Poisson regression with random intercept
}
\description{
MCMC algorithm for Poisson regression with random intercept
}
\usage{
PoisGLMM(y,site,mu0,sig0,a0,b0, mu_beta, sig_beta, delta_alpha, delta_beta, niter)
}
\arguments{
  \item{y}{
	Vector with observed counts
}
  \item{site}{
    vector of length(y) relating each count to a site
}
  \item{mu0}{
    Mean of normal hyperprior for mean(alpha)
}
  \item{sig0}{
    Standard deviation of normal hyperprior for mean(alpha)
}
  \item{a0}{
    Shape parameter for inverse-gamma hyperprior on var(alpha)
}
  \item{b0}{
    Scale parameter for inverse-gamma hyperprior on var(alpha)    
}
  \item{mu_beta}{
    Mean of normal prior for beta
}
  \item{sig_beta}{
    Standard deviation of normal prior for beta
}
  \item{delta_alpha}{
    Tuning parameter for alpha
}
  \item{delta_beta}{
    Tuning parameter for beta
}
  \item{niter}{
    Number of iterations to run the algorithm for
}



}
\value{
  A matrix with one row per iteration and three columns, for the estimated mean and stadard deviation of alpha (mu_alpha and sig_alpha),
  and estimates of beta.
}
\details{
  Starting values are chosen at random within the function:
S<-cbind(runif(M, xl, xu), runif(M,yl,yu))
sigma<-runif(1,0.5,5)
lam0<-runif(1,0.1,1)
psi<-runif(1,0.2,0.8)
z<-rbinom(M,1,psi)
}
\author{
Rahel Sollmann
}
\examples{
n<-10 #number plots per site
i=10 #number sites
set.seed(2013)

#true hyperparameters for random intercept alpha
mua.t<-2
siga.t<-0.1

#random site specific intercepts
alpha=rep(NA, i) 
for (i in 1:n){
alpha[i]=rnorm(1,mua.t,siga.t)
}

#true value for beta
beta<-1.2

#generate plot-level covariate x
x<-rnorm(n*i,0,1)

#generate Poisson mean
ybar=site=matrix(nrow=10, ncol=10)
for (i in 1:n){
ybar[i,]<-exp(alpha[i]+beta*x[c((i*10-9) : (i*10))])
site[i,]<-i
}

#data
ybar=as.vector(t(ybar))
site=as.factor(as.vector(t(site)))
y<-rpois(n*i,ybar)

#set hyperparameters for mua (Normal) and siga (inverse-gamma)
mu0=0
sig0=100
a0=0.01
b0=0.01

#set prior mean and SD for beta
mu_beta=0
sig_beta=100

#tuning parameters for individual alpha and beta
delta_alpha=0.1
delta_beta=0.1

#run Metropolis within Gibbs sampler
Pmod<-PoisGLMM(y,site,mu0, sig0, a0, b0, mu_beta, sig_beta,delta_alpha, delta_beta, 5000)

#use coda to look at output
library(coda)
plot(mcmc(Pmod))
rejectionRate(mcmc(Pmod))
summary(window(mcmc(Pmod), start=1001))
}
