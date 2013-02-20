PoisGLMM <-
function(y=y,site=site,mu0=mu0,sig0=sig0,a0=a0,b0=b0, mu_beta=mu_beta, sig_beta=sig_beta, 
delta_alpha=delta_alpha, delta_beta=delta_beta, niter=niter){

lev<-length(unique(site))     #number of sites
alpha<-rnorm(lev,0,10)#initial values alpha
beta<-runif(1,0,5)#initial value beta
mu_alpha<-mean(alpha)
var_alpha<-var(alpha)
var0<-sig0^2

out<-matrix(nrow=niter, ncol=3)
colnames(out)<-c('mu_alpha','sig_alpha','beta')

for (iter in 1:niter) {

##########update alpha
#initiate counter for acceptance rate of alpha
alphaUps<-0

#loop over sites, update intercepts alpha one at a time; only data at site i contributes information

for (i in 1:lev) { 
alpha.cand<-rnorm(1, alpha[i], delta_alpha)
loglike<- sum(dpois (y[site==i], exp(alpha[i] + beta*x[site==i]), log=TRUE))  
logprior<- dnorm(alpha[i], mu_alpha, sqrt(var_alpha), log=TRUE)
loglike.cand<- sum(dpois (y[site==i], exp(alpha.cand + beta *x[site==i]), log=TRUE))
logprior.cand<- dnorm(alpha.cand,  mu_alpha,sqrt(var_alpha), log=TRUE)
if (runif(1)< exp((loglike.cand+logprior.cand) -(loglike+logprior))) {
alpha[i]<-alpha.cand
alphaUps<-alphaUps+1
}
}

if(iter %% 100 == 0) {  #this lets you check the acceptance rate of alpha at every 100th iteration
            cat("   Acceptance rates\n")
            cat("     alpha =", alphaUps/lev, "\n")
}

###########update beta
beta.cand<-rnorm(1, beta, delta_beta)
avec<-rep(alpha, times=c(rep(10,10)))
loglike<- sum(dpois (y, exp(avec + beta*x), log=TRUE))  
logprior<- dnorm(beta, mu_beta,sig_beta, log=TRUE)
loglike.cand<- sum(dpois (y, exp(avec + beta.cand *x), log=TRUE))  
logprior.cand<- dunif(beta.cand, mu_beta,sig_beta, log=TRUE)
if (runif(1)< exp((loglike.cand+logprior.cand) - (loglike+logprior) )) {
beta<-beta.cand
}

###########update mu_alpha using Gibbs sampling
abar<-mean(alpha)

mun<- (var_alpha/(var_alpha+lev*var0))*mu0 + (lev*var0/(var_alpha+lev* var0))*abar 

varn <- (var_alpha*var0)/ (var_alpha+lev*var0)
mu_alpha<-rnorm(1,mun, sqrt(varn))

#update var_alpha using Gibbs sampling
an<-lev/2 + a0
bn<- 0.5 * (  sum(   (alpha-mu_alpha)^2  )   ) +b0

var_alpha<-1/rgamma(1,shape=an, rate=bn)

out[iter,]<-c(mu_alpha, sqrt(var_alpha), beta)

}

return(out)

}

