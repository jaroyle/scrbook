##### write function to run bear models in JAGS

bear.JAGS<-function(model=c('SCR0', 'SCR0exp', 'SCRt','SCRB','SCRb', 'SCRsex', 'SCRh'), n.chains, n.adapt, n.iter) {

mod<-match.arg(model)

if (mod =="SCRB" | mod == "SCRb") warning("SCRb and SCRB: Detection parameter alpha0 reported on logit scale; all other models report p0 (on real scale)")

##data setup
library(rjags)
library(scrbook)


data(beardata)
ymat<-beardata$bearArray
trapmat<-beardata$trapmat
sex<-beardata$sex
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3 ]
ntraps<-dim(beardata$bearArray)[2]
M=650
nz<-M-nind

#create augmented array
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-ymat
y<-apply(Yaug,1:2, sum)

#center the coordinates of the trap matrix
X=as.matrix(cbind((trapmat[,1]- mean(trapmat[,1])), (trapmat[,2]- mean(trapmat[,2]))))

#set up the state-space

Xl=min(trapmat[,1]- mean(trapmat[,1])) - 20
Xu=max(trapmat[,1]- mean(trapmat[,1])) + 20
Yl=min(trapmat[,2]- mean(trapmat[,2])) - 20
Yu=max(trapmat[,2]- mean(trapmat[,2])) + 20
areaX=(Xl-Xu)*(Yl-Yu)

#set up run specifications
n.chains=n.chains 
n.adapt=n.adapt
n.iter=n.iter

#get mean activity centers for observed bears; create initial values for remaining s
toad<-spiderplot(ymat, X)
Sin<-matrix(NA, ncol=2, nrow=M)
Sin[1:nind,]<-toad$avg.s
Sin[(nind+1):M,]<-cbind(runif(nz, Xl, Xu), runif(nz,Yl,Yu))


#######################################################################################################
###### SCR0 #############################################################################################

if(mod=='SCR0') {

cat("
model {
alpha0~dnorm(0,.1)
logit(p0)<- alpha0
alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dbin(p[i,j],K)
p[i,j]<- z[i]*p0*exp(- alpha1*d[i,j]*d[i,j])
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCR0.txt")

data<-list(y=y,M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','p0','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, 
		sigma=runif(1,2,3),alpha0=runif(1)) }

modelFile= "SCR0.txt"

} #end if

######################################################################################################
###### SCR0exp #############################################################################################

if(mod=='SCR0exp') {

cat("
model {
alpha0~dnorm(0,.1)
logit(p0)<- alpha0
alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dbin(p[i,j],K)
p[i,j]<- z[i]*p0*exp(- alpha1*d[i,j])
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCR0exp.txt")

data<-list(y=y,M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','p0','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, 
		sigma=runif(1,2,3),alpha0=runif(1)) }

modelFile= "SCR0exp.txt"

} #end if

######################################################################################################
###### SCRt #############################################################################################

if(mod=='SCRt') {

cat("
model {

for(k in 1:K){
alpha0[k]~dnorm(0,.1)
logit(p0[k])<- alpha0[k]
}

alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)

for(k in 1:K){
y[i,j,k] ~ dbin(p[i,j,k],1)
p[i,j,k]<- z[i]*p0[k]*exp(- alpha1*d[i,j]*d[i,j])
}
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCRt.txt")

data<-list(y=Yaug,M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','p0','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)), psi=runif(1), s=Sin,
		 sigma=runif(1,2,3),alpha0=runif(K)) }

modelFile= "SCRt.txt"

}#end if


######################################################################################################
###### SCRsex #############################################################################################

if(mod=='SCRsex') {
cat("
model {

psi~dunif(0,1)
pi~dunif(0,1)

for(t in 1:2){
alpha0[t]~dnorm(0,.1)
logit(p0[t])<- alpha0[t]
alpha1[t]<-1/(2*sigma[t]*sigma[t])
sigma[t]~dunif(0, 15)
}

for(i in 1:M){
 z[i] ~ dbern(psi)
 SEX[i]~dbern(pi)
SEX2[i]<-SEX[i] + 1
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)

for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dbin(p[i,j],K)
p[i,j]<- z[i]*p0[SEX2[i]]*exp(-alpha1[SEX2[i]]*d[i,j]*d[i,j])
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCRsex.txt")

SEX<-c(sex-1, rep(NA, nz))
data<-list(y=y,SEX=SEX, M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','p0','N', 'D', 'sigma', 'pi')
inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, SEX=c(rep(NA, nind), rbinom(nz, 1,0.5)),
		pi=runif(1), sigma=runif(2,2,3),alpha0=runif(2)) }

modelFile= "SCRsex.txt"

}#end if

#######################################################################################################
###### SCRh #############################################################################################

if(mod=='SCRh') {

cat("
model {

alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)
mu_p~dnorm(0,.001)
tau_p~dgamma(.001,.001)

for(i in 1:M){
alpha0[i]~dnorm(mu_p,tau_p)
logit(p0[i])<- alpha0[i]

 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
y[i,j] ~ dbin(p[i,j],K)
p[i,j]<- z[i]*p0[i]*exp(- alpha1*d[i,j]*d[i,j])
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCRh.txt")

data<-list(y=y,M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','mu_p', 'tau_p','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, 
		sigma=runif(1,2,3),alpha0=rnorm(M, 1,1), mu_p=runif(1), tau_p=runif(1)) }

modelFile= "SCRh.txt"

} #end if


######################################################################################################
###### SCRB #############################################################################################

if(mod=='SCRB') {

C=array(0, dim=c(M, ntraps, K))
for(j in 1:ntraps){
for(k in 1:7){
b=which(Yaug[,j,k] > 0)
if(length(b) > 0) {C[b,j,(k+1):8]<-1}
}
}
C<-apply(C, c(1,3), sum)
C[C >1] =1

cat("
model {
alpha0~dnorm(0,.1)
alpha2~dnorm(0,.1)
alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)

for(k in 1:K){
logit(p0[i,j,k])<- alpha0 + alpha2*C[i,k]
y[i,j,k] ~ dbin(p[i,j,k],1)
p[i,j,k]<- z[i]*p0[i,j,k]*exp(- alpha1*d[i,j]*d[i,j])
}
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCRB.txt")


data<-list(y=Yaug, M=M, K=K, C=C, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','alpha0','alpha2','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)), psi=runif(1), s=Sin,
		 sigma=runif(1,2,3),alpha0=runif(1), alpha2=runif(1)) }

modelFile= "SCRB.txt"

}#end if

######################################################################################################
###### SCRb #############################################################################################

if(mod=='SCRb') {

C=array(0, dim=c(M, ntraps, K))
for(j in 1:ntraps){
for(k in 1:7){
b=which(Yaug[,j,k] > 0)
if(length(b) > 0) {C[b,j,(k+1):8]<-1}
}
}

cat("
model {
alpha0~dnorm(0,.1)
alpha2~dnorm(0,.1)
alpha1<-1/(2*sigma*sigma)
sigma~dunif(0, 15)
psi~dunif(0,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)
for(j in 1:J){
d[i,j]<- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)

for(k in 1:K){
logit(p0[i,j,k])<- alpha0 + alpha2*C[i,j,k]
y[i,j,k] ~ dbin(p[i,j,k],1)
p[i,j,k]<- z[i]*p0[i,j,k]*exp(- alpha1*d[i,j]*d[i,j])
}
}
}
N<-sum(z[])
D<-N/area
}
",file = "SCRb.txt")


data<-list(y=Yaug, M=M, K=K, C=C, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
parameters<-c('psi','alpha0','alpha2','N', 'D', 'sigma')

inits =  function() {list(z=c(rep(1,nind), rbinom(nz,1,0.5)), psi=runif(1), s=Sin,
		 sigma=runif(1,2,3),alpha0=runif(1), alpha2=runif(1)) }

modelFile= "SCRb.txt"

}#end if



######################################################################################################

#run model
mod.out <- jags.model(modelFile, data, inits, n.chains=n.chains, n.adapt=n.adapt)
out <- coda.samples(mod.out,  parameters, n.iter=n.iter)



return(out)

} #end function