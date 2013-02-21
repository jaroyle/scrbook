multisession_sim<-function(){

set.seed(2013)

grid1<-expand.grid(1:5,1:5)/10

#grid2<-grid1; grid2[,2]<-grid2[,2]+1
#grid3<- grid1; grid3[,2]<-grid3[,2]+2
#grid4<-grid1; grid4[,1]<-grid4[,1]+1
#grid5<-grid1; grid5[,1]<-grid5[,1]+1; grid5[,2]<-grid5[,2]+1
#grid6<-grid1; grid6[,1]<- grid6[,1]+1; grid6[,2]<-grid6[,2]+2
#plot(grids<-rbind(grid1,grid2,grid3,grid4,grid5,grid6))
#Sx<-runif(875,-.5,2.0)
#Sy<-runif(875,-.5,3.0)
#S<-cbind(Sx,Sy)
#plot(S,pch=20)
#points(grids,pch="x")

G<- 20
beta0<- 3
beta1<- .6
p<- .3
K<- 5 #sample occasions for capture-recapture
x<- rnorm(G)
lambda<- exp(beta0+beta1*x)
N<-rpois(G,lambda=lambda)
y<-NULL
for(g in 1:G){
if(N[g]>0)
 y<-c(y,  rbinom(N[g],K,p))
}
g<- rep(1:G,N)

g<-g[y>0]
y<-y[y>0]

M<- 1200
g<-c(g,rep(NA,M-length(g)))
y<-c(y,rep(0,M-length(y)))

cat("
model {
# This will show that psi and b0 
#   are confounded. 
  p~ dunif(0,1)
  beta0~dnorm(0,.1)
  beta1~dnorm(0,.1)
  psi ~ dunif(0,1)
  for(j in 1:G){
    log(lam[j]) <- beta0 + beta1*x[j]
    gprobs[j]<- lam[j]/sum(lam[1:G])
  }
  for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] ~ dbern(psi)
   y[i]~ dbin(mu[i],K)
   mu[i] <- z[i]*p
  }
  N <- sum(z[1:M]) 
}
",file="model1a.txt")
cat("
model {
# This version constrains psi with 
#   the intercept parameter
  p~ dunif(0,1)
  beta0~dnorm(0,.1)
  beta1~dnorm(0,.1)
  psi<- sum(lam[])/M
  for(j in 1:G){
    log(lam[j]) <- beta0 + beta1*x[j]
    gprobs[j]<- lam[j]/sum(lam[1:G])
  }
  for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] ~ dbern(psi)
   y[i]~ dbin(mu[i],K)
   mu[i] <- z[i]*p
  }
  N <- sum(z[1:M]) 
}
",file="model1b.txt")
cat("
model {
# This version constrains psi with 
#   the intercept parameter
  p ~ dunif(0,1)
  beta0 ~ dnorm(0,.1)
  beta1 ~ dnorm(0,.1)
  psi<- sum(lam[])/M
  for(j in 1:G){
    log(lam[j]) <- beta0+ beta1*x[j]
    gprobs[j]<- psi*lam[j]/sum(lam[1:G])
  }
  gprobs[G+1]<- (1-psi)

  for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] <- 1 - (g[i] == (G+1))
    y[i] ~ dbin(mu[i],K)
    mu[i] <- z[i]*p
  }
  N <- sum(z[1:M]) 
}
",file="model2.txt")
data1 <- list(y = y, g=g,M=M,K=K,G=G,x=x)
params1 = c("p", "beta0","beta1","psi","N")
inits = function() {
list(z = as.numeric(y >= 1), p = runif(1), beta0=rnorm(1),beta1=rnorm(1) )  }
library("R2jags")
out1 = jags(data1, inits, params1, model.file = "model1b.txt", 
            working.directory = getwd(), n.chains = 3, n.iter = 12000, n.burnin = 1000, n.thin = 1)

inits = function() {
list(p = runif(1), beta0=rnorm(1),beta1=rnorm(1) )  }

out2 = jags(data1, inits, params1, model.file = "model2.txt", 
            working.directory = getwd(), n.chains = 3, n.iter = 12000, n.burnin = 1000, n.thin = 1)
return(out)
}








set.seed(2013)
ngroups<- 20
# km of water boundary
Xkm<- 3+exp(rnorm(ngroups,2,1))
# area of watershed
Xarea<- runif(ngroups,19,300)
# buildings in watershed
Xbld<- rpois(ngroups, 12*mean(Xarea))
Dbld<- Xbld/Xarea
Dbld<- (Dbld - mean(Dbld))/sqrt(var(Dbld))
beta<- -.2
#lambda<- exp(1 + beta*Dbld)*Xkm

lambda<- exp(log(Xarea) + log(Xkm)*.85)
lam.tot<- sum(lambda)   # this is a huge number


# could do this by putting density directly in lambda
# or else you can renormlized based on a desired value of $\psi$ as follows

Ntotal<-  sum(Xarea)*.1   # desired density
lambda.renorm<- Ntotal*(lambda/lam.tot)   # forces lambdas to sum to Ntotal


cellprobs<-  lambda/sum(lambda)

psi<- .5
M<- 600  # total data augmentation size -- check this to be 
           # larger than sum(lambda)


newcellprobs<- c(cellprobs*(1-psi),psi)



psi<- (1/M)*lam.tot   # data augmentation parameter chosen so that
                      # E[Ntot] = sum(lam.tot)

g<- sample(1:ngroups,size=M,replace=TRUE,prob=cellprobs)
z<- rbinom(M,size=1,prob=psi)