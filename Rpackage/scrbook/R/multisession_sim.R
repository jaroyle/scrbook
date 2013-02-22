multisession_sim <-
function(){

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
# This version absorbs psi into the multinomial probabilities
# it produces a different answer because the unif(0,1) prior on psi
# is highly informative in this case (see Ch 10 of Royle and Dorazio 2008
# for an explanation)
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
out1a = jags(data1, inits, params1, model.file = "model1a.txt", 
            working.directory = getwd(), n.chains = 3, n.iter = 2000, n.burnin = 1000, n.thin = 1)
out1b = jags(data1, inits, params1, model.file = "model1b.txt", 
            working.directory = getwd(), n.chains = 3, n.iter = 2000, n.burnin = 1000, n.thin = 1)

#inits = function() {
#list(p = runif(1), beta0=rnorm(1),beta1=rnorm(1) )  }
#out2 = jags(data1, inits, params1, model.file = "model2.txt", 
#            working.directory = getwd(), n.chains = 3, n.iter = 12000, n.burnin = 1000, n.thin = 1)
return(list(out1a,out1b))
}
