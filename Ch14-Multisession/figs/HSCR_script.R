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





model {
# This will show that psi and b0 
#   are confounded. 
  p~ dunif(0,1)
  b0~dnorm(0,.1)
  b1~dnorm(0,.1)
  psi ~ dunif(0,1)
  for(s in 1:S){
    log(lam[s]) <- b0 + b1*x[s]
    gprobs[s]<- lam[s]/sum(lam[1:S])
  }
  for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] ~ dbern(psi)
   y[i]~ dbin(mu[i],J)
   mu[i] <- z[i]*p
  }
  N <- sum(z[1:M]) 
}
\end{verbatim}
}
\end{minipage}
&
\begin{minipage}{2.75in}
{\small
\begin{verbatim}
model {
# This version constrains psi with 
#   the intercept parameter
  p~ dunif(0,1)
  b0~dnorm(0,.1)
  b1~dnorm(0,.1)
  psi<- sum(lam[])/M
  for(j in 1:K){
    log(lam[j]) <- b0 + b1*x[j]
    gprobs[j]<- lam[j]/sum(lam[1:K])
  }
  for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] ~ dbern(psi)
   y[i]~ dbin(mu[i],J)
   mu[i] <- z[i]*p
  }
  N <- sum(z[1:M]) 
}
\end{verbatim}
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