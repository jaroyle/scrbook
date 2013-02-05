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