
#To run this example code, set the working directory below and place 
# the example dataset in that directory
#Also, be sure to have installed the 'rjags' library from CRAN
 

source("Royle_Converse_Data.R")



################################################################################
#BUGS CODE 

cat("
model {

#PRIORS
#abundance model
for(i in 1:M){
  group.mem[i] ~ dcat(gprobs[])
  z[i] ~ dbern(psi)
}
psi ~ dunif(0,1)
for(i in 1:n.sites){
  b.site[i] ~ dnorm(int.lam,tau.site)
}
int.lam <- 0
sigma.site ~ dunif(0,10)
tau.site<-1/(sigma.site*sigma.site)
for(i in 1:2){
  b.season[i] <- b.seas[i]
  b.seas[i] ~ dunif(-10,10)
}
b.season[3] <- -1*(b.season[1]+b.season[2])

b.fire ~ dunif(-10,10)
b.thin ~ dunif(-10,10)

#observation model
for(i in 1:M){
  cent[i,1] ~ dunif(Xl,Xu)
  cent[i,2] ~ dunif(Yl,Yu)
}

for(s in 1:24){
  bgroup.p[s] ~ dnorm(int.p,tau.p)
}

int.p ~ dunif(-10,10)
sigma.p ~ dunif(0,10)
tau.p <- 1/(sigma.p*sigma.p)

bcap.p ~ dunif(-10,10)

for(i in 1:24){
  b.dist[i] ~ dnorm(int.dist,tau.dist)
}
int.dist ~ dunif(-10,10)
sigma.dist ~ dunif(0,10)
tau.dist <- 1/(sigma.dist*sigma.dist)

#LIKELIHOOD
#abundance model
for(j in 1:n.groups){
  log(lam[j]) <- b.site[site[j]] + b.season[season[j]] + b.thin*thin[j] + b.fire*fire[j]
  gprobs[j] <- lam[j]/sum(lam[1:n.groups])
}
#observation model
for(i in 1:M){
  for(j in 1:n.traps){
    #distance from capture to the center of the home range
    d[i,j] <- pow(pow(cent[i,1]-trap.locs[j,1],2) + pow(cent[i,2]-trap.locs[j,2],2),.5)
  }
  # assumes traplocs are the SAME for all groups 
  for(k in 1:last.cap[i]){
    for(j in 1:n.traps){
      lp[i,k,j] <- (exp(bgroup.p[group.mem[i]] + bcap.p*reencounter[i,k] + 
                    b.dist[group.mem[i]]*d[i,j])*traps.avail[j,group.mem[i]])*z[i]            
      cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,]))
    }
    cp[i,k,n.traps+1] <- 1-sum(cp[i,k,1:n.traps])  # last cell = not captured
    Ycat[i,k] ~ dcat(cp[i,k,])
  }  
}   

#DERIVED PARAMETERS
#The G.N are the total population sizes by group
N.tot <- sum(z[1:M]) 
for(i in 1:M){
  group.out[i] <- group.mem[i]*z[i]
  #This will allow us to count the number of guys in each replicate
  for(j in 1:n.groups){
    g.N[j,i] <- step(0.01*(j-group.out[i])-0.02*(j-group.out[i])*(j-group.out[i])+0.001)
  }
}
for(j in 1:n.groups){
  G.N[j] <- sum(g.N[j,])
}
     
}
",file="replicated_scr_model.txt")

zst<-rep(1,length(group.mem))
gst <- group.mem
gst[is.na(group.mem)]<- sample(1:24,sum(is.na(group.mem)),replace=TRUE)
gst[!is.na(group.mem)]<-NA
inits <- function(){list(z=zst,group.mem=gst,psi=runif(1),b.site=runif(8),
sigma.site=runif(1),b.seas=runif(2),b.thin=runif(1),b.fire=runif(1),
cent=S.st,bgroup.p=runif(24),int.p=runif(1),sigma.p=runif(1),bcap.p=runif(1),
b.dist=runif(24),int.dist=runif(1),sigma.dist=runif(1))
}              

parameters <-c("psi","sigma.site","b.season","b.fire","b.thin","int.p",
"sigma.p","bcap.p","int.dist","sigma.dist","N.tot","G.N")

library("rjags")

time1<-date()
out.1 <- jags.model("replicated_scr_model.txt", data=peromyscus.data, inits=inits, n.chains=3, n.adapt=1000)
time2<-date()
out.2 <- coda.samples(out.1,parameters,n.iter=1000)
time3<-date()

