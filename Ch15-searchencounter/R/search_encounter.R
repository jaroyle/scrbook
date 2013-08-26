basex<-c(0,0,1,1,0)
basey<-c(0,1,1,0,0)

base<-cbind(basex,basey)
plot(basex,basey,pch=" ")
polygon(base)


p1<- cbind(basex,basey+3)
p2<- cbind(basex+1,basey+3)
p3<- cbind(basex+1,basey+2)
p4<-cbind(basex+2,basey+2)
p5<- cbind(basex+1,basey+1)
p6<-cbind(basex+2,basey+1)
p7<-cbind(basex+2,basey)

pp<-rbind(p1,p2,p3,p4,p5,p6,p7,base)

png("snakeline.png",width=7,height=7, units="in", res=400)
plot(pp,pch=" ")
polygon(p1)
polygon(p2)
polygon(p3)
polygon(p4)
polygon(p5)
polygon(p6)
polygon(p7)

line1<-source("line1.R")$value
if(!exists("line1"))
line1<-locator(100)

line1<-cbind(line1$x,line1$y)
lines(line1,lwd=2)
dev.off()


perbox<- 4
N<- 30*perbox
Xl<- -1
Xu<- 4
Yl<-  -1
Yu<- 5


set.seed(2013)



sx<-runif(N,Xl,Xu)
sy<-runif(N,Yl,Yu)
points(sx,sy,pch=20,col="red")

sigma<-.3
beta0<- -.5
beta1<- -1/(2*(.2^2))
X<-regpoints@coords
J<-nrow(X)

K<- 5   ##period study
U<-array(NA,dim=c(N,K,2))
y<-pmat<-matrix(NA,nrow=N,ncol=K)
for(i in 1:N){
for(k in 1:K){
U[i,k,]<-c(rnorm(1,sx[i],sigma),rnorm(1,sy[i],sigma))
dvec<-     sqrt( ( U[i,k,1] - X[,1])^2 + (U[i,k,2] - X[,2])^2  )
loghaz<- beta0 + beta1*dvec
H<- sum(exp(loghaz))
pmat[i,k]<- 1-exp(-H)
y[i,k]<- rbinom(1,1,pmat[i,k])
}
}
Ux<-U[,,1]
Uy<-U[,,2]
Ux[y==0]<-NA
Uy[y==0]<-NA
points(Ux,Uy,pch=20,col="black")

ncap<-apply(y,1,sum)
y<-y[ncap>0,]
Ux<-Ux[ncap>0,]
Uy<-Uy[ncap>0,]


##
#I,t]  ~ bvn( s[i], sigma^2)
#Log(h(u[I,t],x)) = beta0 + beta1*dist(u[I,t],x)
#Total hazard is this:
#H(u[I,t]) = exp(beta0)* sum_{j} exp(beta1*dist(u[I,t],x[j])
#x[j] = point on line
#Then: p[I,t] = 1-exp(-H(u[I,t]))



M<-100
nind<-nrow(y)
y<-rbind(y,matrix(0,nrow=(M-nrow(y)),ncol=ncol(y)))
Namat<-matrix(NA,nrow=(M-nind),ncol=ncol(y))
Ux<-rbind(Ux,Namat)
Uy<-rbind(Uy,Namat)
S<-cbind(runif(M,Xl,Xu),runif(M,Yl,Yu))
for(i in 1:nind){
S[i,]<-c( mean(Ux[i,],na.rm=TRUE),mean(Uy[i,],na.rm=TRUE))
}
Ux.st<-Ux
Uy.st<-Uy
for(i in 1:M){
Ux.st[i,!is.na(Ux[i,])]<-NA
Uy.st[i,!is.na(Uy[i,])]<-NA
Ux.st[i,is.na(Ux[i,])]<-S[i,1]
Uy.st[i,is.na(Uy[i,])]<-S[i,2]
}



### write BUGS model description

cat("
model {

# Priors
beta0~dunif(-25,25)
beta1~dunif(-25,25)
lsigma~dunif(-5,5)
sigma<-exp(lsigma)
tau<-1/(sigma*sigma)
psi~dunif(0,1)

# Likelihood
for(i in 1:M){ # Loop over individuals
 w[i]~dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu) 
 for(k in 1:K){ # Loop over temporal replicates
    u[i,k] ~ dnorm(s[i,1],tau) 
    v[i,k] ~ dnorm(s[i,2],tau) 
    for(j in 1:J){ # Loop over each point defining line segments
      d[i,k,j]<-  pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)
      h[i,k,j]<-exp(beta0+beta1*d[i,k,j])
   }
   H[i,k]<-sum(h[i,k,1:J])
   p[i,k]<- w[i]*(1-exp(-H[i,k]))
   y[i,k] ~ dbern(p[i,k])
 }
}

# Derived quantity
N<-sum(w[])
}
",file="model0.txt")


## Prepare objects to run WinBUGS from R


# Data 
data <- list (y=y, u=Ux,  v=Uy, X=X, K=K, M=M, J=J, Xl=Xl,Xu=Xu,Yl=Yl,Yu=Yu)

#Define function to generate initial values

inits <- function(){
  list("beta0"=beta0-.1,"beta1"=beta1-3,"lsigma"=log(.5),
       "s"=S,w=c(rep(1,nind),rep(0,M-nind)) ,u=Ux.st,v=Uy.st)
}

# List parameters to be estimated
parameters <- c("beta1", "N", "psi", "sigma", "beta0")

# Set MCMC settings
nthin<-1
nc<-3
nb<-1000
ni<-5000

# Load interface package and start WinBUGS (note: this takes about 4.5 h)
library("R2WinBUGS")
library("R2jags")
wbout <- jags(data, inits, parameters, "model0.txt", n.thin=nthin, n.chains=nc,
n.burnin=nb, n.iter=ni, working.dir=getwd())

# Produce summary of posterior distributions (after WinBUGS has been exited manually or change of argument in previous function to debug = FALSE)
print(wbout, dig = 3)








plot(x/100,.5+.5*sin(2*pi*x/100))








