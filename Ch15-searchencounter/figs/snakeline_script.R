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
plot(pp,xlab="Easting",ylab="Northing",pch=" ",cex.axis=1.5,cex=2,cex.lab=1.5)
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

### Part 2:  Use SpatialPoints and other functions to chop the lines up into 
### a regular mesh of points
library(rgeos)
library(sp)
line1<-source("line1.R")

x<-line1$value$x
y<-line1$value$y

line1<-cbind(x,y)

line1<- as.matrix(line1)
points<-SpatialPoints(line1)

sLine<-Line(points)
###sLine<-SpatialLines(sLine)

## should be 250 or higher
regpoints<-sample.Line(sLine,100,type="regular")

plot(line1,type="l")
#original points:
points(points,col="grey")

#regularly spaced points
points(regpoints,col="red",pch=20,lwd=2)



### Part 3


perbox<- 4
N<- 30*perbox
xlim<-c(-1,4)
ylim<-c(-1,5)


set.seed(2014)



sx<-runif(N,xlim[1],xlim[2])
sy<-runif(N,ylim[1],ylim[2])
points(sx,sy,pch=20,col="red")

sigma.move<- .35
sigma<-.4
alpha0<- .8
alpha1<- 1/(2*(sigma^2))
X<-regpoints@coords
J<-nrow(X)

K<- 10   ##period study
U<-array(NA,dim=c(N,K,2))
y<-pmat<-matrix(NA,nrow=N,ncol=K)
for(i in 1:N){
for(k in 1:K){
U[i,k,]<-c(rnorm(1,sx[i],sigma.move),rnorm(1,sy[i],sigma.move))
dvec<-     sqrt( ( U[i,k,1] - X[,1])^2 + (U[i,k,2] - X[,2])^2  )
loghaz<- alpha0 - alpha1*dvec*dvec
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



M<-200
nind<-nrow(y)
y<-rbind(y,matrix(0,nrow=(M-nrow(y)),ncol=ncol(y)))
Namat<-matrix(NA,nrow=(M-nind),ncol=ncol(y))
Ux<-rbind(Ux,Namat)
Uy<-rbind(Uy,Namat)
S<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
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
alpha0~dunif(-25,25)
alpha1~dunif(0,60)
sigma<- sqrt(1/(2*alpha1))
lsigma~dunif(-5,5)
sigma.move<-exp(lsigma)
tau<-1/(sigma.move*sigma.move)
psi~dunif(0,1)

# Likelihood
for(i in 1:M){ # Loop over individuals
 z[i]~dbern(psi)
 s[i,1]~dunif(xlim[1],xlim[2])
 s[i,2]~dunif(ylim[1],ylim[2])
 for(k in 1:K){ # Loop over temporal replicates
    u[i,k] ~ dnorm(s[i,1],tau) 
    v[i,k] ~ dnorm(s[i,2],tau) 
    for(j in 1:J){ # Loop over each point defining line segments
      d[i,k,j]<-  pow(pow(u[i,k]-X[j,1],2) + pow(v[i,k]-X[j,2],2),0.5)
      h[i,k,j]<-exp(alpha0-alpha1*d[i,k,j]*d[i,k,j])
   }
   H[i,k]<-sum(h[i,k,1:J])
   p[i,k]<- z[i]*(1-exp(-H[i,k]))
   y[i,k] ~ dbern(p[i,k])
 }
}

# Derived quantity
N<-sum(z[])
}
",file="model0.txt")

## Prepare objects to run WinBUGS from R
# Data 
data <- list(y=y, u=Ux,  v=Uy, X=X, K=K, M=M, J=J,xlim=xlim,ylim=ylim)

#Define function to generate initial values
inits <- function(){
  list(alpha0=alpha0-.3,alpha1=alpha1-1.5,lsigma=log(.5),
       s=S,z=c(rep(1,nind),rep(0,M-nind)) ,u=Ux.st,v=Uy.st)
}

# List parameters to be estimated
parameters <- c("alpha0","alpha1", "N", "psi", "sigma.move","sigma")

# Set MCMC settings
nthin<-1
nc<-3
nb<-500
ni<-3500

# Load interface package and start WinBUGS (note: this takes about 4.5 h)
library("R2jags")
wbout2 <- jags(data, inits, parameters, "model0.txt", n.thin=nthin, n.chains=nc,
n.burnin=nb, n.iter=ni, working.dir=getwd())








plot(x/100,.5+.5*sin(2*pi*x/100))








