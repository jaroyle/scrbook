SCRpossum <-
function(){

######
###### possum analysis
######

# load the data
library("secr")
data(possum)
# do some standard plots and summaries using secr functions
png("possum.png",width=7,height=7, units="in", res=400)
plot(possummask)
plot(possumCH, tracks = TRUE, add = TRUE)
plot(traps(possumCH), add = TRUE)
lines(possumarea)
dev.off()
summary(possumCH)

#####

x<-as.matrix( traps(possumCH))
mask<-possummask
mask[,1]<-mask[,1]-min(x[,1])
mask[,2]<-mask[,2]-min(x[,2])
x[,1]<-x[,1]-min(x[,1])
x[,2]<-x[,2]-min(x[,2])
ntraps<-nrow(x)
y<-as.matrix(possumCH)
# categorical trap variable, set to ntraps+1 for "not captured"
# secr uses 0 to indicate "not captured"
y[y==0]<-  ntraps+1
png("possum.png",width=8,height=4.5, units="in", res=400)
# this shows a small amount of movement among trapping grids
spider<-spiderplot(y,x)
dev.off()


cat("
model {
psi ~ dunif(0,1)
alpha0 ~ dnorm(0,.1)
sigma ~dunif(0,1000)
alpha1<- 1/(2*sigma*sigma)

for(i in 1:M){
  z[i] ~ dbern(psi)
  S[i,1] ~ dunif(xlim[1],xlim[2])
  S[i,2] ~ dunif(ylim[1],ylim[2])
  for(j in 1:ntraps){
    #distance from capture to the center of the home range
    d[i,j] <- pow(pow(S[i,1]-X[j,1],2) + pow(S[i,2]-X[j,2],2),1)
  }
  for(k in 1:K){
    for(j in 1:ntraps){
      lp[i,k,j] <- exp(alpha0 - alpha1*d[i,j])*z[i]            
      cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,]))
    }
    cp[i,k,ntraps+1] <- 1-sum(cp[i,k,1:ntraps])  # last cell = not captured
    Ycat[i,k] ~ dcat(cp[i,k,])
  }  
}   
N <- sum(z[1:M]) 
A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
D <- N/A
Dha<- D*10000
    
}
",file="model.txt")


nind<-nrow(y)
K<-ncol(y)
M<-300
buff<-100
xlim<-c(min(x[,1])-buff,max(x[,1])+buff)
ylim<-c(min(x[,2])-buff,max(x[,2])+buff)
yaug<-matrix(ntraps+1,nrow=(M-nind),ncol=K)
Ycat<-rbind(y[1:nind,],yaug[1:(M-nind),])
Sst<-rbind(spider$avg.s,cbind(runif(M-nind,xlim[1],xlim[2]),runif(M-nind,ylim[1],ylim[2])))

zst<-c(rep(1,nind),rep(0,M-nind))
inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst) }              

parameters <- c("psi","alpha0","alpha1","sigma","N","D","Dha")
                                                                   
data <- list (X=x,K=K,Ycat=Ycat,M=M,ntraps=ntraps,ylim=ylim,xlim=xlim)         
library("R2jags")
out <- jags(data, inits, parameters, "model.txt", n.thin=1,n.chains=3, 
n.burnin=1000,n.iter=2000,DIC=FALSE)
  
out

}
