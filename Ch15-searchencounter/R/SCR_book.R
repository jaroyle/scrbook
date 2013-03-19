

if(1==2){
data<-read.csv("effort2000.csv")
mask<-data[,8]
mask<-ifelse(mask=="y",1,0)
mask[apply(Ygen,2,sum)>0]<-1
effort<-as.numeric(data[,9])

grid<-data[,2:3]

grid<-as.matrix(grid)
grid<-grid/5000
grid[,1]<-grid[,1]-min(grid[,1])
grid[,2]<-grid[,2]-min(grid[,2])
Dmat<-as.matrix(dist(grid))
plot(grid)

Ygen<- source("Ygen2000.R")$value
}
fisher<-list(
effort=effort, Ygen=Ygen, mask=mask,grid=grid)
save("fisher",file="fisher.Rda")


load("fisher.Rda")
###source("fisher.R")

ni<- 6000
nb<-1000
nthin<-1
nc<-3
library("R2WinBUGS")

effort<-fisher$effort
Ygen<-fisher$Ygen
mask<-fisher$mask
grid<-fisher$grid

nz<- 500
delta<-0  # grid is already buffered
Xl<-min(grid[,1] - delta)
Xu<-max(grid[,1] + delta)
Yl<-min(grid[,2] - delta)
Yu<-max(grid[,2] + delta)



effort<-effort[mask==1]
effort[effort==0]<-1
grid.Ygen<-grid[mask==1,]

Ygen<-Ygen[,mask==1]
ngrid.Ygen<-nrow(grid.Ygen)
ngen<-dim(Ygen)[1]
Ygen<-rbind(Ygen,matrix(0,nrow=nz,ncol=ncol(Ygen)))
ngen<-ngen+nz

sink("modelfile.txt")
cat("
model {
sigma~dunif(0,10)
psi ~ dunif(0,1)
lam0~dgamma(.1,.1)
alpha~dunif(0,1)
#p0~dunif(0,1)
for(i in 1:ngen){
 w[i]~dbern(psi)
 x0g[i]~dunif(Xl,Xu)
 y0g[i]~dunif(Yl,Yu) 
for(j in 1:ngrid.Ygen){
dist2[i,j]<- ( pow(x0g[i]-grid.Ygen[j,1],2)  +  pow(y0g[i]-grid.Ygen[j,2],2) )
#tmp[i,j]<-  p0*w[i]*exp(-dist2[i,j]/sigma)
lam[i,j]<- w[i]*pow(effort[j],alpha)*lam0*exp(-dist2[i,j]/sigma)
tmp[i,j]<-  1-exp(-lam[i,j])


Ygen[i,j] ~ dbin(tmp[i,j],1)
}
}
N<-sum(w[1:ngen])
}
",fill=TRUE)
sink()

data <- list ("Ygen","grid.Ygen","ngen","ngrid.Ygen","Yl","Yu","Xl","Xu","effort")
wst<-c(rep(1,20),rep(0,ngen-20))
inits <- function(){
  list (sigma=runif(1,.4,1),w=wst,psi=runif(1),alpha=runif(1,0,1),lam0=runif(1,0,1))
}
parameters <- c("psi","sigma","N","lam0","alpha") 




out <- bugs (data, inits, parameters, "modelfile.txt", n.thin=nthin,n.chains=nc, n.burnin=nb,n.iter=ni,debug=FALSE,working.directory=getwd())








out



}

