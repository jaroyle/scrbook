wolvESA <-
function(noargs=TRUE){

library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)

####SCRsmy(y3d)



traplocs<-as.matrix(traps[,2:3])
mingridx<-min(traplocs[,1])
mingridy<-min(traplocs[,2])
traplocs[,1]<-traplocs[,1] -min(traplocs[,1])
traplocs[,2]<-traplocs[,2]- min(traplocs[,2])
traplocs<-traplocs/10000 ###units of 10 km
## set the state-space
ntraps<- nrow(traplocs)

Sgrid<-as.matrix(wolverine$grid2)
Sgrid[,1]<-Sgrid[,1]-mingridx
Sgrid[,2]<-Sgrid[,2]-mingridy
Sgrid<-Sgrid/10000 # units of 10 km

MASK<-traps[,4:ncol(traps)]
ndays<-apply(MASK,1,sum)


# posterior means of parameters for the 2km state-space
p0 <- 0.05
sigma <- 0.62


## For each state-space grid point we now compute the probability of
## encounter FOR THAT GRID POINT. This is "one minus the probability of
##     not being captured"
netp<- rep(NA,nrow(Sgrid))
for(i in 1:nrow(Sgrid)){
   d2 <- (traplocs[,1]-Sgrid[i,1])^2 + (traplocs[,2]-Sgrid[i,2])^2
  pvec <-  p0*exp(-(1/(2*sigma*sigma))*d2)
 netp[i] <- 1-prod( (1-pvec)^ndays )  # prob of being caught in any trap
}

y<-netp
x<-Sgrid
par(mar=c(3,3,3,6))
plot(x,pch=" ")
nc <- as.numeric(cut(y, 10))

   # if (col == "gray") {
    #    cc <- seq(3, 17, , 10)/20
   #     cc <- gray(cc)
    #}

cc <- topo.colors(10)
points(x, pch = 20, col = cc[nc], cex = 1)
image.scale(y, col = cc)

cat("   ",fill=TRUE)
cat("   ",fill=TRUE)

cat("Effective sample area (units of state-space grid cells): ", sum(y), fill=TRUE)
cat("Scale this by the area of the grid cell",fill=TRUE)
cat(" ", fill=TRUE)

cat("For the wolverine data using the 2km state-space grid this is : "
  , sum(y)*4, "km^2",fill=TRUE)

}
