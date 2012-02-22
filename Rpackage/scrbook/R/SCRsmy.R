SCRsmy <-
function(y3d){
nind<-dim(y3d)[1]
totcaps<-nperiods<-sprecaps<-rep(NA,nind)
for(i in 1:nind){
 x<- y3d[i,,]
 ntraps<-sum( apply(x,2,sum) >0)
 ncaps<- sum( x )
 nperiods[i]<- sum(apply(x,1,sum)>0)
 sprecaps[i]<- ifelse(ntraps>1,1,0)*ncaps
 totcaps[i]<-sum(x)
}
cat("Total captures: ",sum(totcaps),fill=TRUE)
cat("Spatial recaptures: ",sum(sprecaps),fill=TRUE)
cat("Ordinary capture events: ",sum(nperiods),fill=TRUE)
cat("Captures lost in non-spatial model: ",sum(totcaps)-sum(nperiods),fill=TRUE)

}
