get.traplocs <-
function(ntraps,ssgrid=NULL){
# this function allows you to manually place some traps in the buffer
loc<-NULL
locid<-NULL
if(!is.null(ssgrid)) plot(ssgrid,pch=20)
for(i in 1:ntraps){
x<-locator(1)
d<- sqrt((x$x - pts[,1])^2 + (x$y - pts[,2])^2 )
kp<- d==min(d)

tmp<-pts[kp,]
 loc<-rbind(loc,tmp)
locid<-c(locid,(1:length(kp))[kp])
}
list(loc=loc,locid=locid)
}
