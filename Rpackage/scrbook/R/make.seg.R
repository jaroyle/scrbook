make.seg <-
function(npts){
cat("Please click on the map ", npts, " times", fill=TRUE)
l2<-locator(npts)
l2<-cbind(l2$x,l2$y)
l2<-round(l2,2)
tmp<-NULL
for(i in 1:nrow(l2)){
tmp<-paste(tmp,l2[i,1],l2[i,2])
if(i<nrow(l2))
tmp<-paste(tmp,",")
}
l2b<- paste("LINESTRING(",tmp,")")
l2<- readWKT(l2b)
return(l2)
}
