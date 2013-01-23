library(scrbook)
data(beardata)

if(1==2){
X<-as.matrix(cbind(id=1:nrow(beardata$trapmat),beardata$trapmat))
opps<-matrix(1,nrow=nrow(X),ncol=8)
dimnames(opps)<-list(NULL,1:8)
X<-cbind(X,opps)
a<-SCR23darray(beardata$flat,X)
toad<-spiderplot(a,beardata$trapmat)
# now grab the distance from centroid variable
xcent<-toad$xcent
# see Chapter 3 of the book
}

xr<-range(beardata$trapmat[,1])
yr<-range(beardata$trapmat[,2])
xsize<-xr[2]-xr[1]
ysize<-yr[2]-yr[1]
rsize<- ysize/xsize
rsize*4

png("bear_spiderplot.png",width=4*1.5,height=2.85*1.5,units="in",res=400)
spiderplot(beardata$bearArray,beardata$trapmat)
dev.off()


