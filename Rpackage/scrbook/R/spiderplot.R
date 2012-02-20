spiderplot <-
function(y3d,traplocs){
# y3d = nind x ntraps x nreps
#
# compute average location of capture for each individual

plot(traplocs,pch=20,xlab=" ",ylab=" ")

avg.s<-matrix(NA,nrow=nind,ncol=2)

for(i in 1:nind){
tmp<-NULL
for(j in 1:T){
aa<-y3d[i,,j]
if(sum(aa)>0){
 aa<-  trapmat[aa>0,]
 tmp<-rbind(tmp,aa)
}
}
avg.s[i,]<-c(mean(tmp[,1]),mean(tmp[,2]))
points(avg.s[i,1],avg.s[i,2],pch=20,cex=2,col="red")
for(m in 1:nrow(tmp)){
if(nrow(tmp)>1)
lines(c(avg.s[i,1],tmp[m,1]),c(avg.s[i,2],tmp[m,2]) )
}
}
points(traplocs,pch=20)




Cx<-mean(trapmat[,1])
Cy<-mean(trapmat[,2])
avg.s<-rbind(avg.s,matrix(NA,nrow=nz,ncol=2))
 # distance from center of grid
xcent<- sqrt( (avg.s[,1]-Cx)^2 +  (avg.s[,2]-Cy)^2)

list(xcent=xcent,avg.s=avg.s,center=c(Cx,Cy))

}
