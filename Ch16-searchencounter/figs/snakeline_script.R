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
plot(pp,xlab="easting",ylab="northing",pch=" ")
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


