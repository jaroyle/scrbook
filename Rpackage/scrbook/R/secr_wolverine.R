secr_wolverine <-
function(){

 data("wolverine")
 traps<-as.matrix(cbind(1:37,wolverine$wtraps))   #[,1:3]
dimnames(traps)<-list(NULL,c("trapID","x","y",paste("day",1:165,sep="")))


traps1<-as.data.frame(traps1[,1:3])
traps1$x<-as.numeric(as.character(traps1$x))
traps1$y<-as.numeric(as.character(traps1$y))

# This seems to ignore the trap operation information
trapfile1<-read.traps(data=traps1,detector="proximity")



hold<-rep(NA,nrow(traps))
for(i in 1:nrow(traps)){
hold[i]<-paste(traps[i,4:ncol(traps)],collapse="")
}
traps1<- cbind(traps[,1:3],"usage"=hold)
#
# to make use of the trap operation information you have to do this I think:
#
write.table(traps1, "traps.txt", row.names=FALSE, col.names=FALSE)
trapfile2<-read.traps("traps.txt",detector="proximity") 


wolv.dat<-wolverine$wcaps
dimnames(wolv.dat)<-list(NULL,c("Session","ID","Occasion","trapID"))
wolv.dat<-as.data.frame(wolv.dat)
wolvcapt1<-make.capthist(wolv.dat,trapfile1,fmt="trapID",noccasions=165)
wolvcapt2<-make.capthist(wolv.dat,trapfile2,fmt="trapID",noccasions=165)

gr2<-read.mask(data=gr)
gr<-(as.matrix(wolverine$grid4))
dimnames(gr)<-list(NULL,c("x","y"))
gr4<-read.mask(data=gr)
gr<-(as.matrix(wolverine$grid8))
dimnames(gr)<-list(NULL,c("x","y"))
gr8<-read.mask(data=gr)

wolv.secr0<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000)
wolv.secr2<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr2)
wolv.secr4<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr4)
wolv.secr8<-secr.fit(wolvcapt2,model=list(D~1, g0~1, sigma~1), buffer=20000,mask=gr8)
 

# reported in the book chapter:
wolv.secr2
}
