SCR23darray <-
function(caps,traps,ntraps=NULL,nperiods=NULL){
nind<-max(caps[,2])
ntraps<-nrow(traps)
nperiods<-ncol(traps)-3
per.id<- as.numeric(dimnames(traps)[[2]][4:ncol(traps)])

if( length(per.id) != length(min(per.id):max(per.id)) ){
 x<- 1:nperiods
 names(x)<-as.character(per.id)
 cap.period<- x[as.character(caps[,3])]
}
else{
cap.period<-caps[,3]
}

y<-array(0,c(nind,nperiods,ntraps))

tmp<-cbind(caps[,2],cap.period,caps[,4])
y[tmp]<-1
y
}
