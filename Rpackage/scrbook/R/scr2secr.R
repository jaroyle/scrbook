scr2secr <-
function(scrtraps,type="proximity"){

hold<-rep(NA,nrow(scrtraps))
for(i in 1:nrow(scrtraps)){
hold[i]<-paste(scrtraps[i,4:ncol(scrtraps)],collapse="")
}
traps1<- cbind(scrtraps[,1:3],"usage"=hold)
#
# to make use of the trap operation information you have to do this I think:
#
write.table(traps1, "traps.txt", row.names=FALSE, col.names=FALSE)
trapfile2<-read.traps("traps.txt",detector=type) 
return(trapfile2)

}
