### make code to run bear example in secr for Ch9, detection models; 

secr.bear<-function(){
library("secr")
library("scrbook")
data("beardata") 
traps<-as.matrix(cbind(c(1:dim(beardata$trapmat)[1]),beardata$trapmat*1000))   #[,1:3]
colnames(traps)<-c("trapID","x","y")

traps1<-as.data.frame(traps[,1:3])

trapfile1<-read.traps(data=traps1,detector="proximity")


bear.cap<-make.capthist(as.data.frame(beardata$flat), trapfile1,fmt="trapID",noccasions=8)


#check which ones to retain
bear.0=secr.fit (bear.cap, model=list(D~1, g0~1, sigma~1),buffer = 20000)
bear.0exp=secr.fit (bear.cap, model=list(D~1, g0~1, sigma~1),buffer = 20000,detectfn=2)
bear.t=secr.fit (bear.cap, model=list(D~1, g0~t, sigma~1),buffer = 20000)
bear.B=secr.fit (bear.cap, model=list(D~1, g0~b, sigma~1),buffer = 20000)
bear.b=secr.fit (bear.cap, model=list(D~1, g0~bk, sigma~1),buffer = 20000) 
bear.Bt=secr.fit (bear.cap, model=list(D~1, g0~b+t, sigma~1),buffer = 20000)
bear.sex=secr.fit (bear.cap, model=list(D~session, g0~session, sigma~session),buffer = 20000)
bear.h2=secr.fit(bear.cap, model=list(D~1, g0~h2, sigma~h2),buffer = 20000)

# produce AIC table
AIC.tab<-AIC(bear.0, bear.0exp, bear.t, bear.B, bear.b,bear.Bt ,bear.sex, bear.h2)

return(list(AIC.tab=AIC.tab,bear.0=bear.0, bear.0exp=bear.0exp, bear.t=bear.t, bear.b=bear.b, 
		bear.b=bear.b, bear.Bt=bear.Bt,bear.sex=bear.sex, bear.h2=bear.h2 ))

}#end function call





