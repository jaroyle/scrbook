
png("bear-modelMh-post.png",width=3.5,height=3.5, units="in", res=400)

N<-source("table_of_N_fort_drum.R")$value
plot(N/sum(N),xlab="N (population size)",ylab="posterior density",
cex=1.5)

dev.off()

sp<- smooth.spline(xg[1:80],N[1:80],cv=TRUE)
    
 sp$x[sp$y==max(sp$y)]
[1] 82


sp<- smooth.spline(xg,N,cv=TRUE)
    
Call:
smooth.spline(x = xg, y = N, cv = TRUE)

Smoothing Parameter  spar= 0.09339815  lambda= 8.201724e-09 (17 iterations) 
Equivalent Degrees of Freedom (Df): 121.1825 
Penalized Criterion: 2544481 
PRESS: 5903.4 
> 
 sp$x[sp$y==max(sp$y)]
[1] 82

