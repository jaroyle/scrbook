SCRdesign <-
function(S=S,C=C,ntraps=9,ndesigns=10,nn=19,sigma=2,crit=3){

# state-space for point process


ngrid<-nrow(S)

 
Cd<-round(e2dist(C,C),8)  # distance among candidates

NN2<-NN<-matrix(0,nrow=nrow(Cd),ncol=ncol(Cd))
for(i in 1:nrow(Cd)){
xx<-Cd[i,]
NN[i,] <- (xx>0 & xx<= sort(xx)[nn])
NN2[i,] <- (xx>0 & xx<= sort(xx)[3])
}


## S as defined above is ALL locations in the state-space. Seems like the MC average
## of populations of some size N=100 say, should be the same. Lets test that.


Qfn<-function(X,S,N=100){
# computes expected information and then inverts that
beta0 <- -1.7
# Jan 2015 added the "2" here. Not in book!
beta1 <-  -1*( 1/(2*sigma^2) )
dmat<-e2dist(X,S)^2   # ntraps x nstatespace points
# Design matrix  
lammat<-exp(beta0+beta1*dmat)
lamvec<- exp(beta0 + beta1*dmat[1:length(dmat)])

lamJ<- as.vector( t(lammat)%*%rep(1,nrow(X))  )
pbar<- as.vector( 1-exp(-t(lammat)%*%rep(1,nrow(X))))
print(summary(pbar))
pbar<-mean(pbar)

M1<- rep(1,ntraps*nrow(S))
M2<- dmat[1:length(dmat)]

I11<- (1/nrow(S))*sum(lamvec)
I12<- (1/nrow(S))*sum(lamvec*M2)
I21<- (1/nrow(S))*sum(lamvec*M2)
I22<- (1/nrow(S))*sum(lamvec*M2*M2)
# this is expected information. Probs need to multiply by E[n]
I<- matrix(c(I11,I12,I21,I22),nrow=2,byrow=TRUE)
I<-  N*pbar*I 

# V is the var-cov matrix of the MLE
V<-solve(I) 
Q1<- sum(diag(V))

# matrix^2  * scalar 
sumsJ<- as.vector(  
t( lammat*lammat*(diag(V)[1] + (dmat^2)*diag(V)[2]  -2*dmat*V[1,2]       ) )%*%rep(1,nrow(X))   )
var.pbar<- ((1/nrow(S))^2)*sum(exp(-lamJ)*exp(-lamJ)*sumsJ)

#sumsJ<- as.vector(  
#t( lammat*lammat*(diag(V)[1] + (dmat^2)*diag(V)[2]       ) )%*%rep(1,nrow(X))   )
#var.pbar<- ((1/nrow(S))^2)*sum(exp(-lamJ)*exp(-lamJ)*sumsJ)

part1<- (N*N*var.pbar)  ### /(pbar*pbar)
part2<- N*(1-pbar)/pbar
total<- part1+part2

newpart2<- N*(1-pbar)*(var.pbar + 1)/pbar

## crit = 4 min var(beta-hat)
## crit = 5, min variance N-hat
## crit = 6, maximizes pbar
## crit = 7 minimizes var(pbar)
old<- N*N*var.pbar + newpart2
fixed<-  N*pbar*( (1-pbar) + N*pbar)*( var.pbar/(pbar^4) )
c(part1,newpart2,total,Q1,fixed,1-pbar,var.pbar)

}









Dlist<-list()

Qhistory<-NULL
for(m in 1:ndesigns){
Qbest<-10^10
X.current <- sample( 1:nrow(C),ntraps)
X<-C[X.current,]
Q<-Qfn(X,S)[crit]
Qhistory<-c(Qhistory,Q)
cat("Initial Q: ",Q,fill=TRUE)
if(is.nan(Q)){
Dlist[[m]]<- list(Q=NA, X=X,X.current=X.current)
 next
}


repeat{

for(i in 1:ntraps){
# have to remove X.current
chk<- NN[X.current[i],]
chk[X.current]<-0
x.consider<-(1:ncol(NN))[chk==1]
##print(x.consider)
####if( length(x.consider) ==0) next
##print(x.consider)
qtest<-rep(10^10,length(x.consider))
for(j in 1:length(x.consider)){
Xtest<-X
Xtest[i,]<- C[x.consider[j],]
xxxx<-Qfn(Xtest,S)

qtest[j]<- xxxx[crit]

}
##print(qtest)
print(qtest)
if(any(is.nan(qtest))){
Dlist[[m]]<- list(Q=NA, X=X,X.current=X.current)
 next
}
if(min(qtest)< Q){
Q<-min(qtest)
kp<- qtest==min(qtest)
X.current[i]<-x.consider[kp][1]
X<-C[X.current,]
cat("new Q: ",Q,fill=TRUE)
plot(S,pch=".")
points(X,pch=20)
}

}
cat("Current value after all J points: ", Q,fill=TRUE)
if(Qbest == Q){
break
}
if(Q<Qbest) Qbest<-Q
if(Q>Qbest) cat("ERROR",fill=TRUE)

Qhistory<-c(Qhistory,Q)
}

Dlist[[m]]<- list(Q=Qbest, X=X,X.current=X.current)
m<-m+1


## closes loop over designs
}

Qvec<-rep(NA,length(Dlist))
Xid<-matrix(NA,nrow=ntraps,ncol=length(Dlist))
Xlst<-list()
for(i in 1:length(Dlist)){
Qvec[i]<-Dlist[[i]]$Q
Xid[,i]<-Dlist[[i]]$X.current
Xlst[[i]]<-Dlist[[i]]$X
}

od<-order(Qvec)
tmp<-list()
for(i in 1:length(od)){
tmp[[i]]<-Xlst[[od[i]]]
}
Qvec<-Qvec[od]
Xid<-Xid[,od]
Xlst<-tmp


output<-list(Qvec=Qvec,Xid=Xid,Xlst=Xlst,C=C,S=S,Qhistory=Qhistory)
return(output)


}
