array2secr<-function(y, session, detector=c("count", "proximity")){ #y is 3D array

out<-matrix(nrow=0, ncol=4)

#sets session to 1 for all individuals
#otherwise session can be a vector of 1 and 2's coding for sex
if(missing(session))session<-rep(1,dim(y)[1])

N=dim(y)[1]
nreps=dim(y)[3]

if(detector=="proximity") {

for (i in 1:N) {
for (j in 1:nreps){
trapcap<- which(y[i,,j]>0)
ntrapcap<-length(trapcap)
if(ntrapcap > 0) {
newcap<-matrix(nrow=ntrapcap, ncol=4)
newcap[,1]<-session[i]
newcap[,2]<-i
newcap[,3]<- j
newcap[,4]<-trapcap 
out<-rbind(out, newcap)}
}}
} #end first detector 

if (detector == "count") {
for (i in 1:N) {
for (j in 1:nreps){
trapcap<- which(y[i,,j]>0)
ncap<- y[i,,j][which(y[i,,j]>0)]
capvec<-rep(trapcap, ncap)
ntrapcap<-length(capvec)
if(ntrapcap > 0) {
newcap<-matrix(nrow=ntrapcap, ncol=4)
newcap[,1]<-session[i]
newcap[,2]<-i
newcap[,3]<- j
newcap[,4]<-capvec 
out<-rbind(out, newcap)}
}}
} #end count detector

colnames(out)<-c("Session","ID","Occasion","trapID")
return(out)
}


