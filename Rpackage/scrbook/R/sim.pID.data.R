sim.pID.data <- function(N=N, K=K, sigma=sigma, lam0=lam0, knownID=knownID,X=X,
                     xlims=xlims, ylims=ylims,  obsmod= c("pois", "bern"), nmarked=c("known", "unknown"),rat=1, tel =0, nlocs=0)
{

###add an error message for when there are more tel guys than nmarked
if(tel>knownID) stop ("tel cannot be bigger than knownID")

    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)

    # Home range centers
npts<-dim(X)[1]
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    D <- e2dist(S, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    Y <- array(NA, c(N, npts, K))
for (i in 1:N){
for (j in 1: npts){

if (identical(obsmod, "bern")){
        Y[i,j,] <- rbinom(K,1, lam[i,j])
        } else if (identical(obsmod, "pois"))  {
        Y[i,j,] <- rpois(K,lam[i,j])
        }
    }}

    n <- apply(Y, c(2,3), sum)

        Yknown <- Y[1:knownID,,]

if (identical(nmarked, "unknown")){
iobs<-which(apply(Yknown>0,1,any))
Yobs<-Y[iobs,,]
} else if (identical(nmarked, "known")){
Yobs<-Yknown }

YknownR<-Yobs
counter<-array(0, c(dim(Yobs)[1],dim(X)[1],K ))
for (i in 1:dim(Yobs)[1]){
for (j in 1: dim(X)[1]){
for (k in 1:K){

if (identical(obsmod, "bern")){
if (YknownR[i,j,k] ==1 ) {
IDed<-rbinom(1,1,rat)
if (IDed ==0) { 
YknownR[i,j,k]<-0
counter[i,j,k]<-1} #counter is the number of marked records that cannot be identified to individual level
}
} else if (identical(obsmod, "pois")) {
if (Yobs[i,j,k] > 0 ) {

IDed<-sum(rbinom(Yobs[i,j,k] ,1,rat))
YknownR[i,j,k]<-IDed

if (IDed!=Yobs[i,j,k] ) { 
counter[i,j,k]<-Yobs[i,j,k]-IDed}
}
}


}}}

n<-n-apply(counter, 2:3, sum) #subtract unidentified pictures from n

#generate telemetry locations if tel>0
if (tel>0) {

itel<-sort(sample(1:knownID, tel, replace=F))
locs<-list()
for (i in 1:tel){
lx<-rnorm(nlocs, S[itel[i],1], sigma)
ly<-rnorm(nlocs, S[itel[i],2], sigma)
locs[[i]]<-cbind(lx, ly)
}

} else {
locs<-NULL
itel<-NULL}

    list(n=n,Y=Y, Yknown=Yknown, Yobs=Yobs, YknownR=YknownR, counter=sum(counter), locs=locs,telID=itel)

}



