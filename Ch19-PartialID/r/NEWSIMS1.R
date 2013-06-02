N=100
K=5
sigma=0.65
lam0=0.5
p0=7
#knownID=50

gx<-gy<-seq(-3, 3, 1)
X<-as.matrix(expand.grid(gx, gy))
delta<-2
cs<-sample(1:dim(X)[1], 5, replace=FALSE)
capsites<-X[cs,]

xlims<-ylims<-c(-5, 5)

dat<-sim.pID.data.marking(N=N, K=K, sigma=sigma, lam0=lam0, X=X,capsites=sort(cs), p0=p0,
                     xlims=xlims, ylims=ylims,  obsmod= "pois", nmarked="known",rat=1)

M=300
inits<-function(){list( S=cbind(runif(M, -5,5), runif(M, -5,5)),lam0=runif(1),p0=runif(1,6,8), sigma=runif(1), psi=runif(1)  )}

source("scrPID.marking.R")
mod<-scrPIDm(n=dat$n, X=X, y=dat$Yknown, ymark=dat$Icap, capsites=sort(cs), M=M, obsmod = "pois",nmarked="known", niters=5000, 
    xlims=xlims, ylims=ylims, inits=inits(), delta=c(0.1, 0.1, 0.1, 0.5) ) 


























sim.pID.data.marking <- function(N=N, K=K, sigma=sigma, lam0=lam0, X=X,capsites, p0=p0,
                     xlims=xlims, ylims=ylims,  obsmod= c("pois", "bern"), nmarked=c("known", "unknown"),rat=1, tel =0, nlocs=0)
{

###add an error message for when there are more tel guys than nmarked
#if(tel>knownID) stop ("tel cannot be bigger than knownID")

    obsmod <- match.arg(obsmod)
    nmarked <- match.arg(nmarked)

    # Home range centers
npts<-dim(X)[1]
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    D <- e2dist(S, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))

####capture model

pp<-exp(p0 + (-(D[,capsites]*D[,capsites]))/(2*sigma*sigma)) 
ptot<-rowSums(pp)
pp<-pp / (1 + ptot)
pp<-cbind(pp, 1-rowSums(pp[, 1:length(capsites)]))

Icap<-matrix(NA, nrow=N, ncol=length(capsites)+1)
for (i in 1:N){
Icap[i,]<-rmultinom(1,1, pp[i,])
}

#who's marked and how many?
marked<-Icap[,6]==0
nm<-sum(marked)

##resighting data
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

        Yknown <- Y[marked,,]

if (identical(nmarked, "unknown")){
iobs<-which(apply(Yknown>0,1,any))
Yobs<-Yknown[iobs,,]
} else if (identical(nmarked, "known")){
Yobs<-Yknown }

##part for imperfect individual ID of marks
YknownR<-Yobs
counter<-0
for (i in 1:dim(Yobs)[1]){
for (j in 1: dim(X)[1]){
for (k in 1:K){

if (identical(obsmod, "bern")){
if (YknownR[i,j,k] ==1 ) {
IDed<-rbinom(1,1,rat)
if (IDed ==0) { 
YknownR[i,j,k]<-0
counter<-counter+1} #counter is the number of marked records that cannot be identified to individual level
}
} else if (identical(obsmod, "pois")) {
if (Yobs[i,j,k] > 0 ) {

IDed<-sum(rbinom(Yobs[i,j,k] ,1,rat))
YknownR[i,j,k]<-IDed

if (IDed!=Yobs[i,j,k] ) { 
counter<-counter+(Yobs[i,j,k]-IDed)}
}
}


}}}


#generate telemetry locations if tel>0
if (tel>0) {

itel<-sort(sample(which(marked), tel, replace=F))
locs<-list()
for (i in 1:tel){
lx<-rnorm(nlocs, S[itel[i],1], sigma)
ly<-rnorm(nlocs, S[itel[i],2], sigma)
locs[[i]]<-cbind(lx, ly)
}

} else {
locs<-NULL
itel<-NULL}

    list(n=n,Y=Y, Yknown=Yknown, Yobs=Yobs, YknownR=YknownR, counter=counter, locs=locs,telID=itel, nm=nm, Icap=Icap[Icap[,6]==0,])

}



