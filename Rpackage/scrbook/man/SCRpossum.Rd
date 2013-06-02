\name{SCRpossum}
\alias{SCRpossum}
\title{
JAGS analysis of possum data %%  ~~function to do ... ~~
}
\description{
uses R2jags to do an analysis of Efford's possum data........%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRpossum()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# load the data
library("secr")
data(possum)

x<-as.matrix( traps(possumCH))
# we don't use the mask in this analysis
mask<-possummask
mask[,1]<-mask[,1]-min(x[,1])
mask[,2]<-mask[,2]-min(x[,2])
x[,1]<-x[,1]-min(x[,1])
x[,2]<-x[,2]-min(x[,2])
ntraps<-nrow(x)
y<-as.matrix(possumCH)
# categorical trap variable, set to ntraps+1 for "not captured"
# secr uses 0 to indicate "not captured"
y[y==0]<-  ntraps+1




cat("
model {
psi ~ dunif(0,1)
alpha0 ~ dnorm(0,.1)
sigma ~dunif(0,1000)
alpha1<- 1/(2*sigma*sigma)

for(i in 1:M){
  z[i] ~ dbern(psi)
  S[i,1] ~ dunif(xlim[1],xlim[2])
  S[i,2] ~ dunif(ylim[1],ylim[2])
  for(j in 1:ntraps){
    #distance from capture to the center of the home range
    d[i,j] <- pow(pow(S[i,1]-X[j,1],2) + pow(S[i,2]-X[j,2],2),1)
  }
  for(k in 1:K){
    for(j in 1:ntraps){
      lp[i,k,j] <- exp(alpha0 - alpha1*d[i,j])*z[i]            
      cp[i,k,j] <- lp[i,k,j]/(1+sum(lp[i,k,]))
    }
    cp[i,k,ntraps+1] <- 1-sum(cp[i,k,1:ntraps])  # last cell = not captured
    Ycat[i,k] ~ dcat(cp[i,k,])
  }  
}   
N <- sum(z[1:M]) 
A <- ((xlim[2]-xlim[1]))*((ylim[2]-ylim[1]))
D <- N/A
Dha<- D*10000
    
}
",file="model.txt")

nind<-nrow(y)
K<-ncol(y)
M<-300
# describe the state-space 
buff<-100
xlim<-c(min(x[,1])-buff,max(x[,1])+buff)
ylim<-c(min(x[,2])-buff,max(x[,2])+buff)
# data augmentation
yaug<-matrix(ntraps+1,nrow=(M-nind),ncol=K)
Ycat<-rbind(y[1:nind,],yaug[1:(M-nind),])

# this shows a small amount of movement among trapping grids
spider<-spiderplot(y,x)
Sst<-rbind(spider$avg.s,cbind(runif(M-nind,xlim[1],xlim[2]),runif(M-nind,ylim[1],ylim[2])))

zst<-c(rep(1,nind),rep(0,M-nind))
inits <- function(){list (z=zst,sigma=runif(1,50,100) ,S=Sst) }              

parameters <- c("psi","alpha0","alpha1","sigma","N","D","Dha")
                                                                   
data <- list (X=x,K=K,Ycat=Ycat,M=M,ntraps=ntraps,ylim=ylim,xlim=xlim)         

# this takes 1-2 hours on a high-end desktop
library("R2jags")
out <- jags(data, inits, parameters, "model.txt", n.thin=1, n.chains=3, 
n.burnin=1000,n.iter=2000,DIC=FALSE)

# print a summary
out



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
