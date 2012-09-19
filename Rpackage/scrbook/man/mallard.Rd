\name{mallard}
\alias{mallard}
\docType{data}
\title{
mallard data  %%   ~~ data name/kind ... ~~
}
\description{
more to come %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(mallard)}
\format{
  The format is:
List of 3
 $ bandings  : num [1:405, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
 $ recoveries: num [1:405, 1:5] NA NA NA NA NA NA NA NA NA NA ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:5] "M1996" "M1997" "M1998" "M1999" ...
 $ locs      : num [1:405, 1:2] -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "long" "lat"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(scrbook)
data(mallard)
library("R2WinBUGS")


sink("model.txt")
cat("
model {
 for(t in 1:5){
    for (i in 1:nobs){
       y[i,t] ~ dbin(p[i,t], B[i,t])
       logit(p[i,t]) <- beta0[t] + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,1]*X[i,2]
     }
}
	beta1~dnorm(0,.001)
	beta2~dnorm(0,.001)
	beta3~dnorm(0,.001)
	for(t in 1:5){
 	beta0[t] ~ dnorm(0,.001)  
 }
}
",fill=TRUE)
sink()

data <- list(B=mallard$bandings, y=mallard$recoveries,
X=mallard$locs,nobs=nrow(mallard$locs))
inits <- 	function(){
list(beta0=rnorm(5),beta1=0,beta2=0,beta3=0)
}
parms <- list('beta0','beta1','beta2','beta3')
out <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=FALSE)
print(out,digits=3)

### now run to get estimates of p at each grid cell
### and make a nice image plot

parms <- list('beta0','beta1','beta2','beta3','p')
out2 <- bugs(data,inits, parms,"model.txt",n.chains=3,
 					n.iter=2000,n.burnin=1000,
					n.thin=2, debug=FALSE)

pbar<-apply(out2$sims.list$p,2,mean)
ux<-unique(mallard$locs[,1])
uy<-unique(mallard$locs[,2])

odx<-order(X[,1],X[,2])
I<-matrix(pbar[odx],nrow=length(uy),ncol=length(ux))
par(mfrow=c(2,1))
spatial.plot(X,pbar,col="ss")
image(ux,uy,t(I),col=terrain.colors(10))
spatial.plot(X,pbar,col="notgray",cx=4)




## maybe str(mallard) ; plot(mallard) ...
}
\keyword{datasets}
