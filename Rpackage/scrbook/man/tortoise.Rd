\name{tortoise}
\alias{tortoise}
\docType{data}
\title{
desert tortoise data %%   ~~ data name/kind ... ~~
}
\description{
From Bob Seidl. Info here. %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(tortoise)}
\format{
  A data frame with 146 observations on the following 3 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Transect}}{a factor with levels \code{T107} \code{T109} \code{T111} \code{T113} \code{T117} \code{T119} \code{T120} \code{T121} \code{T123} \code{T124} \code{T128} \code{T130} \code{T134} \code{T136} \code{T142} \code{T152} \code{T153} \code{T154} \code{T155} \code{T156} \code{T157} \code{T158} \code{T160} \code{T161} \code{T162} \code{T163} \code{T165} \code{T166} \code{T167} \code{T168} \code{T170} \code{T172} \code{T173} \code{T174} \code{T179} \code{T180} \code{T182} \code{T183} \code{T188} \code{T189} \code{T190} \code{T191} \code{T193} \code{T198} \code{T202} \code{T211} \code{T213} \code{T214} \code{T215} \code{T216} \code{T218} \code{T221} \code{T222} \code{T224} \code{T231} \code{T237} \code{T238} \code{T240} \code{T247} \code{T262} \code{T301} \code{T302} \code{T304} \code{T307} \code{T312} \code{T314} \code{T316} \code{T317} \code{T322} \code{T324} \code{T325} \code{T326} \code{T327} \code{T328} \code{T331} \code{T333} \code{T335} \code{T336} \code{T337} \code{T339} \code{T340} \code{T343} \code{T346} \code{T355} \code{T359} \code{T365} \code{T366} \code{T367} \code{T368} \code{T370} \code{T371} \code{T373} \code{T375} \code{T376} \code{T383} \code{T385} \code{T390} \code{T393} \code{T397} \code{T399} \code{T407} \code{T408} \code{T409} \code{T411} \code{T412} \code{T418} \code{T419} \code{T421} \code{T422} \code{T424} \code{T431} \code{T433} \code{T447} \code{T450} \code{T455} \code{T468} \code{T469} \code{T475} \code{T477} \code{T483}}
    \item{\code{Dist}}{a numeric vector}
  }
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
data(tortoise)

data.mn<- formatDistData(tortoise,"Dist","Transect",dist.breaks=c(0:32))

library("unmarked")
 umf <- unmarkedFrameDS(y=as.matrix(data.mn), siteCovs=NULL, survey="line",
 dist.breaks=c(0:32), tlength=rep(1000, 120),
 unitsIn="m")
 summary(umf)

png("tortoise.png",width=7,height=7, units="in", res=400)
 hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)
dev.off()

m0<-distsamp(~1~1,umf,keyfun="halfnorm",output="density",unitsOut="ha")

x<-tortoise[,"Dist"]
nind<-sum(!is.na(x))
y<-rep(1,nind)
nz<-700
y<-c(y,rep(0,nz))
x<-c(x,rep(NA,nz))
z<-y

cat("
model{
alpha1~dunif(0,10)
sigma<- sqrt(1/(2*alpha1))
psi~dunif(0,1)

for(i in 1:(nind+nz)){
   z[i]~dbern(psi)    # DA Variables
   x[i]~dunif(0,B)    # B=strip width
   p[i]<-exp(logp[i])   # DETECTION MODEL
   logp[i]<-   -alpha1*(x[i]*x[i])
   mu[i]<-z[i]*p[i]
   y[i]~dbern(mu[i])  # OBSERVATION MODEL
 }

N<-sum(z[1:(nind+nz)])
D<- N/striparea  # area of transects
}
",file="dsamp.txt")
library("R2WinBUGS")
# density to be in units of ind/ha
data<-list(y=y,x=x,nz=nz,nind=nind,B=40,striparea=(120*1000*40*2/10000))
params<-list('alpha1','sigma','N','D','psi')
inits =  function() {list(z=z, psi=runif(1), alpha1=runif(1,0,.02) )}
fit = bugs(data, inits, params, model.file="dsamp.txt",working.directory=getwd(),    
       debug=T, n.chains=3, n.iter=3000, n.burnin=1000, n.thin=2)






}
\keyword{datasets}
