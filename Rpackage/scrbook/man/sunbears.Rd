\name{nybears}
\alias{nybears}
\docType{data}
\title{
Black bear data from Cat Sun and Angela Fuller, USGS and Cornell University%%   ~~ data name/kind ... ~~
}
\description{
Encounter history data from 33 individuals encountered by an array of
103 hair snares. Includes locations of hair snares ("traplocs"), a
covariate "elevation", and additional auxiliary information. 
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(nybears)}
\format{
  The format is:
List of 6
 $ y2d        : num [1:50, 1:49] 0 0 0 0 0 0 0 0 0 0 ...
 $ K          : num 10
 $ traplocs   : num [1:49, 1:2] 5 5 5 5 5 5 5 10 10 10 ...
 $ elevation  : num [1:1600, 1] -1.207 -0.829 -0.056 -1.393 -0.894 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:1600] "1" "2" "3" "4" ...
  .. ..$ : NULL
 $ trap2raster: int [1:49] 165 365 565 765 965 1165 1365 170 370 570 ...
 $ ssgrid     : int [1:1600, 1:2] 1 2 3 4 5 6 7 8 9 10 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "Var1" "Var2"
}
\details{
  

  %%  ~~ If necessary, more details than the __description__ above ~~
}
\source{



  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Royle, J. A., Chandler, R. B., Sun, C. C., & Fuller,
A. K. (2013). Integrating resource selection information with spatial
capture-recapture. Methods in Ecology and Evolution, 4:520-530.

Sun, C.C. (2013). Population Estimation, Genetic Diversity, and
Structure of Black bears in Southwestern New York, Master's Thesis. Cornell University.
  %%  ~~ possibly secondary sources and usages ~~
}
\examples{




## This basically recreates all of the analyses in Royle, Chandler, Sun
##   and Fuller, 2013, Methods in Ecology


library(scrbook)
data(sunbears)
y<-sunbears$y2d
K<-sunbears$K
X<-sunbears$traplocs
z<-sunbears$elevation
t2r<-sunbears$trap2raster
gr<-sunbears$ssgrid
ntel<-sunbears$ntel

###sunbears<-list(y2d=y,K=5,traplocs=X,elevation=as.vector(z2),trap2raster=trap2raster,ssgrid=gr2,ntel=ntel2)



# Basic SCR model with RSF covariate at trap locations.  The covariate is passed 
# to the likelihood as "ztrap"
tmp1<-nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=z[t2r],G=gr,hessian=TRUE)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp2<-nlm(intlik3rsf,c(-4.7,log(1),-.3,5),y=y,K=K,X=X,
ztrap=z[t2r],G=gr,ntel=ntel,zall=as.vector(z),hessian=TRUE)

# Fits SCR model with isotropic Gaussian encounter model
tmp3<- nlm(intlik3rsf,c(-3,log(3),1,0),y=y,K=K,X=X,ztrap=rep(0,length(t2r)),G=gr,hessian=TRUE)

# Fits SCR model with isotropic Gaussian encounter model
tmp4<- nlm(intlik3rsf,c(-3,log(3),.17,4,0),y=y,K=K,X=X,ztrap=z[t2r],G=gr,zall=as.vector(z),hessian=TRUE)

# use telemetry data and activity centers for those are marginalized out of the likelihood
tmp5<-nlm(intlik3rsf,c(-4.7,-.8,-.3,5,0),y=y,K=K,X=X,ztrap=z[t2r],G=gr,ntel=ntel,zall=z,hessian=TRUE)

tmp6<- nlm(intlik3rsf,c(-3.8,-1.2,0,5.4,0),y=y,K=K,X=X,ztrap=rep(0,length(t2r)),G=gr,zall=z,hessian=TRUE)



# Estimates of the 6 models (and SEs)
#
#                 alpha0   log(sigma)    alpha2      log(n0)  beta        -loglik
#SCR+p(x)      -2.8561676 -1.1174638  0.1747187   4.1395461                122.738
#   SE          0.3899063  0.1389833  0.2477921   0.3656961
#SCR           -2.729194   -1.122389   1.000000   4.109886                 122.990   
#   SE          0.3453705   0.1403783             0.3618065
#SCR+D(x)      -2.715320  -1.133076   0.000000   4.113903   1.247256       118.007
#   SE          0.3526155 0.1394352              0.3575286 0.4083330       
#SCR+p(x)+D(x) -2.4838347 -1.1567458 -0.3842881  4.2547317  1.5710664      117.075
#               0.3910420  0.1421062  0.2760693  0.3767896  0.4630096
#SCR+RSF       -3.0676938 -0.8141204 -0.2810946   3.8841581               1271.739
#   SE          0.27218129  0.03640307 0.11759909 0.36255170
#SCR+RSF+D(x)  -3.0701403 -0.8100523 -0.3706265  4.0284284  1.272629      1266.700
#   SE          0.27199799 0.03683849 0.12387969 0.36606116 0.411030



use<- exp(-.371*z)
odds<- use/exp(-.371*mean(z) )
# odds ratio of use of x relative to average pixel  at same distance from s




area<- 4375.36  ## km^2
area.per.pixel<- area/nrow(gr)

Nhat<-exp(4.028)+nrow(y)

Dhat<- exp(1.27*z)
Dhat<-(Dhat/sum(Dhat))*Nhat   # individuals per pixel
Dhat<- (Dhat/area.per.pixel)*100

#png("spaceusage.png",width=7,height=7,units="in",res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr,odds)
#dev.off()

#png("density.png",width=7,height=7, units="in", res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr,Dhat)
#dev.off()

#png("elev_captures.png",width=7,height=7,units="in",res=400)
par(mar=c(3,3,3,6))
spatial.plot(gr,z)
tmp<-X[col(y)[y>0],]
tmp<-tmp + rnorm(prod(dim(tmp)),0,.1)
points(tmp,pch=20)  # traps where captures happened.
#dev.off()





## maybe str(sunbears) ; plot(sunbears)
}
\keyword{datasets}
