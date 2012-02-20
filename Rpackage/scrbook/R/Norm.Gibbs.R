Norm.Gibbs <-
function(y=y,mu_0=mu_0,sigma2_0=sigma2_0,a=a,b=b,niter=niter){

ybar<-mean(y)
n<-length(y)
mu<-1         #mean initial value
sigma2<-1     #sigma2 initial value
an<-n/2 + a   #shape parameter of IvGamma distribution of sigma2

out<-matrix(nrow=niter, ncol=2)
colnames(out)<-c('mu', 'sig')

for (i in 1:niter) {

#update mu according to Eq. 7.2
mu_n<- (sigma2/(sigma2+n*sigma2_0))*mu_0 + (n*sigma2_0/(sigma2 + n*sigma2_0))*ybar 
sigma2_n <- (sigma2*sigma2_0)/ (sigma2 + n*sigma2_0)
mu<-rnorm(1,mu_n, sqrt(sigma2_n))

#update sigma2 according to Eq. 7.3
bn<- 0.5 * (sum((y-mu)^2)) + b
sigma2<-1/rgamma(1,shape=an, rate=bn)
out[i,]<-c(mu,sqrt(sigma2))

}
return(out)
}

