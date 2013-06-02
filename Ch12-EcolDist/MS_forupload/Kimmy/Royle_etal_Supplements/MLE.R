########################################################################
#### Define the likelihood function.
####
	
intlik3ed<-function(start=NULL,y=y,K=NULL,X=traplocs,distmet="ecol",
    covariate,alpha2=NA){  

#start is the starting values for the parameters to be estimated
#y is the data matrix
#K is the number of occasions
#X is the traplocation matrix  with 
#   (one column for x-coords one column for y-coords)
#distmet is either "ecol" or "euclid" for least cost or Euclidean distance

#First, some input checks
if(is.null(K)) return("need sample size")
if(class(covariate)!="RasterLayer") {
 cat("make a raster out of this",fill=TRUE)
 return(NULL)
 }


# Build integration grid. This derives from the covariate raster
# i.e., potential values of s are the mid-point of each raster pixel
nc<-covariate@ncols  #number of columns in the covariate raster
nr<-covariate@nrows  #number of rows in the covariate raster
Xl<-covariate@extent@xmin  #minimum X of the covariate raster
Xu<-covariate@extent@xmax  #maximum X of the covariate raster
Yl<-covariate@extent@ymin  #minimum Y of the covariate raster
Yu<-covariate@extent@ymax  #maximum Y of the covariate raster
SSarea<- (Xu-Xl)*(Yu-Yl)   #Total area of the area over estimating density


###Create potential activity centers evenly distributed across the area
delta<- (Xu-Xl)/nc  #identify offset between locations
xg<-seq(Xl+delta/2,Xu-delta/2,delta) 
yg<-seq(Yl+delta/2,Yu-delta/2,delta) 
npix.x<-length(xg)
npix.y<-length(yg)
area<- (Xu-Xl)*(Yu-Yl)/((npix.x)*(npix.y))  # area of pixels

G<-cbind(rep(xg,npix.y),sort(rep(yg,npix.x)))
nG<-nrow(G) #total number of potential activity centers 
ntraps<- nrow(X)

if(distmet=="euclid")  #if distance metric is Euclidean,
D<- e2dist(X,G)        #calculate distance matrix among traplocs and integration grid

if(distmet=="ecol"){   #if distance metric is Ecological distance 
   if(is.na(alpha2)) alpha2<-exp(start[4]) # if estimating alpha2, use 
                                           #  this starting value
  #the next series of commands calculates the ecological distance matrix 
  cost<- exp(alpha2*covariate)  #create resistance surface
  #find neighbors and calculate conductances among neighbors
  tr1<-transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE) #adjust diag.conductances
  D<-costDistance(tr1CorrC,X,G)  #calculate the ecological distance matrix
  }

if(is.null(start)) start<-c(0,0,0,0) #if no starting values were given, assign as 0s
alpha0<-start[1]
alpha1<-start[2]
n0<-exp(start[3])

#calculate the capture probability at each trap location for each potential 
#activity center, equivalent to each location in the integration grid

probcap<- (exp(alpha0)/(1+exp(alpha0)))*exp(-alpha1*D*D) 

#create integrand matrix (# traplocations X #potential activity centers)

Pm<-matrix(NA,nrow=ntraps,ncol=nG) 
ymat<-y

#Add a zero capture record to the capture history. 
#We weight this last record later to get the total
#estimated number of individuals

ymat<-rbind(y,rep(0,ncol(y)))  
lik.marg<-rep(NA,nrow(ymat)) #create vector to store the marginal likelihood

#calculate the likelihood  
#loop through each individual (plus the all zero-capture record) in capture history

for(i in 1:nrow(ymat)){ 
#calculate the probability of capturing individual i the recorded #of times 
#in each trap (rows) when there are K occasions 
#at each potential activity center (columns).
#These are logged for numerical stability (so very small probabilities recorded)

Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))

#conditional on s likelihood (for individual i for each potential activity center)
#exponentiate to remove log above now that the numbers are larger

lik.cond<- exp(colSums(Pm,na.rm=TRUE) )

#the marginal likelihood is formed by explicitly evaluating the integrand 
#on our integration grid, averaging over the conditional likelihood 
#at all potential activity centers as described above
lik.marg[i]<- sum( lik.cond*(1/nG) )  
}      

#define weights = 1 for all individuals we captured, 
#and estimate weight for the number of individuals that were not captured =n0                 
                     
nv<-c(rep(1,length(lik.marg)-1),n0)  
#calculate combinatorial term to account for (N choose n) ways to  
#realize sample of size n
part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1) #combinatorial term  

###Multiply by weight vector(nv) and sum resulting marginal likelihoods for all individuals
part2<- sum(nv*log(lik.marg))  
out<-  -1*(part1+ part2)  #final negative log likelihood
attr(out,"SSarea")<- SSarea   #save an attribute that has just the
                              #  area of the integration grid
out   # the value of the neg log likelihood

