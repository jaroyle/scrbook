hra <-
function(func,parms,plot=TRUE,xlim,ylim){

# individual centered at (3,3) with 10000 "potential trap" locations
s<-c( (xlim[2]-xlim[1])/2, (ylim[2]-ylim[1])/2 )
x1<- rep(seq(xlim[1],xlim[2],,100),100)
x2<- sort(rep(seq(ylim[1],ylim[2],,100),100))
delta<- min(diff(x1[1:10]))
#area<- ( (6+delta/2)^2)/10000
# plot those traps
X<-cbind(x1,x2)
#plot(X,pch=".")
# compute distances and encounter probabilities
# under some model

D<-sqrt(  ( s[1]-x1)^2  + (s[2]-x2)^2)
#sigma<- .396
#lp<-  -2.0 -(1/(2*sigma*sigma))*D*D
#p<- 1-exp(-exp(lp))
# 2 alternative link functions
#p<-  exp(lp)
#p<- expit(lp)
p<- func(parms,D)

# plot this surface
if(plot){
spatial.plot(X,p)
}
# Imagine that capture probability is perfect and we can do a single sample
# Where an individual is captured in that case is the outcome of a single movement
# of the animal. Therefore where an individual is captured should occur in proportion to 
# the value of "p" like this:

psi<- p/sum(p)

# want 95% of mass -- smallest area that contains 95% of the mass.  If symmetric, should be
# symmetric integration problem.
# sum cells < k, and then vary k

# find area that contains 95% of movements
x0<-.2
repeat{
in.hr<- D<=x0
total<- sum(psi[in.hr])  
cat("Total probability: ",total,fill=TRUE)
if(total>=.95) {
print(x0)    # if condition is met, break
break
}
x0<-x0*1.01  # otherwise increase x0
}
radius<-x0

cat("radius to achieve 95% of area: ", radius, fill=TRUE)
area<- pi*radius^2
cat("home range area: ", area,fill=TRUE)
return(area)
}
