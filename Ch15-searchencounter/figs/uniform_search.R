uniform_search<-function(){

# 100 individuals 4 periods (would be years in a CJS model but phi=1 here)
N<-100
nyear<- 4
Sx<-Sy<-matrix(NA,nrow=N,ncol=nyear)
sigma.move<- .25

# simulate initial coordinates on the square:
Sx[,1]<-runif(N,0,16)
Sy[,1]<-runif(N,0,16)

# increment time and simulate ANGLE and RADIAL DISTANCE of movement:
for(t in 2:nyear){
#phi<-runif(nind,0,360)
#r<- rexp(nind,1)
#dx<- r*cos(phi)
#dy<- r*sin(phi)
## SEE here: NEXT LOCATION IS DERIVED
#Sx[,t]<-Sx[,t-1] + dx
#Sy[,t]<-Sy[,t-1]+dy
Sx[,t]<- rnorm(N,Sx[,t-1],sigma.move)
Sy[,t]<- rnorm(N,Sy[,t-1],sigma.move)
}

# now we generate encounter histories:
f<-rep(NA,N)
Y<-matrix(0,nrow=N,ncol=nyear)
for(i in 1:N){
for(t in 1:nyear){
  # see here: IF individual is IN THE SQUARE we can capture it:
 if( Sx[i,t] > 3 & Sx[i,t]< 13 & Sy[i,t]>3 & Sy[i,t]<13 )
 Y[i,t]<-rbinom(1,1,.5)
}
 # f[i] = 1 period of 1st capture, NA if never captured
 if(sum(Y[i,])>0)
  f[i]<-  min((1:nyear)[Y[i,]==1][1])

}

# subset data. If a guy is never captured, cannot have him in our data set
Y<-Y[!is.na(f),]
Sx<-Sx[!is.na(f),]
Sy<-Sy[!is.na(f),]

Sx[Y==0]<-NA
Sy[Y==0]<-NA

## data augmentation
M<-200
Y<-rbind(Y,matrix(0,nrow=(M-nrow(Y)),ncol=nyear))
Sx<- rbind(Sx,matrix(NA,nrow=(M-nrow(Sx)),ncol=nyear))
Sy<- rbind(Sy,matrix(NA,nrow=(M-nrow(Sy)),ncol=nyear))

# make 3-d array of coordinates
G<-array(NA,dim=c(M,nyear,2))
G[,,1]<-Sx
G[,,2]<-Sy


# BUGS model
sink("model1.mod")
cat("
model {
psi ~ dunif(0,1)
tau ~dgamma(.1,.1)
p0 ~ dunif(0,1)
sigma.move<- sqrt(1/tau)


# Likelihood 
for (i in 1:M){
z[i] ~ dbern(psi)

  G[i,1,1] ~ dunif(0,16)
  G[i,1,2] ~ dunif(0,16)

   for (t in 2:n.occasions){
   
## See here I can only make a model for LOCATION
      G[i,t,1] ~ dnorm(G[i,t-1,1], tau)
      G[i,t,2] ~ dnorm(G[i,t-1,2], tau)
      
# Test whether the actual location is in- or outside the study area. Needs to be done for each grid cell
    }
   for(t in 1:n.occasions){ 
      inside[i,t] <- step(G[i,t,1]-3) * step(13-G[i,t,1]) *step(G[i,t,2]-3) * step(13-G[i,t,2])
      Y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p0 * inside[i,t] * z[i]
      } #t
   } #i
N<- sum(z[])
}
",fill = TRUE)
sink()



# MCMC settings
ni <- 500
nt <- 1
nb <- 300
nc <- 2

# Bundle data
jags.data <- list(Y = Y,  G = G, M=M, n.occasions = ncol(Y))
inits <- function(){list(psi = runif(1, 0, 1),p0 = runif(1, 0, 1), 
z=c(rep(1,N),rep(0,M-N)) , tau = 1) }

library("R2jags")

# Parameters monitored
parameters <- c("p0","tau","N","psi","sigma.move")

# Call JAGS from R
model1 <- jags(jags.data, inits, parameters, "model1.mod", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


return(model1)

}













