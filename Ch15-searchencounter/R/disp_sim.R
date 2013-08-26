# Here is a simple closed population capture-recapture model that allows individuals to move
# about according to some model for RADIAL DISTANCE. They move about a square [3,13] x [3,13]
# and they can be captured IF AND ONLY IF they are located within the square.
# I can simulate the data with a model on RADIAL DISTANCE but I can only fit the
# model if I put a model on LOCATION
# My question to everyone: Do you see a way to model RADIAL DISTANCE within JAGS/BUGS?
# The technical issue is that we have DATA  = LOCATION but we need to specifiy Location as a DERIVED
# variable (except for time = 1)

# here is how to simulate data to see what I'm talking about

# 100 individuals 4 periods (would be years in a CJS model but phi=1 here)
nind<-100
nyear<- 4
Sx<-Sy<-matrix(NA,nrow=nind,ncol=nyear)

# simulate initial coordinates on the square:
Sx[,1]<-runif(nind,3,13)
Sy[,1]<-runif(nind,3,13)

# increment time and simulate ANGLE and RADIAL DISTANCE of movement:
for(t in 2:nyear){
phi<-runif(nind,0,360)
r<- rexp(nind,1)
dx<- r*cos(phi)
dy<- r*sin(phi)
## SEE here: NEXT LOCATION IS DERIVED
Sx[,t]<-Sx[,t-1] + dx
Sy[,t]<-Sy[,t-1]+dy
}

# now we generate encounter histories:
f<-rep(NA,nind)
Y<-matrix(0,nrow=nind,ncol=nyear)
for(i in 1:nind){
for(t in 1:nyear){
  # see here: IF individual is IN THE SQUARE we can capture it:
 if( Sx[i,t] > 3 & Sx[i,t]< 13 & Sy[i,t]>3 & Sy[i,t]<13 )
 Y[i,t]<-rbinom(1,1,.5)
}
 # f[i] = 1 period of 1st capture, NA if never captured
 if(sum(Y[i,])>0)
  f[i]<-  min((1:4)[Y[i,]==1])

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



# Likelihood 
for (i in 1:M){
w[i] ~ dbern(psi)

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
      mu2[i,t] <- p0 * inside[i,t] * w[i]
      } #t
   } #i
N<- sum(w[])
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
inits <- function(){list(psi = runif(1, 0, 1),p0 = runif(1, 0, 1), w=c(rep(1,100),rep(0,M-100)) , tau = 1) }

library("R2jags")

# Parameters monitored
parameters <- c("p0","tau","N","psi")

# Call JAGS from R
model1 <- jags(jags.data, inits, parameters, "model1.mod", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())








sink("model2.mod")
cat("
model {
psi ~ dunif(0,1)
tau ~dgamma(.1,.1)
p0 ~ dunif(0,1)

# Likelihood 
for (i in 1:M){
w[i] ~ dbern(psi)

  G[i,1,1] ~ dunif(3,13)
  G[i,1,2] ~ dunif(3,13)
  #Gmat[i,1,1]<- G1[i,1]
  #Gmat[i,1,2]<- G1[i,2]

   for (t in 2:n.occasions){

r[i,t-1] ~ dexp(tau)
phi[i,t-1] ~ dunif(0,360)
 

    dx[i,t-1]<- r[i,t-1]*cos(phi[i,t-1])
    dy[i,t-1]<- r[i,t-1]*sin(phi[i,t-1])

   # Gmat[i,t,1] <-  G[i,t-1,1] + dx[i,t-1]
   # Gmat[i,t,2] <-  G[i,t-1,2] + dy[i,t-1]


      G[i,t,1] ~ dsum(G[i,t-1,1], dx[i,t-1]) 
      G[i,t,2] ~ dsum(G[i,t-1,2], dy[i,t-1]) 

 ####   G2x = G1x + r*cos(phi)

      # Observation process
      # Test whether the actual location is in- or outside the study area. Needs to be done for each grid cell
    }
   for(t in 1:n.occasions){ 
      inside[i,t] <- step(G[i,t,1]-3) * step(13-G[i,t,1]) *step(G[i,t,2]-3) * step(13-G[i,t,2])
         
      Y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p0 * inside[i,t] * w[i]
      } #t
   } #i
N<- sum(w[])
}
",fill = TRUE)
sink()



# MCMC settings
ni <- 500
nt <- 1
nb <- 300
nc <- 2

# Bundle data
jags.data <- list(Y = Y, G=G,M=M, n.occasions = ncol(Y))
inits <- function(){list(psi = runif(1, 0, 1),p0 = runif(1, 0, 1), w=c(rep(1,100),rep(0,M-100)) , tau = 1) }

library("R2jags")

# Parameters monitored
parameters <- c("p0","tau","N","psi")

# Call JAGS from R
model1 <- jags(jags.data, inits, parameters, "model2.mod", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

















 I'm at TWS this week talking with people about SCR and related things and also
 thinking about the dispersal model you've been fitting (I have your code for 
the shrike data which I am about to run).

 I recall that the "problem" that you had hoped to solve was to be able to put 
some arbitrary distribution on "distance dispersed".  Forgive me if we talked about 
this particular solution that i'm about to propose, and maybe it was ineffective 
in some way, but I don't remember.

 

 Anyhow, you have data (x(t),y(t)) being the coordinates of the individual at time t.  
You can transform these data into, instead of coordinates, distance and angle dispersed,
 say (r,phi) the two sets of data being related by (I think):


  x(t) = r(t)*cos(phi(t))

 y(t) = r(t)*sin(phi(t))

 

(quite possibly I remember my trig stuff incorrectly).

 
Anyhow, if you do this then you can specify any distribution for r(t) phi(t) such as:

 r(t) ~ dexp(theta)

 phi(t) ~ dunif(0,360) 


 

