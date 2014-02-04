
shad.cjs<-function(model=c('NS-CJS', 'MS-CJS', 'S-CJS'), n.chains, n.adapt, n.iter) {

mod<-match.arg(model)

#install required packages
library(rjags)
library(reshape)
data(Ch16shaddata)

##data setup and organization
# Input and format data matrix for spatial CJS
shad <- Ch16shaddata$shad
melted.rkm <- melt(shad, id=c("TagID","RKM")) 
tagid.week <- cast(melted.rkm, TagID ~ RKM ~ value, fill=0, length)

#data for non-spatial CJS
hold=as.matrix(table(shad$TagID, shad$Week))
hold[hold > 1] <- 1
y <- hold

#data for Multistate model

melt.state <- melt(shad, id=c("TagID","state", "Week")) 
shad.state<- cast(melt.state, TagID ~ Week ~ state, fill=0, length)

newy=matrix(NA, 315,12)
for(i in 1:315){
for(j in 1:12){
a=which(shad.state[i,j,] == max(shad.state[i,j,]))
newy[i,j] <- as.numeric(a[1])
if(max(shad.state[i,j,]) == 0) newy[i,j] <- 3
}
}

first<-Ch16shaddata$first


# Constants:
M <- 315       # Number of individuals
T <- 12     # Number of periods (weeks)

#set up trap locations and state space
nantenna <- 7  # weir, 6 antennas
antenna.loc <- c(3.7, 7.7, 13.4, 45.3, 56.4, 72.0, 77.0)  # antenna locations
xl <- 0     # lower boundary, river mouth
xu <- 82    # upper boundary, Atkinson Mill Dam



#######################################################################################################
###### Multistate-CJS  MS-CJS #########################################################################


if(mod=='MS-CJS') {

sink("MultistateCJS.txt")
cat("

model { #model

for(r in 1:2){
   phi[r] ~ dunif(0,1)
   psi[r] ~ dunif(0,1)
     p[r] ~ dunif(0,1)
}

for (i in 1:M){    
    z[i,first[i]] <- y[i, first[i]]

for (t in (first[i]+1):T){
    z[i,t] ~ dcat(ps[z[i,t-1], i, ])
    y[i,t] ~ dcat(po[z[i,t], i, ])
}

ps[1, i, 1] <- phi[1] * (1-psi[1])
ps[1, i, 2] <- phi[1] * psi[1]
ps[1, i, 3] <- 1-phi[1]
ps[2, i, 1] <- phi[2] * (1-psi[2])
ps[2, i, 2] <- phi[2] * psi[2]
ps[2, i, 3] <- 1-phi[2]
ps[3, i, 1] <- 0
ps[3, i, 2] <- 0
ps[3, i, 3] <- 1

po[1, i, 1] <- p[1]
po[1, i, 2] <- 0
po[1, i, 3] <- 1-p[1]
po[2, i, 1] <- 0
po[2, i, 2] <- p[2]
po[2, i, 3] <- 1-p[2]
po[3, i, 1] <- 0
po[3, i, 2] <- 0
po[3, i, 3] <- 1
}

} #model

", fill = TRUE)
sink()

######################################################

#Set up data input
data<-list(y=newy, first=first, M=M, T=T)

#create intial values for z state
z=matrix(NA, M, T)
for(i in 1:M){ 
for(t in (first[i]+1):T){ 
z[i,t] <-newy[i,t]
} 
}
z[z == 3] <- 2

#Set initial values
inits =  function() {list(z=z, phi=runif(2,.6,1), psi=runif(2,0,1), p=runif(2,0.7,1)) }

# Parameters to follow
parameters <- c("psi", "phi", "p")

modelFile="MultistateCJS.txt"                                                  

} #end if

######################################################################################################
###### Non-Spatial Cormack Jolly Seber   NS-CJS#######################################################


if(mod=='NS-CJS') {
sink("ModelNSCJS.txt")
cat("

model { #model
# Priors

phi ~ dunif(0,1)   # Survival (constant)

for(t in 1:T){
p[t] ~ dunif(0, 1)    # detection (constant)
}
    
for (i in 1:M){ #m    

   z[i,first[i]] ~ dbern(1)  # Known to be alive at entry into study area

  for (t in (first[i]+1):T) { #t
         tmp[i,t] <- z[i,t]*p[t]
           y[i,t] ~ dbern(tmp[i,t])
       phiUP[i,t] <- z[i,t-1]*phi     
           z[i,t] ~ dbern(phiUP[i,t])
         
   } # t
  } # m
} #model


", fill = TRUE)
sink()


######################################################

#Set up data input
data<-list(y=y, first=first, M=M, T=T)

z=matrix(NA, M, T)
for(i in 1:M){ 
for(t in first[i]:T){ 
z[i,t] <-1
} 
}

#Set initial values
inits =  function() {list(z=z, phi=runif(1,0,1), p=runif(12,0,1)) }

# Parameters to follow
parameters <- c("phi", "p")
                                                                
modelFile= "modelNSCJS.txt"                                                  

} #end if

######################################################################################################
###### Spatial Cormack Jolly Seber   S-CJS#######################################################


if(mod=='S-CJS') {

sink("SpatialCJS.txt")
cat("

model {
# Priors
sigma ~ dunif(0,80)
sigma2 <- sigma*sigma
phi ~ dunif(0, 1)   # Survival (constant across time)
tauv~dunif(0, 30)
tau<-1/(tauv*tauv)

for(t in 1:T){
lam0[t] ~ dgamma(0.1, 0.1)
}

for (i in 1:M){
  z[i,first[i]] <- 1
  S[i,first[i]] ~ dunif(0,50)

for(j in 1:nantenna) {
  D2[i,j,first[i]] <- pow(S[i,first[i]]-antenna.loc[j], 2)
       lam[i,j,first[i]]<-  lam0[first[i]]*exp(- D2[i,j,first[i]]/(2*sigma2))
       tmp[i,j,first[i]] <- lam[i,j,first[i]]
         y[i,j,first[i]] ~ dpois(tmp[i,j,first[i]])
      }

   for (t in first[i]+1:T) {
          S[i,t] ~ dunif(xl, xu)
         for(j in 1:nantenna) {
                D2[i,j,t] <- pow(S[i,t]-antenna.loc[j], 2)
               lam[i,j,t] <- lam0[t] * exp(-D2[i,j,t]/(2*sigma2))
               tmp[i,j,t] <- z[i,t]*lam[i,j,t]
                 y[i,j,t] ~ dpois(tmp[i,j,t])
 }
    phiUP[i,t] <-  z[i,t-1]*phi
       z[i,t] ~ dbern(phiUP[i,t])
}
}
}
", fill = TRUE)
sink()

######################################################

#Set up data input

#Set up a data input
data<-list(y=tagid.week, first=first, M=M, T=T, nantenna=nantenna, antenna.loc=antenna.loc)

z=matrix(NA, M, T)
St=matrix(NA, M, T)

for(i in 1:M){ 
for(t in first[i]:T){ 
z[i,t] <-1
St[i,t] <- (sum(y[i,,t] * antenna.loc)/sum(y[i,,t])) + 3
if(St[i,t] == "NaN") St[i,t] <- St[i,t-1]
} 
}

#Set initial values
inits =  function() {list(z=z, phi=runif(1,0,1), lam0=runif(12,1,3), tauv=runif(1,0,20), mu=runif(1,1,3), beta=runif(1,0,1)) }

# Parameters to follow
parameters <- c("sigma", "phi", "lam0")

modelFile= "SpatialCJS.txt"                                                  

} #end if

######################################################################################################

#run model
mod.out <- jags.model(modelFile, data, inits, n.chains=n.chains, n.adapt=n.adapt)
out <- coda.samples(mod.out,  parameters, n.iter=n.iter)

return(out)

} #end function
