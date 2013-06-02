
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

