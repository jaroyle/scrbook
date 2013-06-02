#### Design optimisation with MINK data

## Part 1:  Makes up a discrete state-space for the SCR model and some potential
#           design points that exist as clusters. The matrix "clusters" is a cluster
#           identity matrix that maps EACH potential design point to a cluster
#           [note: points can be in multiple clusters]
## Part 2:  This is just an R function SCRdesign() which, with certain arguments
#           specified, will compute the optimal design.
## Part 3:  This shows how to call the function. Execute all of this and check
#           out the cool result.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### SET THE WD HERE #####

setwd("V:/scrDesignCHRIS")

install.packages(c("cluster","raster","rgdal"),repos="http://cran.us.r-project.org")

library(scrbook)
library(cluster)
library(raster)
library(rgdal)
library(MASS)
if(!"water" %in% ls()){
water 	<- shapefile("Streams_clipped.shp")
water2 	<- water[water@data$Study_Area=="Hudson",]
}
#############
## Part 1

 # can make up some data:
  lower=9; upper=21
  delta<- .4 ##1/3   # spacing of design points
  pts<-seq(lower,upper,delta)
  s1<- sort(rep(pts,length(pts)))
  s2<-rep(pts,length(pts))
  S<-cbind(s1,s2)
  inX<- (s1 <=20 & s1 >=10) & (s2 <=20 & s2>=10)
  C<-S[inX,]  # candidate set
  avail<- 1:nrow(C)
  clusters<-NULL
 repeat{
  if(all(is.na(avail))) break
   pt<- avail[!is.na(avail)][1]
   cat(pt)
   d<- as.vector(e2dist( matrix(C[pt,],ncol=2),C ))
   d[is.na(avail)]<-NA
   inclust<- order(d)
   inclust<- inclust[1:min(4)]
   avail[inclust]<-NA
   thisclust<-rep(0,nrow(C))
   thisclust[inclust]<-1
   clusters<-rbind(clusters,thisclust)
 }

 # or can use MINK data
  source("Hudson_data.Rd")
  S <- G
  C <- sites[,c("X","Y")]/1000
  clusters <- matrix(NA,nrow=length(unique(sites$Clust)),ncol=nrow(C))
  rownames(clusters) <- unique(sites$Clust)
  colnames(clusters) <- sites$Trap_ID
   for(i in 1:nrow(clusters)){
    x <- rownames(clusters)[i]
    clusters[i,] <- sites$Trap_ID %in% sites$Trap_ID[sites$Clust==x] + 0
   }
  clust.mids <- data.frame(Cl=unique(sites$Clust),midX=NA,midY=NA)
   for ( i in clust.mids$Cl){
    x<- sites[which(sites$Clust%in%i),]
    clust.mids[which(clust.mids$Cl%in%i),2:3] <- pam(cbind(x$X,x$Y),1)$medoids/1000
   }

  par(mar=c(0,0,0,0),oma=c(0,0,0,0),mfrow=c(1,1))
  plot(water2,col=4,asp=1)
  points(S*1000,pch=".",asp=1)
  points(C*1000,pch=21,bg=1,cex=2)
  points(clust.mids[,-1]*1000,pch=21,cex=1,bg=2)

#############
## PART 2
 #S=S;C=C;clusters=clusters;ntraps=9;ndesigns=10;nn=19;sigma=0.6;beta0=2;crit=3
 SCRdesign<-function(S=S,C=C,clusters=NULL,ntraps=9,ndesigns=10,nn=19,beta0=0.6,sigma=2,crit=3){
 # ntraps = clusters in this case
 ngrid<-nrow(S)

 # if trap locations only i.e NOT CLUSTER CLUSTERS
 if(is.null(clusters)){
   Cd <- round(e2dist(C,C),8)  # distance among candidates
   NN2 <- NN <- matrix(0,nrow=nrow(Cd),ncol=ncol(Cd))
  for(i in 1:nrow(Cd)){
    xx<-Cd[i,]
    NN[i,]  <- (xx>0 & xx<= sort(xx)[nn]) # select nearest nn neighbors
    NN2[i,] <- (xx>0 & xx<= sort(xx)[3])  # select nearest 3 neighbors
  }
 }

 # if using clusters
 if(!is.null(clusters)){
   Cd <- round(e2dist(clust.mids[,-1],clust.mids[,-1]),8)
   NN2 <- NN <- matrix(0,nrow=nrow(Cd),ncol=ncol(Cd))
  for(i in 1:nrow(Cd)){
    xx<-Cd[i,]
    NN[i,] <- (xx>0 & xx<= sort(xx)[nn])
    NN2[i,] <- (xx>0 & xx<= sort(xx)[3])
  }
 }

## S as defined above is ALL locations in the state-space.
## Seems like the MC average of populations of some size N=100 say,
## should be the same. Lets test that.

# X: subset of traps!
X <- C
# S: the state space!
S <- S
# N: n individuals/home range centers
N <- 100

Qfn<-function(X,S,N=100){

# Computes expected information and then inverts that
 beta0 <- beta0
 beta1 <-  -1*( 1/(sigma^2) )
 dmat  <-e2dist(X,S)^2   # ntraps x nstatespace points
sum(is.na(dmat))
# Design matrix
 lammat <- exp(beta0 + beta1 * dmat)                    # capts per trap | s
 lamvec <- exp(beta0 + beta1 * dmat[1:length(dmat)])    # c(lammat)
 lamJ <- as.vector( t(lammat)%*%rep(1,nrow(X)))         # sum of pr(c) across all traps | s
 pbar <- as.vector( 1-exp(-t(lammat)%*%rep(1,nrow(X)))) #
 #or pbar <- as.vector(1 - exp(-lamJ))                  # INCLUSION prob of AC(s)


#                        # visualise this:
#                        par(mfrow=c(1,3))
#                           x <- as.matrix(lammat[1:20,])
#                           pbar <- as.vector( 1-exp(-t(x)%*%rep(1,nrow(x)))) #
#                           rbPal <- colorRampPalette(c('blue','grey95'))
#                           cc<-rbPal(50)[as.numeric(cut(1-pbar,breaks = 50))]
#                           plot(S,pch=".",asp=1,cex=6,col=cc)
#                           points(C[1:20,],pch=".",cex=5,col=2)
#                           legend("top",paste(round(mean(pbar),2)),cex=2,bty="n")
#
#                           x <- as.matrix(lammat[(nrow(C)-20):nrow(C),])
#                           pbar <- as.vector( 1-exp(-t(x)%*%rep(1,nrow(x)))) #
#                           rbPal <- colorRampPalette(c('blue','grey95'))
#                           cc<-rbPal(50)[as.numeric(cut(1-pbar,breaks = 50))]
#                           plot(S,pch=".",asp=1,cex=6,col=cc)
#                           points(C[(nrow(C)-20):nrow(C),],pch=".",cex=5,col=2)
#                           legend("top",paste(round(mean(pbar),2)),cex=2,bty="n")
#
#                           v <- sample(1:nrow(C),20)
#                           x <- as.matrix(lammat[v,])
#                           pbar <- as.vector( 1-exp(-t(x)%*%rep(1,nrow(x)))) #
#                           rbPal <- colorRampPalette(c('blue','grey95'))
#                           cc<-rbPal(50)[as.numeric(cut(1-pbar,breaks = 50))]
#                           plot(S,pch=".",asp=1,cex=6,col=cc)
#                           points(C[v,],pch=".",cex=5,col=2)
#                           legend("top",paste(round(mean(pbar),2)),cex=2,bty="n")
#
 pbar<-mean(pbar)

 M1  <- rep(1,ntraps*nrow(S)) # DM[,1]: the intercept
 M2  <- dmat[1:length(dmat)]  # DM[,2]: X (in Y = b_0 + b_1*X)

 # Elements of the Information matrix
 I11 <- (1/nrow(S))*sum(lamvec)          # Y
 I12 <- (1/nrow(S))*sum(lamvec*M2)       # YX         | Y  | YX  |
 I21 <- (1/nrow(S))*sum(lamvec*M2)      # YX         | YX | YXX |
 I22 <- (1/nrow(S))*sum(lamvec*M2*M2)    # YXX
 # This is expected information. Probs need to multiply by E[n]
 I <- matrix(c(I11,I12,I21,I22),nrow=2,byrow=TRUE)
 I <- N*pbar*I
#  V is the var-cov matrix of the MLE
#  V <- solve(I)
  V <- solve(I,tol=1e-23)
#  V <- qr.solve(I,tol=1e-23)



 Q1 <- sum(diag(V))

 # matrix^2  * scalar
 sumsJ <- as.vector(
            t( lammat*lammat*(diag(V)[1] + (dmat^2) * diag(V)[2] -2*dmat*V[1,2]) ) %*%rep(1,nrow(X)))
 var.pbar <- ((1/nrow(S))^2)*sum(exp(-lamJ)*exp(-lamJ)*sumsJ) # Q4
 part1 <- (N*N*var.pbar)  ### /(pbar*pbar)
 part2 <-  N*(1-pbar)/pbar
 total <-  part1+part2
 newpart2 <- N*(1-pbar)*(var.pbar + 1)/pbar

 ## crit = 4 min var(beta-hat)
 ## crit = 5, min variance N-hat
 ## crit = 6, maximizes pbar
 ## crit = 7 minimizes var(pbar)
 old <- N*N*var.pbar + newpart2
 fixed <-  N*pbar*( (1-pbar) + N*pbar)*( var.pbar/(pbar^4) )  # Q2

 Q1 <- part1
 Q2 <- newpart2
 Q3 <- total
 Q4 <- Q1
 Q5 <- fixed
 Q6 <- 1 - pbar
 Q7 <- var.pbar
 #c(Q1,Q2,Q3,Q4)
 c(Q1,Q2,Q3,Q4,Q5,Q6,Q7)
}


#######



 Dlist <- list()
 Qhistory <- NULL

 for(m in 1:ndesigns){
  Qbest<-10^10

  ###############################
  # Pick a random starting design
  ##

  # for either sites:
  if(is.null(clusters)){
   X.current <- sample( 1:nrow(C),ntraps)
   X<-C[X.current,]
  }
  # or clusters:
  else{
   X.current <- sample(1:nrow(clusters),ntraps)         # pick 'ntraps' clusters
   which.sites <- apply(clusters[X.current,],2,sum)>0   # which sites are is the design
   X <-C[which.sites,]                                  # traps XY
  }
  Q <- Qfn(X,S)[crit]             # evaluate Qfn|X
  Qhistory <- c(Qhistory,Q)       # store Q
  cat("Initial Q: ",Q,fill=TRUE)  # print starting Q
  if(is.nan(Q)){
   Dlist[[m]]<- list(Q=NA, X=X,X.current=X.current)
   next # go back to the top and start the 'm' loop again
  }

  #############################
  # swap sites with neighbours
  ##

  repeat{ # repeat until 'break' is satisfied i.e. [Qbest == Q]

   for(i in 1:ntraps){                      # here we will replace each site with [nn] alternatives
     # chk is a vector of 1s for sites that can be swapped with a focal site/cluster
     chk <- NN[X.current[i],]               # chk is a vector of [1 = a near neighbour]
     chk[X.current] <- 0                    # remove any sites already in the design
     x.consider <- (1:ncol(NN))[chk==1]     # these are the alternatives to consider
     qtest <- rep(10^10,length(x.consider))
    if(length(x.consider)>0){
    for(j in 1:length(x.consider)){
    # now swap each trap with all alternative to find the best one
      Xtest.current <- X.current         # this is just a list of clusters thats all.
      Xtest.current[i] <- x.consider[j]  # switch focal with alternative site
      which.test <- apply(clusters[Xtest.current,],2,sum)>0 # new design
      Xtest <- C[which.test,]  ## new design having substituted cluster i with cluster j
  points(C*1000,pch=21,bg="grey",cex=1)
  points(Xtest*1000,pch=21,bg="red",cex=1)
      xxxx <- Qfn(Xtest,S)     # evaluate Qfn|NEW_X
      qtest[j] <- xxxx[crit]   # store the criteria of interest
    }} else {qtest <- NaN}
    # if there are any NAs, ignore and go back to the top (why??)
    if(any(is.nan(qtest))){
     Dlist[[m]] <- list(Q=NA, X=X,X.current=X.current)
     next
    }
    # if Q is better, then make the change permanent
    if(min(qtest) < Q){
      Q <- min(qtest)                    # best switch
      kp <- qtest==min(qtest)            # which switch resulted in the best swith
      X.current[i] <- x.consider[kp][1]  # make the switch permanent
      # now to deal with clusters
     if(is.null(clusters)){
       X <- C[X.current,]
     }
     else{
       which.sites <- apply(clusters[X.current,],2,sum)>0
       X <- C[which.sites,]
     }
    cat("new Q: ",Q,fill=TRUE)
    #plot(S,pch=".")
    }
   }
  cat("Current value after all J points: ", Q, fill=TRUE)
  if(Qbest == Q){                         #
  break                                   #
  }                                       #
  if(Q<Qbest) Qbest<-Q                    #
  if(Q>Qbest) cat("ERROR",fill=TRUE)      ### checks whether and better switches can be made
  Qhistory<-c(Qhistory,Q)
  }

 Dlist[[m]]<- list(Q = Qbest, # final Q value
                   X = X,     # starting desgin
                   X.current = X.current # final/best design
                   )
 m <- m+1
 }
 # post processing
 Qvec <- rep(NA,length(Dlist))
 Xid <- matrix(NA,nrow=ntraps,ncol=length(Dlist))
 Xlst <- list()

  for(i in 1:length(Dlist)){
   Qvec[i] <- Dlist[[i]]$Q
   Xid[,i] <- Dlist[[i]]$X.current
   Xlst[[i]] <- Dlist[[i]]$X
  }
 od <- order(Qvec)
 tmp <- list()
  for(i in 1:length(od)){
   tmp[[i]] <- Xlst[[od[i]]]
  }
 Qvec <- Qvec[od]
 Xid <- Xid[,od]
 Xlst <- tmp
 output <- list(Qvec = Qvec,Xid = Xid,Xlst = Xlst,C = C,S = S,Qhistory = Qhistory)
 return(output)
 }

 ### PART 3
 ###
 ###

 ## crit=4  = trace(V(alpha))
 ## crit 5 = Var(N)
 ## crit 6 = 1-pbar
 ## crit 7 = var(pbar)

 # We can do this for cluster but to get the optimum number of sites (90) I must
 # loop over several ntrap values (e.g. 32:36)
 plot(water2,col=4,asp=1)
 nclusts <- 29:31
 des.Clust <- list()
 # here I will use:  ntraps == a range of values
 #                   ndesign == 2 (more is better but how many more??)
 #                   nn == 5
 #                   sigma == 0.64
 #                   beta0 == 2.53 **from the poisson encounter frequency**
 #                   crit == 5
 for(i in 1:length(nclusts)){
  des.Clust[[i]]<-SCRdesign(S,C,clusters=clusters,ntraps=nclusts[i],ndesigns=10,nn=20,
                          sigma=0.64,beta0=2.53, crit=5)
  }

 for(j in 1:length(des.Clust)){
  for(i in 1:length(des.Clust[[1]]$Xlst)){
   tmpmerge <- merge(des.Clust[[j]]$Xlst[[i]]*1000,sites,by=c("X","Y"))
   print(length(factor(unique(tmpmerge$Site_ID)))/(31+j))
  }
 }


## Process the data
#plot(des.Clust[[1]]$Qvec,type="l")
#for(i in 1:5) lines(des.Clust[[i]]$Qvec,type="l")
#
#par(mfrow=c(5,4),oma=c(0,0,0,0),mar=c(0,0,0,0))
#for(i in 1:20){
# plot( des.Clust[[5]]$C$X,des.out[[5]]$C$Y,
#       asp=1,pch=21,bg="grey",xlab="",ylab="",axes=F)
# points(des.Clust[[5]]$Xlst[[i]]$X,des.Clust[[5]]$Xlst[[i]]$Y,
#        pch=21,bg=2)
# box(bty="o")
# legend("topleft",paste("Q = ",round(des.Clust[[5]]$Qvec[i],3),sep=""),cex=1.5,pch="",bty="n")
#}
#