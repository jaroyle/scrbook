intlik3 <-
function (start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs, 
          ssbuffer = 2,model="B",predict=FALSE) 
{
    Xl <- min(X[, 1]) - ssbuffer
    Xu <- max(X[, 1]) + ssbuffer
    Yu <- max(X[, 2]) + ssbuffer
    Yl <- min(X[, 2]) - ssbuffer
    SSarea <- (Xu - Xl) * (Yu - Yl)
    if (is.null(K)) 
        return("need sample size")
    xg <- seq(Xl + delta/2, Xu - delta/2, delta)
    yg <- seq(Yl + delta/2, Yu - delta/2, delta)
    npix.x <- length(xg)
    npix.y <- length(yg)
    area <- (Xu - Xl) * (Yu - Yl)/((npix.x) * (npix.y))
    G <- cbind(rep(xg, npix.y), sort(rep(yg, npix.x)))
    nG <- nrow(G)
    D <- e2dist(X, G)
    if (is.null(start)) 
        start <- c(0, 0, 0)
    alpha0 <- start[1]
    alpha1 <- exp(start[2])
    n0 <- exp(start[3])
    if(model=="B")
        probcap <- plogis(alpha0) * exp(-alpha1 * D * D)
    if(model=="P")
       probcap <- exp(alpha0) * exp(-alpha1 * D * D)
        Pm <- matrix(NA, nrow = nrow(probcap), ncol = ncol(probcap))
    ymat <- y
    ymat <- rbind(y, rep(0, ncol(y)))
    lik.marg <- rep(NA, nrow(ymat))
    for (i in 1:nrow(ymat)) {
        if(model=="B")
        Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap[1:length(Pm)], log = TRUE))
        if(model=="P")
        Pm[1:length(Pm)] <- (dpois(rep(ymat[i, ], nG), rep(K, nG)*probcap[1:length(Pm)], log = TRUE))
        lik.cond <- exp(colSums(Pm))
        lik.marg[i] <- sum(lik.cond * (1/nG))
    }
    if(predict==FALSE){
    nv <- c(rep(1, length(lik.marg) - 1), n0)
    part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
    part2 <- sum(nv * log(lik.marg))
    out <- -1 * (part1 + part2)
    attr(out, "SSarea") <- SSarea
    return(out)
    }
    if(predict==TRUE){
        
        posterior<-matrix(NA,nrow=nG,ncol=nrow(ymat))        
    for(i in 1:nrow(ymat)){
    if(model=="B")
        Pm[1:length(Pm)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap[1:length(Pm)], log = TRUE))
    if(model=="P")
        Pm[1:length(Pm)] <- (dpois(rep(ymat[i, ], nG), rep(K, nG)*probcap[1:length(Pm)], log = TRUE))
         
    lik.cond <- exp(colSums(Pm))*(1/nG)
    posterior[,i]<- lik.cond/lik.marg[i]
}

return(cbind(G,posterior))
}


}
