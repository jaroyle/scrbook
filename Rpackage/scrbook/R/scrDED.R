



# Fit poisson SCR model with density and ecological-distance covariates
scrDED <- function(y=y, traplocs=traplocs,
                      den.formula=~1, dist.formula=~1,
                      rasters, start, transform=TRUE, ...) {
    if(!require(raster))
        stop("raster package must be loaded")
    if(!require(gdistance))
        stop("gdistance package must be loaded")
    if(! class(rasters)[1] %in% c("RasterLayer", "RasterStack"))
        stop("rasters should have class RasterLayer or RasterStack")
    if(any(rowSums(y)==0))
        stop("you cannot observer an all 0 capture history")

    # do a check here that trap locations exist in same space as raster.

    dims <- dim(rasters)
    rp <- as.data.frame(rasterToPoints(rasters))
    G <- as.matrix(rp[,1:2]) # state-space pixel centers

    res <- res(rasters)
    area <- prod(res)
    npix <- ncell(rasters)
    SSarea <- area*npix # check for rasterStacks
    ext <- extent(rasters)

    isEuclid <- isTRUE(all.equal(dist.formula, ~1))
    if(isEuclid)
        D <- e2dist(traplocs, G)

    den.vars <- all.vars(den.formula)
    dist.vars <- all.vars(dist.formula)

    if(length(den.vars)>0 && (!den.vars %in% names(rasters)))
        stop("variables in den.formula must occur in names(rasters)")
    if(length(dist.vars)>0 && (!dist.vars %in% names(rasters)))
        stop("variables in dist.formula must occur in names(rasters)")

    Xden <- model.matrix(den.formula, rp)
    Xdist <- model.matrix(dist.formula, rp)

    isInt1 <- colnames(Xden) == "(Intercept)"
    isInt2 <- colnames(Xdist) == "(Intercept)"
    if(any(isInt1))
        Xden <- Xden[,-which(isInt1),drop=FALSE]
    if(any(isInt2))
        Xdist <- Xdist[,-which(isInt2),drop=FALSE]

    if(transform & ncol(Xdist)>0) {
        Xdist <- apply(Xdist, 2, function(x) x-min(x, na.rm=TRUE))
        Xdist <- apply(Xdist, 2, function(x) x/max(x, na.rm=TRUE))
    }

    np.den <- ncol(Xden)
    np.dist <- ncol(Xdist)
    np <- 3+np.den+np.dist

    den.names <- dist.names <- character(0)
    if(np.den>0)
        den.names <- paste("den", 1:np.den, sep="")
    if(np.dist>0)
        dist.names <- paste("dist", 1:np.dist, sep="")

    if(missing(start)) {
        start <- rep(0, np)
        start[2] <- log((ext@xmax-ext@xmin)/3)
        start[3] <- log(nrow(y)/2)
    }
    if(is.null(names(start)))
        names(start) <- c("lam0", "sigma", "n0", den.names, dist.names)

    y.ij <- rbind(y, rep(0, ncol(y)))
    nry <- nrow(y.ij)

    y.ijg <- t(apply(y.ij, 1, rep, npix))

    mu <- rep(1/npix, npix)

    # Negative log-likelihood
    nll <- function(pars) {
        lam0 <- exp(pars[1])
        sigma <- exp(pars[2])
        n0 <- exp(pars[3])
        if(np.den > 0) {
            mu <- exp(Xden %*% pars[4:(3+np.den)])
            mu <- mu/sum(mu) #Gprobs
        }
        if(!isEuclid) {
            cost <- exp(Xdist %*% pars[(4+np.den):np])
            cost <- raster(matrix(cost, dims[1], dims[2], byrow=TRUE))
            extent(cost) <- ext # should add projection too
            tr1 <- transition(cost, transitionFunction = function(x)
                              1/mean(x), directions=8)
            tr1CorrC <- geoCorrection(tr1, type="c",
                                      multpl=FALSE, scl=FALSE)
            D <- costDistance(tr1CorrC, traplocs, G)
        }

        probcap <- as.numeric(lam0*exp(-D*D/(2*sigma*sigma)))
        lik.marg <- rep(NA_real_, nrow(y.ij))
        for(i in 1:nry) {
            Pm <- dpois(y.ijg[i,], probcap, log=TRUE)
            Pm <- matrix(Pm, nrow(D), ncol(D))
            lik.cond <- exp(colSums(Pm))
            lik.marg[i] <- sum(lik.cond*mu)
        }
        nv <- c(rep(1, length(lik.marg)-1), n0)
        part1 <- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
        part2 <- sum(nv*log(lik.marg))
        out <-  -1*(part1+ part2)
        return(out)
    }

    fm <- optim(start, nll, ...)
    return(fm)
}


