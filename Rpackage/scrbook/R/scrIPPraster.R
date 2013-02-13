

# MCMC. SCR model with inhomogenous point process and sex-secific sigma

# Like scrJag6, but habitat is categorical, not ordinal

# FIXME: Make notation consistent with chapter, e.g. w -> z


scrIPPraster <- function(Z, X, niters, habitat,
                         tune=c(5, 0.1, 0.1, 0.1, 100),
                         sex, augment=50, map.thin=1, map.burn=1000)
{
    require(raster)

    Zdims <- dim(Z)
    nind <- Zdims[1]
    R <- Zdims[2]
    T <- Zdims[3]

    M <- nind+augment
    Zaug <- array(0, c(M, R, T))
    Zaug[1:nind,,] <- Z
    for(i in (nind+1):M) {
        Zaug[i,,][is.na(Z[1,,])] <- NA
    }
    Zin <- Z
    Z <- Zaug

    sexAug <- rep(NA, M)
    sexAug[1:nind] <- sex
    rho <- 0.5
    sexAug[(nind+1):(nind+augment/2)] <- 0
    sexAug[(nind+1+augment/2):M] <- 1
    sexIn <- sex
    sex <- sexAug

    # initial values
    # FIXME: Make generic using spsample() or something
    S <- cbind(rnorm(M, 780000, 100), rnorm(M, 7140000, 100))

    plot(habitat)
    points(X, pch="+")
    points(S, pch=16, col=4)

    D <- e2dist1(S, X)
    sigmaM <- sigmaF <- runif(1, 4000, 8000)
    sigma <- rep(NA, M)
    sigma[sex==0] <- sigmaM
    sigma[sex==1] <- sigmaF
    lam0 <- runif(1, 0.0003, 0.001)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    lam[lam==0] <- 1e-50

    beta0 <- rnorm(1, -5, .5) # intercept of intensity function
    beta1 <- rnorm(1, 0, .5)  # covariate effect

    s.habitat <- extract(habitat, S)

    npix <- cellStats(!is.na(habitat), sum)
    pix.area <- prod(res(habitat)) / 1e6 # km^2
    EN <- cellStats(exp(beta0 + beta1*habitat) *
                    pix.area, sum)
    psi <- EN / M
    if(psi > 1)
        stop("psi > 1")

    w <- rbinom(M, 1, psi)
    w[rowSums(Z, na.rm=TRUE)>0] <- 1

    habitatLow <- habitat == -1
    habitatMed <- habitat == 0
    habitatHigh <- habitat == 1
    habitat.area1 <- cellStats(habitatLow, sum) * pix.area
    habitat.area2 <- cellStats(habitatMed, sum) * pix.area
    habitat.area3 <- cellStats(habitatHigh, sum) * pix.area
    habitat.area <- c(habitat.area1, habitat.area2, habitat.area3)

    which.habitat <- factor(extract(habitat, S[w==1,]),
                            levels=-1:1)

    density.habitat <- table(which.habitat) / habitat.area * 100


    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=11)
    colnames(out) <- c("sigmaF", "sigmaM", #"rho",
                       "lam0", "psi",
                       "beta0", "beta1", #"beta2",
                       "N", "rho", "D1", "D2", "D3")

    if(!missing(habitat)) {
        dmap <- habitat
        dmap[!is.na(habitat)] <- 0
    }

#    browser()

    cat("\ninitial values =",
        round(c(sigmaF, sigmaM, lam0, psi, beta0, beta1, #beta2,
                sum(w), density.habitat), 3), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 10 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("current =", round(out[iter-1,], 3), "\n")
            cat("  Acceptance rates\n")
            cat("    S =", Sups/M, "\n")
            cat("    w =", wUps/M, "\n")
        }

        ll <- sum(dpois(Z, lam*w, log=TRUE), na.rm=TRUE)

        # update sigmaM
        sigmaM.cand <- rnorm(1, sigmaM, tune[1])
        if(sigmaM.cand > 0) {
            sigma.cand <- sigma
            sigma.cand[sex==0] <- sigmaM.cand
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE), na.rm=TRUE )
            if(runif(1)<exp( llcand  - ll ) ){
                ll <- llcand
                lam <- lam.cand
                sigmaM <- sigmaM.cand
                sigma <- sigma.cand
            }
        }

        # update sigmaF
        sigmaF.cand <- rnorm(1, sigmaF, tune[1])
        if(sigmaF.cand > 0) {
            sigma.cand <- sigma
            sigma.cand[sex==1] <- sigmaF.cand
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE), na.rm=TRUE )
            if(runif(1)<exp( llcand  - ll ) ){
                ll <- llcand
                lam <- lam.cand
                sigmaF <- sigmaF.cand
                sigma <- sigma.cand
            }
        }

        # update lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE), na.rm=TRUE )
            if(runif(1) < exp( llcand - ll ) ) {
                lam0<-lam0.cand
                lam<-lam.cand
                ll <- llcand
            }
        }


        # update beta0
        EN <- cellStats(exp(beta0 + beta1*habitat) *
                        pix.area, "sum")
        beta0.cand <- rnorm(1, beta0, tune[3])
        if(beta0.cand > -10) {
        EN.cand <- cellStats(exp(beta0.cand + beta1*habitat) *
                             pix.area,"sum")
        ll.beta <- dbinom(sum(w), M, EN/M, log=TRUE) +
            sum((beta0 + beta1*s.habitat +
                 log(pix.area) - log(EN))*w)

        if(EN.cand <= M) {
            ll.beta.cand <- dbinom(sum(w), M, EN.cand/M, log=TRUE) +
                sum((beta0.cand + beta1*s.habitat +
                     log(pix.area) - log(EN.cand))*w)

            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
                beta0 <- beta0.cand
                EN <- EN.cand
                ll.beta <- ll.beta.cand
            }
        }
        }

        # update beta1
        beta1.cand <- rnorm(1, beta1, tune[4])
        if(beta1.cand > -10) {
        EN.cand <- cellStats(exp(beta0 + beta1.cand*habitat) * pix.area,
                             "sum")
        if(EN.cand <= M) {
            ll.beta.cand <- dbinom(sum(w), M, EN.cand/M, log=TRUE) +
                sum((beta0 + beta1.cand*s.habitat +
                     log(pix.area) - log(EN.cand))*w)
            if(runif(1) < exp(ll.beta.cand - ll.beta) )  {
                beta1 <- beta1.cand
                EN <- EN.cand
                ll.beta <- ll.beta.cand
            }
        }
        }


        # update psi
        psi <- EN / M

        # update w
        wUps <- 0
        seen <- apply(Z>0, 1, any, na.rm=TRUE)
        for(i in 1:M) {
            if(seen[i])
                next
            wcand<-w
            if(w[i]==0) {
                wcand[i] <- 1
                ll.w <- 0
                ll.w.cand <- sum(dpois(Z[i,,], lam[i,]*wcand[i], log=TRUE),
                                 na.rm=TRUE) #+
            } else {
                wcand[i] <- 0
                ll.w <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE),
                            na.rm=TRUE) #+
                ll.w.cand <- 0
            }
            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp((ll.w.cand+prior.cand) - (ll.w+prior))) {
                w <- wcand
                wUps <- wUps+1
            }
        }

        # update S
        Sups <- 0
        for(i in 1:M) {
            Scand <- c(rnorm(1, S[i,1], tune[5]),
                       rnorm(1, S[i,2], tune[5]))
#            if(is.na(over(SpatialPoints(matrix(Scand, 1)), poly)))
#                next
            if(is.na(extract(habitat, matrix(Scand,1))))
                next
            dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
            lam.cand <- lam
            lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma[i]^2))
            if(w[i]==0)
                ll.S <- ll.S.cand <- 0
            else {
                ll.S <- sum(dpois(Z[i,,], lam[i,], log=TRUE),
                            na.rm=TRUE)
                ll.S.cand <- sum(dpois(Z[i,,], lam.cand[i,], log=TRUE),
                                 na.rm=TRUE)
            }
            #ln(prior), denominator is constant
            s.habitat.cand <- extract(habitat, matrix(Scand,1))
            prior.S <- beta0 + beta1*s.habitat[i]
            prior.S.cand <- beta0 + beta1*s.habitat.cand

           if(runif(1)< exp((ll.S.cand+prior.S.cand) - (ll.S+prior.S))) {
               s.habitat[i] <- s.habitat.cand
                S[i,] <- Scand
                lam <- lam.cand
                D[i,] <- dtmp
                Sups <- Sups+1
            }
        }

        if(map.burn <= iter)
            points(S[w==1,], cex=0.4, pch=16, col=rgb(0,0,1,0.1))

        # Density in each habitat type
        which.habitat <- factor(extract(habitat, S[w==1,]),
                                levels=-1:1)
        density.habitat <- table(which.habitat) / habitat.area * 100

        if((map.thin>0) & (map.burn <= iter)) {
            if(iter %% map.thin == 0) {
                real.guys <- S[w==1,]
                cells <- cellFromXY(habitat, real.guys)
                val <- extract(dmap, real.guys)
                dmap[cells] <- val + 1
            }
        }

        out[iter,] <- c(sigmaF, sigmaM, #rho,
                        lam0, psi, beta0, beta1, #beta2,
                        sum(w), sum(w*(1-sex))/sum(w), density.habitat)
    }
    if(map.thin>0) {
        dmap <- dmap / ((niters-map.burn)/map.thin)
    }

    last <- list(S=S, lam=lam, w=w, sex=sex, sigma=sigma, D=D)
    list(out=out, last=last, map=dmap)
}





































