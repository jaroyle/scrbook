




scrQUAD <- function(y, X, M, raster,
                  niters, xlims, ylims, tune=c(0.5, 0.1, 2)) {

    J <- nrow(y)
    K <- ncol(y)
    s <- cbind(runif(M, xlims[1], xlims[2]),
               runif(M, ylims[1], ylims[2]))
    u <- array(NA, c(M, 2, K))
    tau <- runif(1, 2, 4)
    for(k in 1:k) {
        u[,,k] <- cbind(rnorm(M, s[,1], tau), rnorm(M, s[,2], tau))
    }

#    D <- e2dist(S, X)
#    sigma <- runif(1, 0.5, 1.5)
#    lam0 <- runif(1, 0.7, 0.9) # This is p0 when obsmod="bern"
#    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    psi <- runif(1, 0.9, 0.99)
    w <- rbinom(M, 1, psi)
    N0 <- sum(w)
    ymx <- max(y, na.rm=TRUE)
    if(N0 < ymx) {
        is0 <- which(w==0)
        w[sample(is0, ymx-N0)] <- 1
    }

    N <- y+3

    out <- matrix(NA,nrow=niters,ncol=4)
    colnames(out) <- c("tau", "p", "psi", "N")

    cat("\nstarting values =", c(tau, p, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

        # update tau
        tau.cand <- rnorm(1, tau, tune[1])
        if(tau.cand > 0) {
            # recycle over K
            llx <- dnorm(u[,1,], s[,1], tau, log=TRUE)
            lly <- dnorm(u[,2,], s[,2], tau, log=TRUE)
            ll <- sum(llx+lly)

            llx.c <- dnorm(u[,1,], s[,1], tau.cand, log=TRUE)
            lly.c <- dnorm(u[,2,], s[,2], tau.cand, log=TRUE)
            ll.c <- sum(llx.c + lly.c)

            if(runif(1) < exp(ll.c - ll)) {
                tau <- tau.cand
            }
        }

        # update p
        p.cand <- rnorm(1, p, tune[2])
        if(p.cand>=0 & p.cand<=1) {
            ll <- sum(dbinom(y, N, p, log=TRUE))
            ll.c <- sum(dbinom(y, N, p, log=TRUE))
            if(runif(1) < exp(ll.c - ll)) {
                p <- p.cand
            }
        }

        # update w
        wUps <- 0
        w.cand <- w
        N.cand <- N
        for(i in 1:M) {
            w.cand[i] <- if(w[i]==0) 1 else 0
            for(k in 1:K) {
                cells <- cellFromXY(raster, u[,,k] * w.cand)
                counts <- table(cells)
                counts.in <- counts[names(count) %in% rownames(N)]
                N.cand[names(counts.in),k] <- counts.in
            }
            ll <- sum(dbinom(y, N, p, log=TRUE))
            ll.c <- sum(dbinom(y, N.cand, p, log=TRUE))
            if(runif(1) < exp(ll.c - ll)) {
                w[i] <- w.cand[i]
                N <- N.cand
                wUps <- wUps+1
            }
        }


        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update s
        s.ups <- 0
        for(i in 1:M) {   # note this is "M" in general
            s.cand <- c(rnorm(1, s[i,1], tune[3]),
                       rnorm(1, s[i,2], tune[3]))
            inbox <- s.cand[1]>=xlims[1] & s.cand[1]<=xlims[2] &
                     s.cand[2]>=ylims[1] & s.cand[2]<=ylims[2]
            if(!inbox)
                next
            llx <- dnorm(u[i,1,], s[i,1], tau, log=TRUE)
            lly <- dnorm(u[i,2,], s[i,2], tau, log=TRUE)
            ll <- sum(llx+lly)

            llx.c <- dnorm(u[i,1,], s.cand[1], tau, log=TRUE)
            lly.c <- dnorm(u[i,2,], s.cand[2], tau, log=TRUE)
            ll.c <- sum(llx.c+lly.c)

            if(runif(1) < exp(ll.c - ll)) {
                s[i,] <- s.cand
                s.ups <- s.ups+1
            }
        }

        # update u
        u.ups <- 0
        u.cand <- u #array(NA, c(M, 2, K))
        N.cand <- N
        for(i in 1:M) {
            for(k in 1:K) {
                u.cand[i,,k] <- c(rnorm(1, s[i,1], tune[3]),
                                  rnorm(1, s[i,2], tune[3]))
                cells <- cellFromXY(raster, u[,,k] * w)
                counts <- table(cells)
                counts.in <- counts[names(count) %in% rownames(N)]
                N.cand[names(counts.in),k] <- counts.in
            }
            ll <- sum(dbinom(y, N, p, log=TRUE))
            ll.c <- sum(dbinom(y, N.cand, p, log=TRUE))

            priorx <- dnorm(u[i,1,], s[i,1], tau, log=TRUE)
            priory <- dnorm(u[i,2,], s[i,2], tau, log=TRUE)
            prior <- sum(priorx + priory)

            priorx.c <- dnorm(u.cand[i,1,], s[i,1], tau, log=TRUE)
            priory.c <- dnorm(u.cand[i,2,], s[i,2], tau, log=TRUE)
            prior.c <- sum(priorx.c + priory.c)

            if(runif(1) < exp((ll.c+prior.c) - (ll + prior))) {
                u[i,,] <- u.cand[i,,]
                u.ups <- u.ups+1
                N <- N.cand
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     Z =", round(Zups/(J*K), 2), "\n")
            cat("     s =", s.ups/M, "\n")
        }

        out[iter,] <- c(tau, p, psi, sum(w) )

    }
    return(out)
}


