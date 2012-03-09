




scrQUAD <- function(y, X, M, raster,
                  niters, xlims, ylims, tune=c(0.5, 0.1, 2, 0.2)) {

    J <- nrow(y)
    K <- ncol(y)

    psi <- runif(1, 0.7, 0.9)
#    w <- rbinom(M, 1, psi)
    w <- rep(0L, M)
#    w[1:nrow(s)] <- 1
    w[] <- 1

#    N <- ifelse(y==0, y, y+2)
    N <- matrix(0, nrow(y), ncol(y))
    rownames(N) <- rownames(y) # Important

    tau <- runif(1, 1, 2)
#    tau <- 2
    p <- runif(1, 0.4, 0.8)

    # Cheat by using original true s matrix
    cheat <- nrow(s)
    s <- rbind(s, matrix(runif((M-cheat)*2, xlims[1], xlims[2]),
                         M-cheat, 2, byrow=TRUE))
#    s <- cbind(runif(M, xlims[1], xlims[2]),
#               runif(M, ylims[1], ylims[2]))

    plot(s, xlim=xlims, ylim=ylims, col=4, pch=16)
    points(X, pch="+")

    cells <- matrix(NA, M, K)

    # cheat by using original u
#    u <- array(NA, c(M, 2, K))
    u2 <- array(NA, c(M, 2, K))
    u2[1:cheat,,] <- u
    for(i in (cheat+1):M) {
        for(k in 1:K) {
            u2[i,,k] <- rnorm(2, s[i,], tau)
        }
    }
    u <- u2

    for(k in 1:K) {
        cat("finding starting values for u[,,", k, "]...\n", sep="")
        points(u[,,k], cex=0.5, pch=16, col=3)
        while(any(N[,k] < y[,k])) {
            N[,k] <- 0
#            u[,,k] <- cbind(rnorm(M, s[,1], tau), rnorm(M, s[,2], tau))
#            points(u[,,k], cex=0.5, pch=16, col=3)

            cells[,k] <- cellFromXY(raster, u[,,k])
            cells[,k] <- ifelse(w==1, cells[,k], -1*cells[,k])

            counts <- table(cells[,k])
            counts.in <- counts[names(counts) %in% rownames(N)]
            N[names(counts.in),k] <- counts.in
        }
    }
    if(any(N < y))
        stop("doh")

    ## for(j in 1:J) {
    ##     dist <- sqrt((s[,1]-X[j,1])^2 + (s[,2]-X[j,2])^2)
    ##     prob <- exp(-dist^2/exp(2*tau^2)) * w
    ##     for(k in 1:K) {
    ##         if(N[j,k]==0)
    ##             next
    ##         where <- sample(1:M, N[j,k], prob=prob)
    ##         u[where,,k] <- X[j,]
    ##     }
    ## }


    llu <- llu.c <- llu1 <- llu2 <- llu1.c <- llu2.c <- matrix(NA, M, K)
    lly <- lly.c <- matrix(NA, J, K)

    out <- matrix(NA,nrow=niters,ncol=4)
    colnames(out) <- c("tau", "p", "psi", "N")

    cat("\nstarting values =", c(tau, p, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 10 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

        for(k in 1:K) {
            # recycle over K
            llu1[,k] <- dnorm(u[,1,k], s[,1], tau, log=TRUE)
            llu2[,k] <- dnorm(u[,2,k], s[,2], tau, log=TRUE)
            llu[,k] <- llu1[,k] + llu2[,k]
        }


        # update tau
        tau.cand <- rnorm(1, tau, tune[1])
        if(tau.cand > 0) {
            for(k in 1:K) {
                llu1.c[,k] <- dnorm(u[,1,k], s[,1], tau.cand, log=TRUE)
                llu2.c[,k] <- dnorm(u[,2,k], s[,2], tau.cand, log=TRUE)
                llu.c[,k] <- llu1.c[,k] + llu2.c[,k]
            }
            if(runif(1) < exp(sum(llu.c) - sum(llu))) {
                tau <- tau.cand
                llu1 <- llu1.c
                llu2 <- llu2.c
                llu <- llu.c
            }
        }

        # update p
        lly[] <- dbinom(y, N, p, log=TRUE)

        p.cand <- rnorm(1, p, tune[2])
        if(p.cand>=0 & p.cand<=1) {
            lly.c[] <- dbinom(y, N, p.cand, log=TRUE)
            if(runif(1) < exp(sum(lly.c) - sum(lly))) {
                p <- p.cand
                lly <- lly.c
            }
        }

        # update w
        wUps <- 0
        w.cand <- w
        N.cand <- N
        for(i in 1:M) {
            w.cand[i] <- if(w[i]==0) 1 else 0
            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.c <- dbinom(w.cand[i], 1, psi, log=TRUE)
            for(k in 1:K) {
                N.cand[,k] <- 0
                cells.cand <- cells
                cells.cand[i,k] <- cellFromXY(raster, matrix(u[i,,k], 1))
                if(w.cand[i]==0)
                    cells.cand[i,k] <- cells.cand[i,k] * -1
                counts <- table(cells.cand[,k])
                counts.in <- counts[names(counts) %in% rownames(N)]
                N.cand[names(counts.in),k] <- counts.in
            }
            if(any(N.cand < y))
                next
            lly.c[] <- dbinom(y, N.cand, p, log=TRUE)
            if(runif(1) < exp((sum(lly.c)+prior.c) - (sum(lly)+prior))) {
                w[i] <- w.cand[i]
                N <- N.cand
                lly <- lly.c
                cells[i,] <- cells.cand[i,]
                wUps <- wUps+1
            }
        }


        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update s
        s.ups <- 0
        for(i in 1:M) {
            s.cand <- c(rnorm(1, s[i,1], tune[3]),
                        rnorm(1, s[i,2], tune[3]))
            inbox <- s.cand[1]>=xlims[1] & s.cand[1]<=xlims[2] &
                     s.cand[2]>=ylims[1] & s.cand[2]<=ylims[2]
            if(!inbox)
                next

            llu1.c[i,] <- dnorm(u[i,1,], s.cand[1], tau, log=TRUE)
            llu2.c[i,] <- dnorm(u[i,2,], s.cand[2], tau, log=TRUE)
            llu.c[i,] <- llu1.c[i,] + llu2.c[i,]

            if(runif(1) < exp(sum(llu.c[i,]) - sum(llu[i,]))) {
                s[i,] <- s.cand
                llu1[i,] <- llu1.c[i,]
                llu2[i,] <- llu2.c[i,]
                llu[i,] <- llu.c[i,]
                s.ups <- s.ups+1
            }
        }

        points(s, cex=0.5, col=4)

        # update u
        u.ups <- 0
        u.cand <- u #array(NA, c(M, 2, K))
        N.cand <- N
#        cells.cand <- cells
        for(i in 1:M) {
            for(k in 1:K) {
#                u.cand[i,,k] <- c(rnorm(1, s[i,1], tune[4]),
#                                  rnorm(1, s[i,2], tune[4]))
                u.cand[i,,k] <- c(rnorm(1, s[i,1], tau),
                                  rnorm(1, s[i,2], tau))
#                cat(" u curr =", u[i,,k], "\n")
#                cat("  u cand =", u.cand[i,,k], "\n")

                if(w[i]==0) {
                    u[i,,k] <- u.cand[i,,k]
                    llu1[i,k] <- dnorm(u.cand[i,1,k], s[i,1], tau,
                                       log=TRUE)
                    llu1[i,k] <- dnorm(u.cand[i,2,k], s[i,2], tau,
                                       log=TRUE)
                    llu[i,k] <- llu1[i,k] + llu2[i,k]
#                    cat("    accepted\n")
                    points(u[i,1,k], u[i,2,k], cex=0.5, col=3)
                    next
                } else {
                    # These are the priors
                    # Are they need since u was drawn from prior?
                    llu1.c[i,k] <- dnorm(u.cand[i,1,k], s[i,1], tau,
                                         log=TRUE)
                    llu1.c[i,k] <- dnorm(u.cand[i,2,k], s[i,2], tau,
                                         log=TRUE)
                    llu.c[i,k] <- llu1.c[i,k] + llu2.c[i,k]

                    N.cand[,k] <- 0
#                    inout <- w
#                    inout[w==0] <- -1
#                    cells <- cellFromXY(raster, u[,,k] * inout)
                    cells.cand <- cells
                    cells.cand[i,k] <- cellFromXY(raster,
                                                  matrix(u.cand[i,,k], 1))
                    counts <- table(cells.cand[,k])
                    counts.in <- counts[names(counts) %in% rownames(N)]
                    N.cand[names(counts.in),k] <- counts.in
                    if(any(N.cand[,k] < y[,k]))
                        next

                    lly.c[,k] <- dbinom(y[,k], N.cand[,k], p, log=TRUE)


                    if(runif(1) < exp((sum(lly.c[,k]) + llu.c[i,k]) -
                                      (sum(lly[,k]) + llu[i,k]))) {
                        u[i,,k] <- u.cand[i,,k]
#                        cat("    accepted\n")
                        u.ups <- u.ups+1
                        if(any(N.cand[,k] < y[,k]))
                            stop("N.cand[,k] < n[,k]")
                        N[,k] <- N.cand[,k]
                        llu1[,k] <- llu1.c[,k]
                        llu2[i,k] <- llu2.c[i,k]
                        llu[i,k] <- llu.c[i,k]
                        lly[,k] <- lly.c[,k]
                        cells[i,k] <- cells.cand[i,k]
                points(u[i,1,k], u[i,2,k], cex=0.5, col=3)
                    }
                }
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     s =", s.ups/M, "\n")
            cat("     u =", u.ups/(M*K), "\n")
        }

        out[iter,] <- c(tau, p, psi, sum(w) )


    }
    return(out)
}


