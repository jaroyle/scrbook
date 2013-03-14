

den <- function(N=100,
                Sx=c(0, 100), Sy=c(0, 100), # state-space dims
                Bx=c(30, 70), By=c(30, 70), # inner box dims
                T=10, tau=5, plot=FALSE) {
    AreaS <- (Sx[2]-Sx[1])*(Sy[2]-Sy[1])
    AreaB <- (Bx[2]-Bx[1])*(Bx[2]-Bx[1])
    D <- N/AreaS
    s <- cbind(runif(N, Sx[1], Sx[2]),
               runif(N, Sy[1], Sy[2]))
    if(plot) {
        op <- par(mai=c(0, 0, 0, 0))
        plot(0, type="n", xlim=Sx, ylim=Sy, axes=FALSE, frame=FALSE)
        rect(Sx[1], Sy[1], Sx[2], Sy[2])
        rect(Bx[1], By[1], Bx[2], By[2])
        points(s, pch=16)
    }
    DB <- rep(NA, T)
    for(t in 1:T) {
        u <- cbind(rnorm(N, s[,1], tau),
                   rnorm(N, s[,2], tau))
        if(plot)
            points(u, col=t)
        NB <- sum((u[,1] > Bx[1]) & (u[,1] < Bx[2]) &
                  (u[,2] > By[1]) & (u[,2] < By[2]))
        DB[t] <- NB / AreaB
    }
    if(plot)
        par(op)
    return(list(D=D, DB=DB, DBbar=mean(DB)))
}



den(N=100, T=100)


nsim <- 1000
simout1 <- list()
N1 <- 10
for(i in 1:nsim) {
    simout1[[i]] <- den(N=N1, T=10)
}

N1/1e4
mean(sapply(simout1, function(x) x$DBbar))


# Density in the box, averaged over T time periods for each simulation
hist(sapply(simout1, "[[", "DBbar"))
abline(v=N1/1e4, lwd=2)
