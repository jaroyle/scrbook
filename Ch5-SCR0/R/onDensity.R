
# Simulate N activity centers in S
# and simulate T movement outcomes
# and compute instantanceous density the inner box B
den <- function(N=100,
                Sx=c(0, 100), Sy=c(0, 100), # state-space dims
                Bx=c(30, 70), By=c(30, 70), # inner box dims
                T=10, sigma=5, plot=FALSE) {
    AreaS <- (Sx[2]-Sx[1])*(Sy[2]-Sy[1])
    AreaB <- (Bx[2]-Bx[1])*(Bx[2]-Bx[1])
    D <- N/AreaS
    s <- cbind(runif(N, Sx[1], Sx[2]), # activity centers
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
        u <- cbind(rnorm(N, s[,1], sigma), # movement outcomes
                   rnorm(N, s[,2], sigma))
        if(plot)
            points(u, col=t)
        NB <- sum((u[,1] > Bx[1]) & (u[,1] < Bx[2]) &
                  (u[,2] > By[1]) & (u[,2] < By[2]))
        DB[t] <- NB / AreaB # Instantaneous density in box
    }
    if(plot)
        par(op)
    return(list(D=D, DB=DB, DBbar=mean(DB)))
}



den(N=20, T=100, plot=TRUE)


nsim <- 1000
simout1 <- list()
N1 <- 100
for(i in 1:nsim) {
    simout1[[i]] <- den(N=N1, T=1)
}

N1/1e4
mean(sapply(simout1, function(x) x$DBbar))


# Density in the box, averaged over T time periods for each simulation
hist(sapply(simout1, "[[", "DBbar"))
abline(v=N1/1e4, lwd=2)


# Figures

png("../figs/DvDB.png", width=6, height=6, units="in", res=400)
den(N=100, T=10, plot=TRUE)
dev.off()
system("open ../figs/DvDB.png")



png("../figs/DvDBsim.png", width=6, height=6, units="in", res=400)
hist(sapply(simout1, "[[", "DBbar"), main="",
     xlab="Density in inner box")
abline(v=N1/1e4, lwd=2, col=4)
dev.off()
system("open ../figs/DvDBsim.png")
