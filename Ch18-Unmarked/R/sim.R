
tr <- seq(0.3, 0.7, length=5)
X <- cbind(rep(tr, each=length(tr)),
           rep(tr, times=length(tr)))    # trap coords
set.seed(4589)
xlim <- c(0, 1); ylim <- c(0, 1)         # S is the unit square
A <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1]) # area of S
mu <- 50                                 # density (animals/unit area)
N <- rpois(1, mu*A)                      # Generate N as Poisson deviate
s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))

