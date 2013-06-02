\name{RSFsim}
\alias{RSFsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  simulate RSF+SCR situation
  %%  ~~function to do ... ~~
}
\description{
simulates RSF + SCR fits models et....
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RSFsim()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

    library(scrbook)
    set.seed(1234)
    gr <- as.matrix(expand.grid(1:40, 1:40))
    Dmat <- as.matrix(dist(gr))
    V <- exp(-Dmat/5)
    z <- t(chol(V)) \%*\% rnorm(1600)
    spatial.plot(gr, z)
    alpha0 <- -2
    sigma <- 2
    beta <- 1
    Ntel <- 4
    nsim <- 10
    Nfixes <- 20
    N <- 100
    simout0 <- simout1 <- simout2 <- simout3 <- matrix(NA, nrow = nsim, 
        ncol = 6)
    for (sim in 1:nsim) {
        cat("sim: ", sim, fill = TRUE)
        Sid <- sample(1:1600, N, replace = TRUE)
        S <- gr[Sid, ]
        poss.tel <- S[, 1] > 5 & S[, 1] < 35 & S[, 2] > 5 & S[, 
            2] < 35
        tel.guys <- sample(Sid[poss.tel], Ntel)
        sid <- tel.guys
        stel <- gr[sid, ]
        if (1 == 2) {
            stest <- source("stest.R")$value
            Ntel <- 8
            Nfixes <- 120
            sid <- stest$id
            stel <- stest$s
        }
        n <- matrix(NA, nrow = Ntel, ncol = 1600)
        i <- 1
        par(mfrow = c(3, 3))
        lammat <- matrix(NA, nrow = Ntel, ncol = 1600)
        for (i in 1:Ntel) {
            d <- Dmat[sid[i], ]
            lam <- exp(1 - (1/(2 * sigma * sigma)) * d * d + 
                beta * z)
            n[i, ] <- rmultinom(1, Nfixes, lam/sum(lam))
            par(mar = c(3, 3, 3, 6))
            lammat[i, ] <- lam
            img <- matrix(lam, nrow = 40, ncol = 40, byrow = FALSE)
            image(1:40, 1:40, rot(img), col = terrain.colors(10))
        }
        X <- cbind(sort(rep(seq(5, 35, 5), 7)), rep(seq(5, 35, 
            5), 7))
        ntraps <- nrow(X)
        raster.point <- rep(NA, nrow(X))
        for (j in 1:nrow(X)) {
            raster.point[j] <- (1:1600)[(X[j, 1] == gr[, 1]) & 
                (X[j, 2] == gr[, 2])]
        }
        points(X, pch = 20, cex = 2)
        D <- e2dist(S, X)
        Zmat <- matrix(z[raster.point], nrow = N, ncol = ntraps, 
            byrow = TRUE)
        loglam <- alpha0 - (1/(2 * sigma * sigma)) * D * D + 
            beta * Zmat
        p <- 1 - exp(-exp(loglam))
        K <- 10
        y <- matrix(NA, nrow = N, ncol = ntraps)
        for (i in 1:N) {
            y[i, ] <- rbinom(ntraps, K, p[i, ])
        }
        cap <- apply(y, 1, sum) > 0
        y <- y[cap, ]
        sbar <- (n \%*\% gr)/as.vector(n \%*\% rep(1, nrow(gr)))
        tmp1 <- optim(c(-3, log(3), 1, 0, 1), intlik3rsf, y = y, 
            K = K, X = X, ztrap = z[raster.point], G = gr)
        simout0[sim, ] <- c(nrow(y), tmp1$estimate)
        tmp2 <- optim(c(-3, log(3), 1, 0, 1), intlik3rsf, y = y, 
            K = K, X = X, ztrap = z[raster.point], G = gr, ntel = n, 
            zall = as.vector(z))
        simout1[sim, ] <- c(nrow(y), tmp2$estimate)
        tmp3 <- optim(c(-3, log(3), 1, 0, 1), intlik3rsf, y = y, 
            K = K, X = X, ztrap = z[raster.point], G = gr, ntel = n, 
            zall = as.vector(z), stel = sbar)
        simout2[sim, ] <- c(nrow(y), tmp3$estimate)
        tmp4 <- optim(c(-3, log(3), 1, 0, 1), intlik3rsf, y = NULL, 
            K = K, X = X, ztrap = z[raster.point], G = gr, ntel = n, 
            zall = as.vector(z))
        simout3[sim, ] <- c(NA, tmp4$estimate)
    }
  
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
