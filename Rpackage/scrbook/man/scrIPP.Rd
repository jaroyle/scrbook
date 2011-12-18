\name{scrIPP}
\alias{scrIPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scrIPP(Z, X, M, niters, xlims, ylims, tune = rep(0.1, 4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{niters}{
%%     ~~Describe \code{niters} here~~
}
  \item{xlims}{
%%     ~~Describe \code{xlims} here~~
}
  \item{ylims}{
%%     ~~Describe \code{ylims} here~~
}
  \item{tune}{
%%     ~~Describe \code{tune} here~~
}
}
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
function(Z, X, M, niters, xlims, ylims, tune=rep(0.1, 4))
{
    Zdims <- dim(Z)
    R <- Zdims[2]
    T <- Zdims[3]

    # initial values
    S <- cbind(runif(M,xlims[1],xlims[2]),runif(M,ylims[1],ylims[2]))
    D <- e2dist1(S, X)
    sigma <-runif(1, .3, .6)
    lam0 <- runif(1, 4, 6)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))

    psi <- runif(1, 0.4, 0.6)
    beta1 <- rnorm(1, 0)

    w <- rbinom(M, 1, psi)
    w[rowSums(Z)>0] <- 1

    # matrix to hold samples
    out <- matrix(NA, nrow=niters, ncol=5)
    colnames(out) <- c("sigma", "lam0", "psi", "beta1", "N")

    cat("\ninitial values =", c(sigma, lam0, psi, beta1, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter \%\% 100 ==0) {
            cat("iter", iter, format(Sys.time(), "\%H:\%M:\%S"), "\n")
            cat("current =", out[iter-1,], "\n")
            cat("  Acceptance rates\n")
            cat("    S =", Sups/M, "\n")
            cat("    w =", wUps/M, "\n")
        }

        ll<- sum(dpois(Z, lam*w, log=TRUE))

        # update sigma
        sigma.cand <- rnorm(1, sigma, tune[1])
        if(sigma.cand > 0) {
            lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE) )
            if(runif(1)<exp( llcand  - ll ) ){
                ll <- llcand
                lam <- lam.cand
                sigma <- sigma.cand
            }
        }

        # update lam0
        lam0.cand <- rnorm(1, lam0, tune[2])
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
            llcand<- sum(dpois(Z, lam.cand*w, log=TRUE) )
            if(runif(1) < exp( llcand - ll ) ) {
                lam0<-lam0.cand
                lam<-lam.cand
                ll <- llcand
            }
        }

        # update w
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i])
                next
            wcand<-w
            if(w[i]==0) {
                wcand[i] <- 1
                ll.w <- 0
                ll.w.cand <- sum(dpois(Z[i,,], lam[i,]*wcand[i], log=TRUE))
            } else {
                wcand[i] <- 0
                ll.w <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                ll.w.cand <- 0
            }
            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
            if(runif(1) < exp((ll.w.cand+prior.cand) - (ll.w+prior))) {
                w <- wcand
                wUps <- wUps+1
            }
        }

        # update psi
        psi <- rbeta(1, 1+sum(w), 1+M-sum(w))

        # update beta1
        D1 <- int2d(beta1, delta=.05)
        beta1.cand <- rnorm(1, beta1, tune[3])
        D1.cand <- int2d(beta1.cand, delta=0.05)
        ll.beta1 <- sum(  beta1*cov(S[,1],S[,2]) - log(D1) )
        ll.beta1.cand <- sum( beta1.cand*(S[,1]+S[,2]) - log(D1.cand) )
        if(runif(1) < exp(ll.beta1.cand - ll.beta1) )  {
          beta1<-beta1.cand
          }


        # update S
        Sups <- 0
        for(i in 1:M) {
            Scand <- c(rnorm(1, S[i,1], tune[4]),
                       rnorm(1, S[i,2], tune[4]))
            inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
                     Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
            if(!inbox)
                next
            dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
            lam.cand <- lam
            lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
            if(w[i]==0)
                ll.S <- ll.S.cand <- 0
            else {
                ll.S <- sum(dpois(Z[i,,], lam[i,], log=TRUE) )
                ll.S.cand <- sum(dpois(Z[i,,], lam.cand[i,], log=TRUE) )
            }
            #ln(prior), denominator is constant
            prior.S <- beta1*cov(S[i,1], S[i,2]) # - log(D1)
            prior.S.cand <- beta1*(Scand[1] + Scand[2]) # - log(D1)

           if(runif(1)< exp((ll.S.cand+prior.S.cand) - (ll.S+prior.S))) {
                S[i,] <- Scand
                lam <- lam.cand
                D[i,] <- dtmp
                ##psi[i] <- psi.cand
                Sups <- Sups+1
            }
        }


        out[iter,] <- c(sigma,lam0,psi,beta1,sum(w) )
    }
    last <- list(S=S, lam=lam, w=w)
    list(out=out, last=last)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
