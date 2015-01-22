\name{SCRdesign}
\alias{SCRdesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes an SCR design.%%  ~~function to do ... ~~
}
\description{
Computes an SCR design given the state-space and number of traps to
deploy. Note: the function is experimental. The math has not been
checked and the theory is tentative. Progress will be made in 2015. 
}
\usage{
SCRdesign(S = S, C = C, ntraps = 9, ndesigns = 10, nn = 19, sigma = 2, crit = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
    The state-space of the SCR model, i.e., "study area" within which
    traps will be deployed, including a suitable habitat buffer where
    individuals exposed to capture might live.  S should be an M x 2
    matrix of coordinates.
    %%     ~~Describe \code{S} here~~
}
  \item{C}{
    The potential (candidate) trap locations. Also called the design space.
    Normally this is a strict subset of S. 
}
  \item{ntraps}{
 Number of traps to be deployed.
}
  \item{ndesigns}{
 Number of random designs to optimize over. More is better in terms of
 finding better designs but increases the run time.  %%     ~~Describe \code{ndesigns} here~~
}
  \item{nn}{
 Number of nearest neighbors to use in the optimization over the design
 space. %%     ~~Describe \code{nn} here~~
}
  \item{sigma}{
    The function uses a "half-normal" hazard rate model of the form:
    Pr(encounter) = 1-exp(-lambda0*g(||x-s||)) where
    g(||x-s||) is the kernel of a bivariate normal density. sigma is the
    parameter. 
}
  \item{crit}{
    Which criterion to use? In the SCR book there are 4: Q1, Q2, Q3 and
    Q4. These are numbered 4-7 in this function because some partial
    things are used in the debugging. The definition of things is as follows:

## crit=4  = trace(V(alpha))  (crit 1 from the book)
## crit 5 = Var(N)            (crit 2 from the book)
## crit 6 = 1-pbar            (crit 3 from the book)
## crit 7 = var(pbar)         (crit 4 from the book)


}
}
\details{
  The function uses a coordinate exchange algorithm to sequentially
  improve the objective function by swapping one point at a time out of
  the design. 
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with components:

   Qvec = the ordered criteria values of the design found for each
   starting design

   Xid = ID of each design point in the best design over all

   Xlst a list of all designs found (one for each starting design)

   C

   S

   Qhistory
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book chapter 10.

Royle, J. A., & Nychka, D. (1998). An algorithm for the construction of spatial coverage designs with implementation in SPLUS. Computers & Geosciences, 24(5), 479-488.

Nychka, D., Yang, Q., & Royle, J. A. (1997). Constructing spatial designs using regression subset selection. Statistics for the Environment, 3, 131-154.
}
\author{
Andy Royle, aroyle@usgs.gov
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Define a state-space S and a design space C
lower=9
upper=21

delta<- .2 ##1/3   # spacing of design points
pts<-seq(lower,upper,delta)
s1<- sort(rep(pts,length(pts)))
s2<-rep(pts,length(pts))
S<-cbind(s1,s2)

inX<- (s1 <=20 & s1 >=10) & (s2 <=20 & s2>=10)
C<-S[inX,]  # candidate set

dbar<-function(X){
 a<-as.matrix(dis(X))
 diag(a)<-10^10
  mean(apply(a,1,min))
}

library(scrbook)
## crit=4  = trace(V(alpha))
## crit 5 = Var(N) 
## crit 6 = 1-pbar
## crit 7 = var(pbar)
a1.1<-SCRdesign(S,C,ntraps=11,ndesigns=10,nn=15,sigma=1,crit=4)
a2.1<-SCRdesign(S,C,ntraps=11,ndesigns=10,nn=15,sigma=1,crit=5)

# Probably use crit=6 right now
a3.1<-SCRdesign(S,C,ntraps=11,ndesigns=10,nn=15,sigma=1,crit=6)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
